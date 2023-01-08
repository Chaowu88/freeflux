#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '02/16/2022'




import re
from functools import lru_cache, reduce
from collections import OrderedDict, deque
from collections.abc import Iterable
from itertools import chain
import numpy as np
import pandas as pd
from sympy import Integer
from multiprocessing import Pool
from .metabolite import Metabolite
from .reaction import Reaction
from .emu import EMU
from ..io.inputs import read_model_from_file
from ..analysis.simulate import Simulator
from ..analysis.inst_simulate import InstSimulator
from ..analysis.fit import Fitter
from ..analysis.inst_fit import InstFitter
from ..optim.optim import Optimizer




class Model():
    '''
    Model class is the fundamental class of FreeFlux hosting operations for 13C MFA.

    Model can be built by adding reactions one by one or reading set of reactions from .tsv or .xlsx file.
    
    Parameters
    ----------
    name: str
        Model name.

    Attributes
    ----------
    name: str
        Model name.
    metabolites_info: dict
        Metabolite ID => list of Metabolites.
    reactions_info: OrderedDict
        Reaction ID => Reaction.
    metabolites: list
        Metabolite IDs, in alphabetical order.
    metabolites_with_atoms: list
        IDs of metabolite with atom assignment, in alphabetical order. 
    end_substrates: list
        Initial substrates of the model, in alphabetical order. 
    end_products: list
        Final products of the model, in alphabetical order. 
    reactions: list
        Reaction IDs, in order of addition.
    n_metabolites: int
        # of metabolites.
    n_reactions: int
        # of reactions.
    _full_net_stoichiometric_matrix: df
        Complete stoichiometric matrix for net reaction with all metabolites in rows, 
        net reactions in columns.
    _full_total_stoichiometric_matrix: df
        Complete stoichiometric matrix for total reaction with all metabolites in rows, 
        total reactions in columns.
    metabolite_adjacency_matrix: df
        Metabolite adjacency matrix (MAM). Metabolites with atoms are in index and 
        columns (no duplicates). List of Reactions are in cells if reactions exists 
        between (sub, pro), [] otherwise.    
    net_fluxes_bounds: dict
        Reaction ID => [lb, ub] by setting.   
    net_fluxes_range: dict
        Reaction ID => estimeted [lb, ub]. All required net fluxes are included.  
    netfluxids: list
        Net flux IDs, alias of reactions. 
    concentrations: ser
        Concentrations, metabolite ID => float.    
    concentrations_bounds: dict
        Metabolite ID => [lb, ub] by setting.    
    concentrations_range: dict
        Metabolite ID => [lb, ub]. All required concentrations are included.    
    concids: list
        Concentration IDs (used only in computation process).    
    total_fluxes: ser
        Total fluxes, fluxe ID (e.g., 'v1_f' or 'v2') => float.   
    totalfluxids: list
        Total flux IDs.    
    target_EMUs: list
        Target EMU IDs.   
    timepoints: list
        Sorted time points for MDV simulation.   
    substrate_MDVs: dict
        Substrate EMU => MDV.   
    substrate_MDVs_der_p: dict
        Substrate EMU => derivatives of substrate MDV w.r.t. variables 
        in shape of (len(MDV), # of vars).
        # of vars = # of free fluxes for steady state MFA;
        # of vars = # of free fluxes + # of concentrations for INST MFA.
    initial_matrix_Xs: dict
        Size => initial MDVs of EMU in matrix X, i.e., the natural MDVs.
        The initial MDV matrix has the same shape of matrix X.   
    initial_matrix_Ys: dict
        Size => initial MDVs of EMU in matrix Y, i.e., either natural MDVs or labeled MDVs.
        The initial MDV matrix has the same shape of matrix Y.        
    initial_matrix_Xs_der_p, initial_matrix_Ys_der_p: dict
        Size => 3-D array in shape of (# of vars), X(Y).shape[0], X(Y).shape[1])
        which is the initial MDV derivatives of EMUs in matrix X(Y) w.r.t. variables. 
        # of vars = # of free fluxes + # of concentrations for INST MFA.  
    initial_sim_MDVs: dict
        EMU ID => {t0 => MDV}. MDV of target EMUs at t0.
    EAMs: dict of df
        Size => EMU adjacency matrix (EAM). Cells are symbolic expression of flux.    
    matrix_As, matrix_Bs: dict
        Size => [lambdified matrix A(B), [flux IDs], [EMUs]].     
    matrix_As_der_p, matrix_Bs_der_p: dict
        Size => 3-D array in shape of (# of vars, A(B).shape[0], A(B).shape[1])
        which is the derivatives of matrix A(B) w.r.t. variables.
        # of vars = # of free fluxes for steady state MFA;
        # of vars = # of free fluxes + # of concentrations for INST MFA.    
    matrix_Ms: dict
        Size => [lambdified matrix M, [metabolite IDs]].    
    matrix_Ms_der_p: dict
        Size => 3-D array in shape of (# of vars, M.shape[0], M.shape[1]),
        which is the derivatives of matrix M w.r.t. variables.
        # of vars = # of free fluxes for steady state MFA;
        # of vars = # of free fluxes + # of concentrations for INST MFA.    
    labeling_strategy: dict
        Metabolite ID => [labeling_pattern(s), percentage(s), purity(s)].    
    measured_MDVs: dict
        EMU ID (metabolite ID + '_' + atom NOs) => [means of MDV, SDs of MDV].    
    measured_MDVs_inv_cov: array
        Inversed covariance matrix of measured MDVs with variances on the diagnol, 
        other elements are zero.    
    measured_fluxes: dict
        Flux ID (e.g., 'v1_f' or 'v2') => [mean, SD].    
    measured_fluxes_inv_cov: array
        Inversed covariance matrix of measured fluxes with variances on the diagnol, 
        other elements are zero.    
    measured_fluxes_der_p: array
        Derivative of measured fluxes w.r.t. variables in shape of 
        (# of measured fluxes, # of vars),
        # of vars = # of free fluxes for steady state MFA;
        # of vars = # of free fluxes + # of concentrations for INST MFA.
    measured_inst_MDVs: dict
        EMU ID (metabolite ID + '_' + atom NOs) => {timepoint => [means of MDV, SDs of MDV]}.    
    measured_inst_MDVs_inv_cov: array
        Inversed covariance matrix of measured concatenated MDVs with variances on the diagnol, 
        other elements are zero. Timepoints are concatenated except t0.     
    null_space: 2-D array
        Null space of total stoichiometric matrix.    
    transform_matrix: 2-D array
        Transform matrix letting vnet = transform_matrix*v.
    '''
    
    def __init__(self, name = 'unnamed'):
        '''
        Parameters
        ----------
        name: str
            Model name.
        '''
        
        self.name = name
        self.reactions_info = OrderedDict()
        
        self.target_EMUs = []
        self.timepoints = []
        self.substrate_MDVs = {}
        self.substrate_MDVs_der_p = {}
        self.EAMs = {}
        self.matrix_As = {}
        self.matrix_Bs = {}
        self.matrix_Ms = {}
        self.matrix_As_der_p = {}
        self.matrix_Bs_der_p = {}
        self.matrix_Ms_der_p = {}
        
        self.null_space = None
        self.transform_matrix = None
        self.initial_matrix_Xs = {}
        self.initial_matrix_Ys = {}
        self.initial_matrix_Xs_der_p = {}
        self.initial_matrix_Ys_der_p = {}
        self.initial_sim_MDVs = {}
        
        self.labeling_strategy = {}
        self.measured_MDVs = {}
        self.measured_fluxes = {}
        self.measured_inst_MDVs = {}
        self.measured_MDVs_inv_cov = None
        self.measured_inst_MDVs_inv_cov = None
        self.measured_fluxes_inv_cov = None
        self.measured_fluxes_der_p = None
        self.unbalanced_metabolites = set()
        
        self.net_fluxes_bounds = {}
        self.net_fluxes_range = OrderedDict()
        self.concentrations_bounds = {}
        self.concentrations_range = OrderedDict()
        self.total_fluxes = pd.Series(dtype = float)
        self.concentrations = pd.Series(dtype = float)
     

    def add_reactions(self, reactions):
        '''
        Parameters
        ----------
        reactions: Reaction or list of Reaction
        '''
        
        if not isinstance(reactions, list):
            reactions = [reactions]
        
        new_reactions_info = OrderedDict([(rxn.id, rxn) for rxn in reactions])
        self.reactions_info.update(new_reactions_info)
        
        for rxn in reactions:
            if rxn.host_models is None:
                rxn.host_models = set([self])
            else:
                rxn.host_models.add(self)
    

    def remove_reactions(self, reactions):
        '''
        Parameters
        ----------
        reactions: Reaction or list of Reaction
        '''
        
        if not isinstance(reactions, list):
            reactions = [reactions]
            
        for rxn in reactions:
            self.reactions_info.pop(rxn.id)
            
        for rxn in reactions:
            rxn.host_models.remove(self)
            if not rxn.host_models:
                rxn.host_models = None
                break    
        
        
    def read_from_file(self, file):
        '''
        Parameters
        ----------
        file: file path
            tsv or excel file with reactions with fields "reaction_ID", "substrate_IDs(atom)", 
            "product_IDs(atom)" and "reversibility".

            Header line starts with "#", and will be skiped.
        '''
        
        dataRaw = read_model_from_file(file)
        
        data = pd.DataFrame()
        for col, ser in dataRaw.iteritems():
            if ser.dtype == object:
                data[col] = ser.str.replace(r'\s+', '')
            else:
                data[col] = ser
        
        for rxn, (subsStr, prosStr, rev) in data.iterrows():
            v = Reaction(rxn, reversible = bool(int(rev)))
            
            for stoy, metab, atoms in re.findall(r'([0-9\.]+|)([.\w]+)\(?([a-z0-9\.,]+|)\)?', subsStr):
                stoy = float(stoy) if stoy else 1.0
                atoms = atoms.split(',') if atoms else None
                v.add_substrates(Metabolite(metab, atoms), stoy)
    
            for stoy, metab, atoms in re.findall(r'([0-9\.]+|)([.\w]+)\(?([a-z0-9\.,]+|)\)?', prosStr):
                stoy = float(stoy) if stoy else 1.0
                atoms = atoms.split(',') if atoms else None
                v.add_products(Metabolite(metab, atoms), stoy)
            
            self.add_reactions(v)
    
    
    @property
    def metabolites_info(self):
        '''
        Returns
        -------
        metabsInfo: dict
            metabolite ID => list of Metabolites
        '''
        
        metabsInfo = {}
        for _, rxn in self.reactions_info.items():
            for subid, sub in rxn.substrates_info['metab'].iteritems():
                metabsInfo.setdefault(subid, []).append(sub)
            for proid, pro in rxn.products_info['metab'].iteritems():
                metabsInfo.setdefault(proid, []).append(pro)
        
        metabsInfo = {metabid: list(set(metabs)) for metabid, metabs in metabsInfo.items()}
        
        return metabsInfo
        
    
    @property
    def metabolites(self):
        
        metabids = []
        for _, rxn in self.reactions_info.items():
            metabids = metabids + rxn.substrates + rxn.products
            
        return sorted(set(metabids))
    
    
    @property
    def metabolites_with_atoms(self):
        
        metabids = []
        for _, rxn in self.reactions_info.items():
            metabids = metabids + rxn.substrates_with_atoms + rxn.products_with_atoms
            
        return sorted(set(metabids))
    
    
    @property
    def reactions(self):
        
        return list(self.reactions_info.keys())
    
    
    @property
    def n_metabolites(self):
    
        return len(self.metabolites)
        
        
    @property
    def n_reactions(self):
    
        return len(self.reactions)
        
        
    @property
    def _full_net_stoichiometric_matrix(self):
        '''
        Returns
        -------
        netS: df
            Complete stoichiometric matrix for net reaction with all metabolites in rows and 
            net reactions in columns.
        '''
        
        netS = pd.DataFrame(0, index = self.metabolites, columns = self.reactions)
        for rxnid, rxn in self.reactions_info.items():
            
            for sub in rxn.substrates:
                value = rxn.substrates_info.loc[sub, 'stoy']
                if isinstance(value, pd.Series):
                    netS.loc[sub, rxnid] = -value.sum()   
                else:
                    netS.loc[sub, rxnid] = -value   
            
            for pro in rxn.products:
                value = rxn.products_info.loc[pro, 'stoy']
                if isinstance(value, pd.Series):
                    netS.loc[pro, rxnid] = value.sum()   
                else:
                    netS.loc[pro, rxnid] = value   
        
        return netS
    
    
    @property
    def _full_total_stoichiometric_matrix(self):
        '''
        Returns
        -------
        totalS: df
            Complete stoichiometric matrix for total reaction with all metabolites in rows and
            total reactions in columns.
        '''
        
        netS = self._full_net_stoichiometric_matrix
        
        totalS = pd.DataFrame(index = netS.index)
        for rxn, col in netS.iteritems():
            if self.reactions_info[rxn].reversible:
                totalS[rxn+'_f'] = col
                totalS[rxn+'_b'] = -col
            else:
                totalS[rxn] = col
        
        return totalS


    @property
    def end_substrates(self):
    
        totalS = self._full_total_stoichiometric_matrix
        
        nLessthan0 = totalS[totalS < 0].count(axis = 1)
        nGreaterthan0 = totalS[totalS > 0].count(axis = 1)
        
        return sorted(totalS.index[(nLessthan0 > 0) & (nGreaterthan0 == 0)].tolist())
        
        
    @property
    def end_products(self):
    
        totalS = self._full_total_stoichiometric_matrix
        
        nLessthan0 = totalS[totalS < 0].count(axis = 1)
        nGreaterthan0 = totalS[totalS > 0].count(axis = 1)
        
        return sorted(totalS.index[(nGreaterthan0 > 0) & (nLessthan0 == 0)].tolist())
        

    def get_net_stoichiometric_matrix(self, exclude_metabs = None, include_ends = False):
        '''
        Parameters
        ----------
        exclude_metabs: list or set
            Excluded metabolites.
        include_ends: bool
            Whether to include end metabolites (i.e., initial substrates and final products).
        
        Returns
        netS: df
            Net stoichiometric matrix with balanced metabolites in rows and net reactions in columns
        '''
        
        if exclude_metabs is None:
            exclude_metabs = []
        else:
            exclude_metabs = list(exclude_metabs)
            
        netS = self._full_net_stoichiometric_matrix
        
        if include_ends:
            return netS.loc[netS.index.difference(exclude_metabs), :]
        else:
            return netS.loc[netS.index.difference(exclude_metabs + self.end_substrates + self.end_products), :]
        
        
    def get_total_stoichiometric_matrix(self, exclude_metabs = None, include_ends = False):
        '''
        Parameters
        ----------
        exclude_metabs: list or set
            Excluded metabolites.
        include_ends: bool
            Whether to include end metabolites (initial substrates and final products).
        
        Returns
        totalS: df
            Total stoichiometric matrix with balanced metabolites in rows and total reactions in columns
        '''
        
        if exclude_metabs is None:
            exclude_metabs = []
        else:
            exclude_metabs = list(exclude_metabs)    
        
        totalS = self._full_total_stoichiometric_matrix
        
        if include_ends:
            return totalS.loc[totalS.index.difference(exclude_metabs), :]
        else:
            return totalS.loc[totalS.index.difference(exclude_metabs + self.end_substrates + self.end_products), :]
    
    
    @property
    @lru_cache()
    def netfluxids(self):
        
        return self.reactions
        
        
    @property
    @lru_cache()
    def totalfluxids(self):
        
        totalfluxids = []
        for rxnid, rxn in self.reactions_info.items():
            if rxn.reversible:
                totalfluxids.extend([rxnid+'_f', rxnid+'_b'])
            else:
                totalfluxids.append(rxnid)
        
        return totalfluxids


    @property
    @lru_cache()
    def concids(self):
        
        if self.EAMs:
            metabids = []
            for _, EAM in self.EAMs.items():
                metabids += [emu.metabolite_id for emu in EAM.columns]
            return sorted(set(metabids))
        else:
            raise ValueError('network not decomposed')
        
    
    @property
    @lru_cache()
    def metabolite_adjacency_matrix(self):
        '''
        Returns
        -------
        MAM: df
            Metabolite adjacency matrix (MAM). Metabolites with atoms are in index and columns (no duplicates). 
            List of Reactions are in cells if reactions exists between sub (index) and pro (columns), [] otherwise
        '''
        
        MAM = pd.DataFrame(index = self.metabolites_with_atoms, columns = self.metabolites_with_atoms)
        MAM = MAM.applymap(lambda x: [])
        for _, rxn in self.reactions_info.items():
            for sub in rxn.substrates_with_atoms:
                for pro in rxn.products_with_atoms:
                    if rxn.reversible:
                        MAM.loc[sub, pro].append(rxn)
                        MAM.loc[pro, sub].append(rxn)
                    else:
                        MAM.loc[sub, pro].append(rxn)
        
        return MAM
    
    
    def __repr__(self):
    
        headStr =  '%s %s (%s metabolites, %s reactions)' % (self.__class__.__name__, 
                                                             self.name if self.name else 'unknown', 
                                                             self.n_metabolites, 
                                                             self.n_reactions)
        bodyStr = '\n'.join([str(rxn) for rxn in self.reactions_info.values()])

        return headStr + '\n' + bodyStr 


    def _BFS(self, iniEMU):
        '''
        Parameters
        ----------
        iniEMU: EMU
            Starting EMU of the decomposition. 
            Metabolite of iniEMU can be any Metabolite instance with the same id.
        
        Returns
        -------
        EAMsInfo: dict of list
            Size => list of [EMU, [precursor EMUs], symbolic expression of flux].
        '''
        
        MAM = self.metabolite_adjacency_matrix
        
        EAMsInfo = {}
        
        searched = []
        toSearch = deque()
        toSearch.appendleft(iniEMU)
        while toSearch:
            
            currentEMU = toSearch.pop()
            searched.append(currentEMU)
                    
            formingRxns = list(set(chain(*[cell for cell in MAM[currentEMU.metabolite_id] if cell != []])))
            for formingRxn in formingRxns:
                
                if formingRxn.reversible:
                    if currentEMU.metabolite_id in formingRxn.products_with_atoms:
                        asProMetabs = formingRxn.products_info['metab'][currentEMU.metabolite_id]
                        direction = 'forward'
                        flux = formingRxn.fflux
                    else:
                        asProMetabs = formingRxn.substrates_info['metab'][currentEMU.metabolite_id]
                        direction = 'backward'
                        flux = formingRxn.bflux
                else:
                    asProMetabs = formingRxn.products_info['metab'][currentEMU.metabolite_id]
                    direction = 'forward'
                    flux = formingRxn.flux
                
                if isinstance(asProMetabs, pd.Series):
                    offset = 1 / asProMetabs.size
                    asProMetabs = list(asProMetabs)
                else:
                    offset = 1.0
                    asProMetabs = [asProMetabs]
                
                
                for asProMetab in asProMetabs:
                    currentEMU = EMU(currentEMU.id, asProMetab, currentEMU.atom_nos)
                    preEMUsInfo = formingRxn._find_precursor_EMUs(currentEMU, direction = direction)
                    
                    for preEMUs, coe in preEMUsInfo:
                        for preEMU in preEMUs:
                            if preEMU not in searched and preEMU not in toSearch:
                                toSearch.appendleft(preEMU)
                            
                        EAMsInfo.setdefault(currentEMU.size, []).append([currentEMU, preEMUs, offset * coe * flux])
                            
        return EAMsInfo
    

    def _get_original_EAMs(self, iniEMU):
        '''
        Parameters
        ----------
        iniEMU: EMU
            Starting EMU of the decomposition.
            
        Returns
        -------
        EAMs: dict of df
            Size => original EMU adjacency matrix (EAM), cells are symbolic expression of flux。
        '''
        
        EAMsInfo = self._BFS(iniEMU)
        
        EAMs = {}
        for size, EMUsInfo in EAMsInfo.items():
            
            nonSourceEMUs = set([EMUInfo[0] for EMUInfo in EMUsInfo])
            SourceEMUs = sorted(set([tuple(EMUInfo[1]) if len(EMUInfo[1]) > 1 else EMUInfo[1][0] for EMUInfo in EMUsInfo]) 
                                - nonSourceEMUs)
            
            EAM = pd.DataFrame(Integer(0), index = sorted(nonSourceEMUs) + sorted(SourceEMUs), columns = sorted(nonSourceEMUs))
            for emu, preEMUs, flux in EMUsInfo:
                
                col = emu
                
                if len(preEMUs) == 1:
                    idx = preEMUs[0]
                else:
                    idx = tuple(preEMUs)
                
                EAM.loc[[idx], col] += flux
            
            EAMs[size] = EAM
        
        return EAMs
            
    
    def _replace_list_item(self, iterable, toReplace, value):
        '''
        Parameters
        ----------
        iterable: iterable
            Iterable can be nested.
        toReplace: scalar
            Value to be replaced.
        value: scalar
            Value that replaces toReplace.
            
        Returns
        -------
        newLst: tuple
            New tuple with toReplace recursively replaced.
        '''
        
        newLst = []
        for item in iterable:
            
            if not isinstance(item, Iterable):
                newItem = value if item == toReplace else item
            else:
                newItem = self._replace_list_item(item, toReplace, value)
            
            newLst.append(newItem)
        
        return tuple(newLst)
        
        
    def _uniquify_dataFrame_index(self, df):
        '''
        Parameters
        ----------
        df: df
            DataFrame to be uniquify.
            
        Returns
        -------
        uniqueDf: df
            DataFrame with duplicate rows combined (summated).
        '''
        
        sortedDf = df.sort_index()
    
        sortDfIndex, idx = np.unique(sortedDf.index, return_index = True)
    
        uniqueDf = pd.DataFrame(np.add.reduceat(sortedDf.values, idx), index = sortDfIndex, columns = sortedDf.columns)
        
        return uniqueDf
        
    
    def _lump_linear_EMUs(self, EAMs, iniEMU):
        '''
        Parameters
        ----------
        EAMs: dict of df
            Size => original EMU adjacency matrix (EAM), cells are symbolic expression of flux.
        iniEMU: EMU
            Starting EMU of the decomposition.
        
        Returns
        -------
        lumpedEAMs: dict of df
            Size => lumped EMU adjacency matrix (EAM), cells are symbolic expression of flux.
        '''
        
        lumpedEAMs = {}
        orderedSizes = sorted(EAMs.keys(), reverse = True)
        for i, size in enumerate(orderedSizes):
            
            lumpedEAM = EAMs[size].copy(deep = 'all')
            for emu in lumpedEAM.columns:
                
                preEMUs = lumpedEAM.index[lumpedEAM[emu] != 0]
                if preEMUs.size == 1:
                    preEMU = preEMUs[0]
                    
                    if emu != iniEMU and (preEMU not in lumpedEAM.columns or lumpedEAM.loc[emu, preEMU] == 0):
                        lumpedEAM.drop(emu, axis = 1, inplace = True)
                        lumpedEAM.index = self._replace_list_item(lumpedEAM.index, emu, preEMU)
                        
                        for j in range(i):
                            largerEAM = lumpedEAMs[orderedSizes[j]]
                            largerEAM.index = self._replace_list_item(largerEAM.index, emu, preEMU)
                        
                        lumpedEAM = self._uniquify_dataFrame_index(lumpedEAM)
                        upper = self._uniquify_dataFrame_index(lumpedEAM.loc[lumpedEAM.columns, :])
                        lower = self._uniquify_dataFrame_index(lumpedEAM.loc[lumpedEAM.index.difference(lumpedEAM.columns), :])
                        lumpedEAM = pd.concat((upper, lower))
            
            lumpedEAMs[size] = lumpedEAM
            
        return lumpedEAMs
        
        
    def _combine_equivalent_EMUs(self, EAMs):
        '''
        Parameters
        ----------
        EAMs: dict of df
            Size => original EMU adjacency matrix (EAM), cells are symbolic expression of flux.
            
        Returns
        -------
        combinedEAMs: dict df
            Size => EAM with equivalent EMUs combined, cells are symbolic expression of flux.
        '''
        
        combinedEAMs = {}
        for size, EAM in EAMs.items():
            
            combinedEAM = EAM.copy(deep = 'all')
            combined = []
            for emu in combinedEAM.columns:
                if emu not in combined:
                    
                    equivEMU = emu.equivalent
                    if equivEMU in combinedEAM.columns:   

                        combinedEAM.loc[:, emu] = combinedEAM.loc[:, [emu, equivEMU]].sum(axis = 1) / 2
                        combinedEAM.drop(equivEMU, axis = 1, inplace = True)
                        
                        combinedEAM.loc[emu, :] = combinedEAM.loc[[emu, equivEMU], :].sum()
                        combinedEAM.drop(equivEMU, inplace = True)
                        
                        combined.append(equivEMU)
                    
            combinedEAMs[size] = combinedEAM

        return combinedEAMs
                    
        
    def get_emu_adjacency_matrices(self, iniEMU, lump = True):
        '''
        Parameters
        ----------
        iniEMU: EMU
            Starting EMU of the decomposition.
        lump: bool
            Whether to lump linear EMUs.
        
        Returns
        -------
        EAMs: dict of df
            Size => EMU adjacency matrix (EAM) after lumping of linear EMUs and combination
            of equivalent EMUs. Index and columns are EMUs, cells are symbolic expression of flux.

        Notes
        -----
        EMUs in sequential reactions can not be lumped in transient MFA.
        '''
        
        import platform
        if platform.system() == 'Linux':
            import os
            os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
            
        oriEAMs = self._get_original_EAMs(iniEMU)
        
        if lump:
            lumpedEAMs = self._lump_linear_EMUs(oriEAMs, iniEMU)
            combinedEAMs = self._combine_equivalent_EMUs(lumpedEAMs)
        else:
            combinedEAMs = self._combine_equivalent_EMUs(oriEAMs)
        
        return combinedEAMs
        
    
    def _merge_EAMs(self, EAM1, EAM2):
        '''
        Parameters
        ----------
        EAM1: df
            EMU adjacency matrix (EAM) to merge
        EAM2: df
            EMU adjacency matrix (EAM) to merge
            
        Returns
        -------
        mergedEAM: df
            merged EAM
        '''
        
        nonSourceEMUunion = EAM2.columns.union(EAM1.columns)
        sourceEMUunion = EAM2.index.difference(EAM2.columns).union(EAM1.index.difference(EAM1.columns))
        
        mergedEAM = pd.DataFrame(Integer(0), index = nonSourceEMUunion.append(sourceEMUunion), columns = nonSourceEMUunion)
        mergedEAM.loc[EAM1.index, EAM1.columns] = EAM1
        mergedEAM.loc[EAM2.index, EAM2.columns] = EAM2
        
        return mergedEAM
    
    
    def _merge_all_EAMs(self, *EAMsAll):
        '''
        Parameters
        ----------
        EAMsAll: tuple of EAMs
            EAMs is dict of DataFrame, i.e., Size => EMU adjacency matrix (EAM).

        Returns
        -------    
        mergedEAMs: dict of df
            Size => merged EAM     
        '''
        
        mergedEAMs = {}
        maxsize = max([max(EAMs) for EAMs in EAMsAll])
        for size in range(1, maxsize+1):
            EAMCurrentSize = list(filter(lambda x: isinstance(x, pd.DataFrame), [EAMs.get(size, 0) for EAMs in EAMsAll]))
            if EAMCurrentSize:
                mergedEAMs[size] = reduce(self._merge_EAMs, EAMCurrentSize)
        
        return mergedEAMs            
        
    
    def _decompose_network(self, metabolites, atom_nos, lump = True, n_jobs = 1):
        '''
        Parameters
        ----------
        metabolites: list of str
            List of metabolite IDs from which initial EMU will be generated to start the decomposition.
        atom_nos: list of str
            Atom NOs of corresponding metabolites, len(atom_nos) should be equal to len(metabolites).
        lump: bool
            Whether to lump linear EMUs.    
        n_jobs: int
            # of jobs to run in parallel.
        
        Returns
        -------
        mergedEAMs: dict of df
            Size => merged EMU adjacency matrix (EAM).

        Notes
        -----
        EMUs in sequential reactions can not be lumped in transient MFA.    
        '''
        
        emus = []
        for metabid, atomNOs in zip(metabolites, atom_nos):
            emus.append(EMU(metabid+'_'+atomNOs, Metabolite(metabid), atomNOs))
            
        EAMsAll = []
        if n_jobs == 1:
            for emu in emus:
                EAMs = self.get_emu_adjacency_matrices(emu, lump)
                EAMsAll.append(EAMs)
        else:
            pool = Pool(processes = n_jobs)
            
            for emu in emus:
                EAMs = pool.apply_async(func = self.get_emu_adjacency_matrices, args = (emu, lump))
                EAMsAll.append(EAMs)
            
            pool.close()    
            pool.join()
            
            EAMsAll = [EAMs.get() for EAMs in EAMsAll]
        
        mergedEAMs = self._merge_all_EAMs(*EAMsAll)
        
        return mergedEAMs    
        
        
    def decompose_network(self, ini_emus, lump = True, n_jobs = 1):
        '''
        Parameters
        ----------
        ini_emus: dict
            Metabolite ID => atom NOs or list of atom NOs. Atom NOs can be int list or str, 
            e.g., {'Ala': [[1,2,3], '23'], 'Ser': '123'}
        lump: bool
            Whether to lump linear EMUs.
        n_jobs: int
            # of jobs to run in parallel.
            
        Returns
        -------
        mergedEAMs: dict of df
            Size => merged EMU adjacency matrix (EAM)

        Notes
        -----
        EMUs in sequential reactions can not be lumped in transient MFA.      
        '''
        
        metabids = []
        atom_nos = []
        for metabid, atomNOs in ini_emus.items():
            if isinstance(atomNOs, list):
                if any(isinstance(item, Iterable) for item in atomNOs):
                    for atomnos in atomNOs:
                        if not isinstance(atomnos, str):
                            atomnos = ''.join(map(str, atomnos))
                        metabids.append(metabid)
                        atom_nos.append(atomnos)
                else:
                    atomNOs = ''.join(map(str, atomNOs))
                    metabids.append(metabid)
                    atom_nos.append(atomnos)
            else:
                metabids.append(metabid)
                atom_nos.append(atomNOs)            
        
        return self._decompose_network(metabids, atom_nos, lump = lump, n_jobs = n_jobs)        
    
    
    def optimizer(self):
    
        return Optimizer(self)
        
        
    def simulator(self, kind):
        '''
        Parameters
        ----------
        kind: {"ss", "inst"}
            * If "ss", simulation at isotopic steady state is performed.
            * If "inst", simulation at isotopically nonstationary state is performed.
        '''
        
        if kind.lower() == 'ss':
            return Simulator(self)
        elif kind.lower() == 'inst':
            return InstSimulator(self)
        else:
            raise ValueError('kind should be either "ss" or "inst"')
            
            
    def fitter(self, kind):
        '''
        Parameters
        ----------
        kind: {"ss", "inst"}
            * If "ss", fitting at isotopic steady state is performed.
            * If "inst", fitting at isotopically nonstationary state is performed.
        '''
        
        if kind.lower() == 'ss':
            return Fitter(self)
        elif kind.lower() == 'inst':
            return InstFitter(self)
        else:
            raise ValueError('kind should be either "ss" or "inst"')
