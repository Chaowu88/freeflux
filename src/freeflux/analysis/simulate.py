#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '03/30/2020'




from collections.abc import Iterable
from functools import partial
from ..core.mdv import MDV
from ..io.inputs import read_preset_values_from_file
from ..io.results import SimResults
from ..utils.utils import Calculator
from ..utils.context import Context




class Simulator():
    '''
    Parameters
    ----------
    model: Model
        Freeflux Model.
    '''
    
    def __init__(self, model):
        '''
        Parameters
        ----------
        model: Model
            Freeflux Model.
        '''
        
        self.model = model
        self.calculator = Calculator(self.model)
        self.contexts = []
        
    
    def __enter__(self):
        
        self.contexts.append(Context())

        return self
        
        
    def __exit__(self, type, value, traceback):
        
        context = self.contexts.pop()
        context.undo()
    

    def set_target_EMUs(self, target_emus):
        '''
        Parameters
        ----------
        target_emus: dict
            Metabolite ID => atom NOs or list of atom NOs. Atom NOs can be int list or str, 
            e.g.,, {'Ala': [[1,2,3], [2,3]], 'Ser': '123'}.
        '''
        
        emuids = []
        for metabid, atomNOs in target_emus.items():
            if isinstance(atomNOs, list):
                if any(isinstance(item, Iterable) for item in atomNOs):
                    for atomnos in atomNOs:
                        if not isinstance(atomnos, str):
                            atomnos = ''.join(map(str, atomnos))
                        emuid = metabid+'_'+atomnos
                        emuids.append(emuid)
                        self.model.target_EMUs.append(emuid)
                else:
                    atomNOs = ''.join(map(str, atomNOs))
                    emuid = metabid+'_'+atomNOs
                    emuids.append(emuid)
                    self.model.target_EMUs.append(emuid)
            else:
                emuid = metabid+'_'+atomNOs
                emuids.append(emuid)
                self.model.target_EMUs.append(emuid)
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_target_EMUs, emuids))
        
    
    def _unset_target_EMUs(self, emuids):
        
        for emuid in emuids:
            if emuid in self.model.target_EMUs:
                self.model.target_EMUs.remove(emuid)
    
        
    def set_labeling_strategy(self, labeled_substrate, labeling_pattern, percentage, purity):
        '''
        Use this method for every substrate tracer.
        
        Parameters
        ----------
        labeled_substrate: str
            Metabolite ID.
        labeling_pattern: str or list of str
            Labeling pattern of substrate, '0' for unlabeled atom, '1' for labeled atom, 
            e.g., '100000' for 1-13C glucose. 
            
            List if tracer with multiple labeling patterns are used. 
            
            Natural substrate (with all '0's) don't need to be explicitly set.
            
            If str, labeling_pattern should not be natural substrate.
        percentage: float or list of float
            Molar percentage (in range of [0,1]) of corresponding tracer. 
            Sum of percentage should be <= 1, and the rest will be considered as natural substrate.
            
            List if tracer with multiple labeling patterns are used. 
            
            * If list, len(percentage) should be equal to len(labeling_pattern).
            * If float, labeling_pattern should not be natural substrate.
        purity: float or list of float
            Labeled atom purity (in range of [0,1]) of corresponding tracer.
            
            List if tracer with multiple labeling patterns are used.

            * If list, len(purity) should be equal to len(labeling_pattern).
            * If float, labeling_pattern should not be natural substrate.
        '''
    
        self.model.labeling_strategy[labeled_substrate] = [labeling_pattern, percentage, purity]
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_labeling_strategy, labeled_substrate))
    
        
    def _unset_labeling_strategy(self, labeled_substrate):
        '''
        Parameters
        ----------
        labeled_substrate: str
            Metabolite ID.
        '''
        
        if labeled_substrate in self.model.labeling_strategy:
            self.model.labeling_strategy.pop(labeled_substrate)
    
    
    def set_flux(self, fluxid, value):
        '''
        Set metabolic flux value.
        
        Parameters
        ----------
        fluxid: str 
            Flux IDs, i.e., reaction ID + '_f' or '_b' for reversible reaction, 
            and reaction ID for irreversible reaction.
        value: float
            Flux value.
        '''
        
        self.model.total_fluxes[fluxid] = value
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_fluxes, fluxid))
        
    
    def set_fluxes_from_file(self, file):
        '''
        Read metabolic flux values from file.
        
        Parameters
        ----------
        file: file path
            tsv or excel file, fields are flux ID and value, 
            flux ID is reaction ID + '_f' or '_b' for reversible reaction, 
            and reaction ID for irreversible reaction.
        '''
        
        fluxes = read_preset_values_from_file(file)
        
        for fluxid, value in fluxes.iteritems():
            self.model.total_fluxes[fluxid] = value
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_fluxes, fluxes.index.tolist()))
        
                
    def _unset_fluxes(self, fluxids):
        '''
        Parameters
        ----------
        fluxids: str or list of str
            Flux ID(s).
        '''
        
        if not isinstance(fluxids, Iterable):
            fluxids = [fluxids]
        
        for fluxid in fluxids: 
            if fluxid in self.model.total_fluxes:
                self.model.total_fluxes.drop(fluxid, inplace = True)
    
    
    def _decompose_network(self, n_jobs, lump = True):
        '''
        Parameters
        ----------
        n_jobs: int
            # of jobs to run in parallel.
        lump: bool
            Whether to lump linear EMUs.
        '''
        
        if not self.model.target_EMUs:
            raise ValueError('call set_target_EMUs first')
        
        if not self.model.EAMs:
            if n_jobs <= 0:
                raise ValueError('n_jobs should be a positive value')    
            else:
                metabids = []
                atom_nos = []
                for emuid in self.model.target_EMUs:
                    metabid, atomNOs = emuid.split('_')
                    metabids.append(metabid)
                    atom_nos.append(atomNOs)
                
                EAMs = self.model._decompose_network(metabids, atom_nos, lump = lump, n_jobs = n_jobs)
                for size, EAM in EAMs.items():
                    self.model.EAMs[size] = EAM
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_decomposition)
            
            
    def _unset_decomposition(self):
        
        self.model.EAMs.clear()
            
    
    def _lambdify_matrix_As_and_Bs(self):
    
        if not self.model.matrix_As or not self.model.matrix_Bs:
            self.calculator._lambdify_matrix_As_and_Bs()
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_matrix_As_and_Bs)        
    
    
    def _unset_matrix_As_and_Bs(self):
        
        self.model.matrix_As.clear()
        self.model.matrix_Bs.clear()
            
    
    def _calculate_substrate_MDVs(self, extra_subs = None):
        '''
        Parameters
        ----------
        extra_subs: str or list of str
            Metabolite ID(s). Additional metabolites that are considered as substrates.
        '''

        if extra_subs is not None and not isinstance(extra_subs, list):
            extra_subs = [extra_subs]

        if not self.model.substrate_MDVs:
            self.calculator._calculate_substrate_MDVs(extra_subs)
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_substrate_MDVs)    
    
            
    def _unset_substrate_MDVs(self):
    
        self.model.substrate_MDVs.clear()    
    
    
    def prepare(self, n_jobs = 1):
        '''
        Parameters
        ----------
        n_jobs: int
            If n_jobs > 1, decomposition job will run in parallel.
        '''
        
        self._decompose_network(n_jobs)
        self._lambdify_matrix_As_and_Bs()
        self._calculate_substrate_MDVs()    
    
    
    def _check_dependencies(self):

        if self.model.total_fluxes.empty:
            raise ValueError('call set_flux or set_fluxes_from_file first')
        if not self.model.target_EMUs:
            raise ValueError('call set_target_EMUs first')
        if not self.model.labeling_strategy:
            raise ValueError('call labeling_strategy first')
        
        if any([not self.model.EAMs, 
                not self.model.substrate_MDVs, 
                not self.model.matrix_As, 
                not self.model.matrix_Bs]):
            raise ValueError('call prepare first')


    def simulate(self):
        
        self._check_dependencies()
        
        simMDVs = self.calculator._calculate_MDVs()
        targetMDVs = {emuid: MDV(simMDVs[emuid]) for emuid in self.model.target_EMUs}
        
        return SimResults(targetMDVs)
    