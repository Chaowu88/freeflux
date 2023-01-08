#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '05/18/2022'




from collections import ChainMap
from collections.abc import Iterable
from functools import reduce
from copy import deepcopy
import warnings
warnings.filterwarnings('ignore', category = RuntimeWarning) 
import numpy as np
from numpy.random import normal
import pandas as pd
from scipy.linalg import null_space, pinv2, expm
from sympy import symbols, lambdify, Matrix, derive_by_array
from multiprocessing import Pool
from ..core.mdv import MDV, get_natural_MDV, get_substrate_MDV, conv, diff_conv



    
class Calculator():
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
        
    
    def _set_timepoints(self):
        
        ts = []
        for _, tpoint_infos in self.model.measured_inst_MDVs.items():
            ts.extend(tpoint_infos.keys())
        ts = sorted(set(ts))
        
        for t in ts:
            self.model.timepoints.append(t)
            
    
    def _calculate_null_space(self):
        
        S = self.model.get_total_stoichiometric_matrix(self.model.unbalanced_metabolites)
        self.model.null_space = null_space(S.values)
    
    
    def _calculate_transform_matrix(self):
        
        transMat = pd.DataFrame(0, index = self.model.netfluxids, columns = self.model.totalfluxids)
        for rxnid, rxn in self.model.reactions_info.items():
        
            if rxn.reversible:
                transMat.loc[rxnid, rxnid+'_f'] = 1
                transMat.loc[rxnid, rxnid+'_b'] = -1
            else:
                transMat.loc[rxnid, rxnid] = 1
        
        self.model.transform_matrix = transMat.values
    
    
    def _calculate_substrate_MDVs(self, extra_subs):
        
        if extra_subs is None:
             extra_subs = []

        for size in self.model.matrix_Bs:
            for sourceEMU in self.model.matrix_Bs[size][2]:
                
                if not isinstance(sourceEMU, Iterable):
                    sourceEMU = (sourceEMU,)
                    
                for emu in sourceEMU:
                    metabid = emu.metabolite_id
                    if metabid in self.model.end_substrates + extra_subs:
                        
                        if metabid in self.model.labeling_strategy:
                            atom_nos = emu.atom_nos
                            labeling_pattern, percentage, purity = self.model.labeling_strategy[metabid]
                            self.model.substrate_MDVs[emu] = get_substrate_MDV(atom_nos, labeling_pattern, percentage, purity)
                        else:
                            natoms = emu.size
                            self.model.substrate_MDVs[emu] = get_natural_MDV(natoms)                        
    
    
    def _calculate_substrate_MDV_derivatives_basic(self, nvars, extra_subs):
        
        if extra_subs is None:
             extra_subs = []

        substrate_MDVs_der = {}
        for size in self.model.matrix_Bs:
            for sourceEMU in self.model.matrix_Bs[size][2]:
                
                if not isinstance(sourceEMU, Iterable):
                    sourceEMU = (sourceEMU,)
                    
                for emu in sourceEMU:
                    metabid = emu.metabolite_id
                    if metabid in self.model.end_substrates + extra_subs:
                        substrate_MDVs_der[emu] = np.zeros((emu.size+1, nvars))
                        
        return substrate_MDVs_der                
                        
                        
    def _calculate_substrate_MDV_derivatives_u(self, extra_subs):
        
        substrate_MDVs_der = self._calculate_substrate_MDV_derivatives_basic(self.model.null_space.shape[1], extra_subs)
        
        return substrate_MDVs_der


    def _calculate_substrate_MDV_derivatives_c(self, extra_subs):
       
        substrate_MDVs_der = self._calculate_substrate_MDV_derivatives_basic(len(self.model.concids), extra_subs)
        
        return substrate_MDVs_der
        
        
    def _calculate_substrate_MDV_derivatives_p(self, kind, extra_subs = None):
        '''
        Parameters
        ----------
        kind: {"ss", "inst"}
            * "ss" if isotopic steady state.
            * "inst" if isotopically nonstationary state.
        extra_subs: list or None
        '''
        
        if kind == 'ss':
            substrate_MDVs_der_u = self._calculate_substrate_MDV_derivatives_u(extra_subs)
            for emu in substrate_MDVs_der_u:
                self.model.substrate_MDVs_der_p[emu] = substrate_MDVs_der_u[emu]

        elif kind == 'inst':
            substrate_MDVs_der_u = self._calculate_substrate_MDV_derivatives_u(extra_subs)
            substrate_MDVs_der_c = self._calculate_substrate_MDV_derivatives_c(extra_subs)
            for emu in substrate_MDVs_der_u:
                self.model.substrate_MDVs_der_p[emu] = np.concatenate((substrate_MDVs_der_u[emu], 
                                                                       substrate_MDVs_der_c[emu]), 
                                                                       axis = 1)
    
    
    def _calculate_measured_fluxes_inversed_covariance_matrix(self):
        
        sds = np.array([sd for _, sd in self.model.measured_fluxes.values()])
        self.model.measured_fluxes_inv_cov = np.diag(1/sds**2) 
    
   
    def _calculate_measured_fluxes_derivative_v(self):
        
        measfluxids = list(self.model.measured_fluxes.keys())
        measured_fluxes_der = np.array(Matrix(measfluxids).jacobian(Matrix(self.model.totalfluxids))).astype(float)
        
        return measured_fluxes_der
        
        
    def _calculate_measured_fluxes_derivative_c(self):
        
        nmeasfluxes = len(self.model.measured_fluxes)
        nmetabs = len(self.model.concids)
        measured_fluxes_der = np.zeros((nmeasfluxes, nmetabs))
        
        return measured_fluxes_der
        
        
    def _calculate_measured_fluxes_derivative_p(self, kind):
        '''
        Parameters
        ----------
        kind: {"ss", "inst"}
            * "ss" if isotopic steady state.
            * "inst" if isotopically nonstationary state.
        '''
        
        if kind == 'ss':
            measured_fluxes_der_v = self._calculate_measured_fluxes_derivative_v()
            measured_fluxes_der_u = measured_fluxes_der_v@self.model.null_space
            self.model.measured_fluxes_der_p = measured_fluxes_der_u

        elif kind == 'inst':
            measured_fluxes_der_v = self._calculate_measured_fluxes_derivative_v()
            measured_fluxes_der_u = measured_fluxes_der_v@self.model.null_space
            measured_fluxes_der_c = self._calculate_measured_fluxes_derivative_c()
            self.model.measured_fluxes_der_p = np.concatenate((measured_fluxes_der_u, measured_fluxes_der_c), axis = 1)

    
    def _generate_random_fluxes(self):
        
        self.ori_measured_fluxes = deepcopy(self.model.measured_fluxes)
        
        for fluxid, [mean, sd] in self.model.measured_fluxes.items():
            while True:
                meanNew = normal(mean, sd)
                if mean-3*sd <= meanNew <= mean+3*sd:
                    break
            self.model.measured_fluxes[fluxid][0] = meanNew
        
        
    def _reset_measured_fluxes(self):
        
        self.model.measured_fluxes = deepcopy(self.ori_measured_fluxes)
        
    
    def _calculate_measured_MDVs_inversed_covariance_matrix(self):
        
        sds = np.concatenate([sd for _, sd in self.model.measured_MDVs.values()])
        self.model.measured_MDVs_inv_cov = np.diag(1/sds**2)
        
        
    def _calculate_measured_inst_MDVs_inversed_covariance_matrix(self):
        
        sds = []
        for _, t_mdvs in self.model.measured_inst_MDVs.items():
            for t in sorted(t_mdvs):
                if t != 0:
                    sds.append(t_mdvs[t][1])
        sds = np.concatenate(sds)    
        
        self.model.measured_inst_MDVs_inv_cov = np.diag(1/sds**2)
        
    
    def _generate_random_MDVs(self):
        
        self.ori_measured_MDVs = deepcopy(self.model.measured_MDVs)
        
        for emuid, [means, sds] in self.model.measured_MDVs.items():
            mdv = MDV(normal(means, sds))
            self.model.measured_MDVs[emuid][0] = mdv.value


    def _reset_measured_MDVs(self):
        
        self.model.measured_MDVs = deepcopy(self.ori_measured_MDVs)
        

    def _generate_random_inst_MDVs(self):
        
        self.ori_measured_inst_MDVs = deepcopy(self.model.measured_inst_MDVs)

        for emuid, instMDVs in self.model.measured_inst_MDVs.items():
            for t, [means, sds] in instMDVs.items():
                if t != 0:
                    mdv = MDV(normal(means, sds))
                    self.model.measured_inst_MDVs[emuid][t][0] = mdv.value


    def _reset_measured_inst_MDVs(self):

        self.model.measured_inst_MDVs = deepcopy(self.ori_measured_inst_MDVs)
    

    @staticmethod
    def _calculate_matrix_A_and_B_derivatives_v(EAM, fluxids):
        '''
        Parameters
        ----------
        EAMstr: df
            EMU adjacency matrix of some size.
        fluxids: list
            Total fluxes IDs.
            
        Returns
        -------
        matrix_A_der, matrix_B_der: 3-D array
            Derivatives of A(B) w.r.t. total fluxes in shape of (len(total fluxes), A(B).shape[0], A(B).shape[1]).
        '''
        
        import platform
        if platform.system() == 'Linux':
            import os
            os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
            
        preAB = EAM.copy(deep = 'all').values
        colSum = preAB.sum(axis = 0)
        for i in range(colSum.size):
            preAB[i, i] = -colSum[i]
        
        A = preAB[:preAB.shape[1], :].T
        B = -preAB[preAB.shape[1]:, :].T
        
        matA = Matrix(A)
        matB = Matrix(B)
        
        matrix_A_der = np.array(derive_by_array(matA, symbols(fluxids)), dtype = float)
        matrix_B_der = np.array(derive_by_array(matB, symbols(fluxids)), dtype = float)
        
        return matrix_A_der, matrix_B_der
    
    
    def _calculate_matrix_As_and_Bs_derivatives_u(self, n_jobs):
        '''
        Parameters
        ----------
        n_jobs: int
            # of jobs to run in parallel.

        Returns
        -------
        matrix_A_der: 3-D array
            Derivatives of A(B) w.r.t. free fluxes in shape of (len(free fluxes), A(B).shape[0], A(B).shape[1]).
        matrix_B_der: 3-D array
            Derivatives of A(B) w.r.t. free fluxes in shape of (len(free fluxes), A(B).shape[0], A(B).shape[1]).        
        '''
            
        matrix_ABs_der = {}
        if n_jobs == 1:
            for size, EAM in self.model.EAMs.items():
                ABder = self._calculate_matrix_A_and_B_derivatives_v(EAM, self.model.totalfluxids)
                matrix_ABs_der[size] = ABder
        else:
            pool = Pool(processes = n_jobs)
            
            for size, EAM in self.model.EAMs.items():
                ABder = pool.apply_async(func = self._calculate_matrix_A_and_B_derivatives_v, 
                                         args = (EAM, self.model.totalfluxids))
                matrix_ABs_der[size] = ABder
            
            pool.close()    
            pool.join()
            
            matrix_ABs_der = {size: ABder.get() for size, ABder in matrix_ABs_der.items()}
        
        matrix_As_der = {}
        matrix_Bs_der = {}
        for size, ABder in matrix_ABs_der.items():
            Ader = ABder[0].swapaxes(0,1).swapaxes(1,2)
            Bder = ABder[1].swapaxes(0,1).swapaxes(1,2)
            
            Ader = Ader@self.model.null_space
            Bder = Bder@self.model.null_space
        
            Ader = Ader.swapaxes(1,2).swapaxes(0,1)
            Bder = Bder.swapaxes(1,2).swapaxes(0,1)
            
            matrix_As_der[size] = Ader
            matrix_Bs_der[size] = Bder
            
        return matrix_As_der, matrix_Bs_der
        
            
    def _calculate_matrix_As_and_Bs_derivatives_c(self):
        '''
        Returns
        -------
        matrix_A_der: 3-D array
            Derivatives of A(B) w.r.t. total fluxes in shape of (len(concs), A(B).shape[0], A(B).shape[1]).
        matrix_B_der: 3-D array
            Derivatives of A(B) w.r.t. total fluxes in shape of (len(concs), A(B).shape[0], A(B).shape[1]).    
        '''
        
        nmetabs = len(self.model.concids)
        
        matrix_As_der = {}
        matrix_Bs_der = {}
        for size, EAM in self.model.EAMs.items():
            nEMUsout = EAM.shape[1]
            nEMUsin = EAM.shape[0] - nEMUsout
            
            matrix_As_der[size] = np.zeros((nmetabs, nEMUsout, nEMUsout))
            matrix_Bs_der[size] = np.zeros((nmetabs, nEMUsout, nEMUsin))
        
        return matrix_As_der, matrix_Bs_der
        
        
    def _calculate_matrix_As_and_Bs_derivatives_p(self, kind, n_jobs):
        '''
        Parameters
        ----------
        kind: {"ss", "inst"}
            * "ss" if isotopic steady state.
            * "inst" if isotopically nonstationary state.
        n_jobs: int
            # of jobs to run in parallel.
        '''
        
        if kind == 'ss':
            matrix_As_der_u, matrix_Bs_der_u = self._calculate_matrix_As_and_Bs_derivatives_u(n_jobs)
            for size in self.model.EAMs:
                self.model.matrix_As_der_p[size] = matrix_As_der_u[size]
                self.model.matrix_Bs_der_p[size] = matrix_Bs_der_u[size]

        elif kind == 'inst':
            matrix_As_der_u, matrix_Bs_der_u = self._calculate_matrix_As_and_Bs_derivatives_u(n_jobs)
            matrix_As_der_c, matrix_Bs_der_c = self._calculate_matrix_As_and_Bs_derivatives_c()
            for size in self.model.EAMs:
                self.model.matrix_As_der_p[size] = np.concatenate((matrix_As_der_u[size], matrix_As_der_c[size]), axis = 0)
                self.model.matrix_Bs_der_p[size] = np.concatenate((matrix_Bs_der_u[size], matrix_Bs_der_c[size]), axis = 0)
    
    
    def _lambdify_matrix_As_and_Bs(self):
       
        for size, EAM in self.model.EAMs.items():
            
            preAB = EAM.copy(deep = 'all')
            for emu in EAM.columns:
                preAB.loc[emu, emu] = -preAB[emu].sum() 
        
            A = preAB.loc[preAB.columns, :].T
            B = -preAB.loc[preAB.index.difference(preAB.columns), :].T
            
            matA = Matrix(A)
            matB = Matrix(B)
            
            fluxidsA = list(map(str, matA.free_symbols))
            fluxidsB = list(map(str, matB.free_symbols))
            
            lambA = lambdify(fluxidsA, matA, modules = 'numpy')
            lambB = lambdify(fluxidsB, matB, modules = 'numpy')
            
            self.model.matrix_As[size] = [lambA, fluxidsA, A.columns.tolist()]
            self.model.matrix_Bs[size] = [lambB, fluxidsB, B.columns.tolist()]
    
    def _calculate_matrix_Ms_derivatives_u(self):
        
        nfreefluxes = self.model.null_space.shape[1]
        
        matrix_Ms_der = {}
        for size, EAM in self.model.EAMs.items():
            nEMUsout = EAM.shape[1]
            matrix_Ms_der[size] = np.zeros((nfreefluxes, nEMUsout, nEMUsout))
        
        return matrix_Ms_der
        
    
    def _calculate_matrix_Ms_derivatives_c(self):
        
        matrix_Ms_der = {}
        for size, EAM in self.model.EAMs.items():            
            matM = Matrix(np.diag(symbols([emu.metabolite_id for emu in EAM.columns])))
            matrix_M_der = np.array(derive_by_array(matM, symbols(self.model.concids)), dtype = float)
            matrix_Ms_der[size] = matrix_M_der
            
        return matrix_Ms_der
        
        
    def _calculate_matrix_Ms_derivatives_p(self):
        
        matrix_Ms_der_u = self._calculate_matrix_Ms_derivatives_u()
        matrix_Ms_der_c = self._calculate_matrix_Ms_derivatives_c()
        
        for size in self.model.EAMs:
            self.model.matrix_Ms_der_p[size] = np.concatenate((matrix_Ms_der_u[size], matrix_Ms_der_c[size]), axis = 0)
        
    
    def _lambdify_matrix_Ms(self):
        
        for size, EAM in self.model.EAMs.items():
            matM = Matrix(np.diag(symbols([emu.metabolite_id for emu in EAM.columns])))
            metabids = list(map(str, matM.free_symbols))
            lambM = lambdify(metabids, matM, modules = 'numpy')
            self.model.matrix_Ms[size] = [lambM, metabids]
        
    
    def _calculate_initial_matrix_Xs(self):
        
        for size in self.model.matrix_As:
            nEMUs = len(self.model.matrix_As[size][2])
            iniX = np.vstack([get_natural_MDV(size).value] * nEMUs)
            self.model.initial_matrix_Xs[size] = iniX
            
            
    def _calculate_initial_matrix_Ys(self):
        
        for size in self.model.matrix_Bs:
            iniY = []
            for sourceEMU in self.model.matrix_Bs[size][2]:
                if not isinstance(sourceEMU, Iterable):
                    sourceMDV = self.model.substrate_MDVs[sourceEMU]
                else:
                    mdvs = []
                    for emu in sourceEMU:
                        if emu not in self.model.substrate_MDVs:
                            mdv = get_natural_MDV(emu.size)
                        else:
                            mdv = self.model.substrate_MDVs[emu]
                        mdvs.append(mdv)
                    sourceMDV = reduce(conv, mdvs)
                iniY.append(sourceMDV)
            iniY = np.array(iniY)
            self.model.initial_matrix_Ys[size] = iniY
        
    
    def _calculate_initial_matrix_Xs_derivatives_p(self):
        
        nvars = self.model.null_space.shape[1] + len(self.model.concids)
        for size, iniX in self.model.initial_matrix_Xs.items():
            Xshape = iniX.shape
            iniXder = np.zeros((nvars, *Xshape))
            self.model.initial_matrix_Xs_der_p[size] = iniXder
        
        
    def _calculate_initial_matrix_Ys_derivatives_p(self):
        
        nvars = self.model.null_space.shape[1] + len(self.model.concids)
        for size, iniY in self.model.initial_matrix_Ys.items():
            Yshape = iniY.shape
            iniYder = np.zeros((nvars, *Yshape))
            self.model.initial_matrix_Ys_der_p[size] = iniYder
        
    
    def _build_initial_sim_MDVs(self):
        
        for size in sorted(self.model.matrix_As):
            productEMUs = self.model.matrix_As[size][2]
            iniX = self.model.initial_matrix_Xs[size]        
            for productEMU, iniMDV in zip(productEMUs, iniX):
                if productEMU.id in self.model.target_EMUs:
                    self.model.initial_sim_MDVs[productEMU.id] = {0: MDV(iniMDV)}
                    
        
    def _calculate_MDVs(self):
        '''
        This method simulate MDVs at isotopically steady state.
        
        Returns
        -------
        simMDVs: dict
            EMU ID => MDV (in array).
        '''
        
        simMDVs = {}
        for size in sorted(self.model.matrix_As):
            
            lambA, fluxidsA, productEMUs = self.model.matrix_As[size]
            lambB, fluxidsB, sourceEMUs = self.model.matrix_Bs[size]
            
            A = lambA(*self.model.total_fluxes[fluxidsA])
            B = lambB(*self.model.total_fluxes[fluxidsB])
            
            Y = []
            for sourceEMU in sourceEMUs:
                if not isinstance(sourceEMU, Iterable):
                    sourceMDV = self.model.substrate_MDVs[sourceEMU]
                else:
                    mdvs = [ChainMap(simMDVs, self.model.substrate_MDVs)[emu] for emu in sourceEMU]
                    sourceMDV = reduce(conv, mdvs)
                Y.append(sourceMDV)    
            Y = np.array(Y)
            
            X = pinv2(A, check_finite = True)@B@Y
        
            simMDVs.update(zip(productEMUs, X))
            
        simMDVs = {emu.id: mdv for emu, mdv in simMDVs.items()}
        
        return simMDVs
                    
    
    def _calculate_MDVs_and_derivatives_p(self):
        '''
        This method simulate MDVs and their derivatives at isotopically steady state.
        
        Returns
        -------
        simMDVs: dict
            EMU ID => MDV (in array).
        simMDVsDer: dict
            EMU ID => 2-D array in shape of (len(MDV), len(u)).
        '''
        
        simMDVs = {}
        simMDVsDer = {}
        for size in sorted(self.model.matrix_As):
            
            lambA, fluxidsA, productEMUs = self.model.matrix_As[size]
            lambB, fluxidsB, sourceEMUs = self.model.matrix_Bs[size]
            
            A = lambA(*self.model.total_fluxes[fluxidsA])
            B = lambB(*self.model.total_fluxes[fluxidsB])
            
            Ainv = pinv2(A, check_finite = True)
            
            Ader = self.model.matrix_As_der_p[size]
            Bder = self.model.matrix_Bs_der_p[size]
            
            Y = []
            Yder = []
            for sourceEMU in sourceEMUs:
                if not isinstance(sourceEMU, Iterable):
                    sourceMDV = self.model.substrate_MDVs[sourceEMU]
                    sourceMDVder = self.model.substrate_MDVs_der_p[sourceEMU]
                else:
                    mdvs = []
                    mdvs_mdvders = []
                    for emu in sourceEMU:
                        mdv = ChainMap(simMDVs, self.model.substrate_MDVs)[emu]
                        mdvder = ChainMap(simMDVsDer, self.model.substrate_MDVs_der_p)[emu]
                        mdvs.append(mdv)
                        mdvs_mdvders.append([mdv, mdvder])
                    sourceMDV = reduce(conv, mdvs)
                    sourceMDVder = reduce(diff_conv, mdvs_mdvders)[1]
                Y.append(sourceMDV)
                Yder.append(sourceMDVder)
            Y = np.array(Y)
            Yder = np.array(Yder).swapaxes(1,2).swapaxes(0,1)
            
            X = Ainv@B@Y
            Xder = Ainv@(Bder@Y + B@Yder - Ader@X)
            Xder = Xder.swapaxes(0,1).swapaxes(1,2)
            
            simMDVs.update(zip(productEMUs, X))
            simMDVsDer.update(zip(productEMUs, Xder))
            
        simMDVs = {emu.id: mdv for emu, mdv in simMDVs.items()}
        simMDVsDer = {emu.id: mdvDer for emu, mdvDer in simMDVsDer.items()}
        
        return simMDVs, simMDVsDer
    
    
    def _calculate_inst_MDVs(self):
        '''
        This method simulate MDVs at isotopically nonstationary state.
        
        Returns
        -------
        simInstMDVs: dict
            EMU ID => {t => MDV (in array)} (starting from t1).
        '''
        
        simInstMDVs = {}
        Ys = {}
        Xs = {}

        t1 = 0.0
        for size in sorted(self.model.matrix_As):
            Y_t1 = self.model.initial_matrix_Ys[size]
            Ys.setdefault(t1, {})[size] = Y_t1
            
            X_t1 = self.model.initial_matrix_Xs[size]
            Xs.setdefault(t1, {})[size] = X_t1

        for t in self.model.timepoints:
            if t != 0.0:
                t0 = t1
                t1 = t            
                deltat = t1 - t0
            
                for size in sorted(self.model.matrix_As):
                    lambA, fluxidsA, productEMUs = self.model.matrix_As[size]
                    lambB, fluxidsB, sourceEMUs = self.model.matrix_Bs[size]
                    lambM, metabids =  self.model.matrix_Ms[size]
                    
                    A = lambA(*self.model.total_fluxes[fluxidsA])
                    B = lambB(*self.model.total_fluxes[fluxidsB])
                    M = lambM(*self.model.concentrations[metabids])
                    Minv = pinv2(M, check_finite = True)
                    
                    F = Minv@A
                    Finv = pinv2(F, check_finite = True)
                    I = np.eye(*F.shape)
                    Phi = expm(F*deltat)
                    Gamma = (Phi - I)@Finv
                    Omega = (Gamma/deltat - I)@Finv
                    
                    X_t0 = Xs[t0][size]
                    
                    Y_t0 = Ys[t0][size]
                    G_t0 = Minv@B@Y_t0
                    
                    Y_t1 = []
                    for sourceEMU in sourceEMUs:
                        if not isinstance(sourceEMU, Iterable):
                            sourceMDV = self.model.substrate_MDVs[sourceEMU]
                        else:
                            mdvs = []
                            for emu in sourceEMU:
                                mdv = ChainMap(simInstMDVs, self.model.substrate_MDVs)[emu]
                                if isinstance(mdv, dict):
                                    mdv = mdv[t1]
                                mdvs.append(mdv)
                            sourceMDV = reduce(conv, mdvs)
                        Y_t1.append(sourceMDV)    
                    Y_t1 = np.array(Y_t1)
                    G_t1 = Minv@B@Y_t1
                    
                    X_t1 = Phi@X_t0 - Gamma@G_t0 - Omega@(G_t1 - G_t0)
                    
                    Ys.setdefault(t1, {})[size] = Y_t1
                    Xs.setdefault(t1, {})[size] = X_t1
                    
                    for productEMU, mdv_t1 in zip(productEMUs, X_t1):
                        simInstMDVs.setdefault(productEMU, {}).update({t1: mdv_t1})
                    
        simInstMDVs = {emu.id: mdvs for emu, mdvs in simInstMDVs.items()}
        
        return simInstMDVs        


    def _calculate_inst_MDVs_and_derivatives_p(self):
        '''
        This method simulate MDVs and their derivatives at isotopically nonstationary state.
        
        Returns
        -------
        simInstMDVs: dict
            EMU ID => {t => MDV (in array)} (starting from t1).
        simInstMDVsDer: dict
            EMU ID => {t => 2-D array in shape of (len(MDV), len(u)+len(c))} (starting from t1).
        '''
        
        simInstMDVs = {}
        simInstMDVsDer = {}
        Ys = {}
        Xs = {}
        Yders = {}
        Xders = {}
        
        t1 = 0.0
        for size in sorted(self.model.matrix_As):
            Y_t1 = self.model.initial_matrix_Ys[size]
            Ys.setdefault(t1, {})[size] = Y_t1
            
            X_t1 = self.model.initial_matrix_Xs[size]
            Xs.setdefault(t1, {})[size] = X_t1
            
            Yder_t1 = self.model.initial_matrix_Ys_der_p[size]
            Yders.setdefault(t1, {})[size] = Yder_t1
            
            Xder_t1 = self.model.initial_matrix_Xs_der_p[size]
            Xders.setdefault(t1, {})[size] = Xder_t1
        
        for t in self.model.timepoints:
            if t != 0.0:
                t0 = t1
                t1 = t            
                deltat = t1 - t0
            
                for size in sorted(self.model.matrix_As):
                    lambA, fluxidsA, productEMUs = self.model.matrix_As[size]
                    lambB, fluxidsB, sourceEMUs = self.model.matrix_Bs[size]
                    lambM, metabids =  self.model.matrix_Ms[size]
            
                    A = lambA(*self.model.total_fluxes[fluxidsA])
                    B = lambB(*self.model.total_fluxes[fluxidsB])
                    M = lambM(*self.model.concentrations[metabids])
                    Minv = pinv2(M, check_finite = True)
                    
                    Ader = self.model.matrix_As_der_p[size]
                    Bder = self.model.matrix_Bs_der_p[size]
                    Mder = self.model.matrix_Ms_der_p[size]
                    Minvder = -Minv@Mder@Minv
                    
                    F = Minv@A
                    Finv = pinv2(F, check_finite = True)
                    I = np.eye(*F.shape)
                    Phi = expm(F*deltat)
                    Gamma = (Phi - I)@Finv
                    Omega = (Gamma/deltat - I)@Finv
                    
                    X_t0 = Xs[t0][size]
                    
                    Y_t0 = Ys[t0][size]
                    G_t0 = Minv@B@Y_t0
                    
                    Xder_t0 = Xders[t0][size]
                    
                    Yder_t0 = Yders[t0][size]
                    H_t0 = Minvder@A@X_t0 + Minv@Ader@X_t0 - Minv@B@Yder_t0 - Minvder@B@Y_t0 - Minv@Bder@Y_t0
                    
                    Y_t1 = []
                    Yder_t1 = []
                    for sourceEMU in sourceEMUs:
                        if not isinstance(sourceEMU, Iterable):
                            sourceMDV = self.model.substrate_MDVs[sourceEMU]
                            sourceMDVder = self.model.substrate_MDVs_der_p[sourceEMU]
                        else:
                            mdvs = []
                            mdvs_mdvders = []
                            for emu in sourceEMU:
                                mdv = ChainMap(simInstMDVs, self.model.substrate_MDVs)[emu]
                                if isinstance(mdv, dict): mdv = mdv[t1]
                                mdvs.append(mdv)    
                                mdvder = ChainMap(simInstMDVsDer, self.model.substrate_MDVs_der_p)[emu]
                                if isinstance(mdvder, dict): mdvder = mdvder[t1]
                                mdvs_mdvders.append([mdv, mdvder])
                            sourceMDV = reduce(conv, mdvs)
                            sourceMDVder = reduce(diff_conv, mdvs_mdvders)[1]
                        Y_t1.append(sourceMDV)
                        Yder_t1.append(sourceMDVder)
                    Y_t1 = np.array(Y_t1)
                    Yder_t1 = np.array(Yder_t1).swapaxes(1,2).swapaxes(0,1)
                    
                    G_t1 = Minv@B@Y_t1
                    X_t1 = Phi@X_t0 - Gamma@G_t0 - Omega@(G_t1 - G_t0)
                    
                    H_t1 = Minvder@A@X_t1 + Minv@Ader@X_t1 - Minv@B@Yder_t1 - Minvder@B@Y_t1 - Minv@Bder@Y_t1
                    Xder_t1 = Phi@Xder_t0 + Gamma@H_t0 + Omega@(H_t1 - H_t0)
                    
                    Ys.setdefault(t1, {})[size] = Y_t1
                    Xs.setdefault(t1, {})[size] = X_t1
                    Yders.setdefault(t1, {})[size] = Yder_t1
                    Xders.setdefault(t1, {})[size] = Xder_t1
                    
                    for productEMU, mdv_t1 in zip(productEMUs, X_t1):
                        simInstMDVs.setdefault(productEMU, {}).update({t1: mdv_t1})
                        
                    Xder_t1 = Xder_t1.swapaxes(0,1).swapaxes(1,2)    
                    for productEMU, mdvder_t1 in zip(productEMUs, Xder_t1):
                        simInstMDVsDer.setdefault(productEMU, {}).update({t1: mdvder_t1})
                        
        simInstMDVs = {emu.id: mdvs for emu, mdvs in simInstMDVs.items()}                
        simInstMDVsDer = {emu.id: mdvders for emu, mdvders in simInstMDVsDer.items()}
        
        return simInstMDVs, simInstMDVsDer
