#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '06/14/2022'




from functools import partial
from collections.abc import Iterable
from copy import deepcopy
from math import ceil
import numpy as np
from multiprocessing import Pool
from ..io.inputs import read_measurements_from_file, read_initial_values
from ..io.results import InstFitResults, InstFitMCResults
from .inst_simulate import InstSimulator
from .fit import Fitter
from ..solver.nlpsolver import InstMFAModel
from ..utils.progress import Progress




class InstFitter(Fitter, InstSimulator):
    '''
    Estimated fluxes are in the unit of umol/gCDW/s if concentrations in the unit of 
    umol/gCDW and timepoints in the unit of s.
    '''
    
    def set_measured_MDVs(self, fragmentid, timepoints, means, sds):
        '''
        Set measured MDVs at various timepoints.
        
        Parameters
        ----------
        fragmentid: str
            Metabolite ID + '_' + atom NOs, e.g., 'Glu_12345'.
        timepoints: float or list of float
            Timepoint(s).
        means: array or list of array
            Mean of measured MDV(s). len(means) should be equal to len(timepoints).
        sds: array or list of array
            Standard deviation of measured MDV(s). len(sds) should be equal to len(timepoints).
        '''
        
        if not isinstance(timepoints, Iterable):
            timepoints = [timepoints]
            means = [means]
            sds = [sds]
            
        for timepoint, mean, sd in zip(timepoints, means, sds):
            self.model.measured_inst_MDVs.setdefault(fragmentid, {})[timepoint] = [np.array(mean), np.array(sd)]
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_measured_MDVs, {fragmentid: timepoints}))
            
    
    def set_measured_MDVs_from_file(self, file):
        '''
        Read measured MDVs at various timepoints from file.
        
        Parameters
        ----------
        file: file path
            Path of tsv or excel file with fields "fragment_ID", "time", "mean" and "sd".
            "fragment_ID" is metabolite ID + '_' + atom NOs, e.g., 'Glu_12345';
            "time" is timepoint when MDVs are measured (while some timepoints could be missing);
            "mean" and "sd" are the mean and standard deviation of MDV with element seperated by ','.

            Header line starts with "#", and will be skiped.
        '''
        
        measMDVs = read_measurements_from_file(file, inst_data = True)
        
        fragmentid_tpoints = {}
        for [emuid, timepoint], [mean, sd] in measMDVs.iterrows():
            timepoint = float(timepoint)
            self.model.measured_inst_MDVs.setdefault(emuid, {})[timepoint] = [np.array(list(map(float, mean.split(',')))), 
                                                                              np.array(list(map(float, sd.split(','))))]
            fragmentid_tpoints.setdefault(emuid, []).append(timepoint)
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_measured_MDVs, fragmentid_tpoints))
    
    
    def _unset_measured_MDVs(self, fragmentid_tpoints):
        '''
        Parameters
        ----------
        fragmentid_tpoints: dict
            measured MDV ID => list of timepoints.
        '''
        
        for fragmentid, timepoints in fragmentid_tpoints.items():
            if fragmentid in self.model.measured_inst_MDVs:
                for timepoint in timepoints:
                    if timepoint in self.model.measured_inst_MDVs[fragmentid]:
                        self.model.measured_inst_MDVs[fragmentid].pop(timepoint)
                        
                if not self.model.measured_inst_MDVs[fragmentid]:
                    self.model.measured_inst_MDVs.pop(fragmentid)
                    
            
    def set_concentration_bounds(self, metabid, bounds):
        '''
        Set lower and upper bounds for concentration in unit of umol/gCDW.
        
        Parameters
        ----------
        metabid: str or 'all'
            Metabolite ID. If 'all', all concentrations will be set to the range.
        bounds: 2-list
            [lower bound, upper bound]. Lower bound is not allow to equal upper bound.
        '''
        
        bounds = list(map(float, bounds))

        metabids = []
        if bounds[0] < bounds[1]:   # lower bound not allow to equal upper bound
            if metabid == 'all':
                for metabid in self.model.metabolites:
                    self.model.concentrations_bounds[metabid] = [max(0.0, bounds[0]), bounds[1]]
                    metabids.append(metabid)
            elif metabid in self.model.metabolites:
                self.model.concentrations_bounds[metabid] = [max(0.0, bounds[0]), bounds[1]]
                metabids = [metabid]
            else:
                raise ValueError('concentration range set to nonexistent metabolite %s' % metabid)
        else:
            raise ValueError('concentration lower bound should be less than upper bound')    
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_concentration_bounds, metabids))
            
            
    def _set_default_concentration_bounds(self):
        '''
        This method assign bounds of [0.01, 100] (umol/gCDW) for concentrations
        not set by set_concentration_bounds
        '''
        
        defBnds = [0.0001, 1]
        
        metabids = []
        for metabid in self.model.metabolites:
            if metabid not in self.model.net_fluxes_bounds:
                self.model.concentrations_bounds[metabid] = defBnds
                metabids.append(metabid)
                
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_concentration_bounds, metabids))        
        
    
    def _unset_concentration_bounds(self, metabids):
        '''
        Parameters
        ----------
        metabids: str or list of str
            Metabolite ID(s).
        '''
        
        if not isinstance(metabids, Iterable):
            metabids = [metabids]
        
        for metabid in metabids:
            if metabid in self.model.concentrations_bounds:
                self.model.concentrations_bounds.pop(metabid)
                
                
    def _decompose_network(self, n_jobs):
        '''
        Parameters
        ----------
        n_jobs: int
            # of jobs to run in parallel.
        '''
        
        if not self.model.measured_inst_MDVs:
            raise ValueError('call set_measured_MDV or set_measured_MDVs_from_file first')
        
        if not self.model.EAMs:
            if n_jobs <= 0:
                raise ValueError('n_jobs should be a positive value')    
            else:
                self.model.target_EMUs = list(self.model.measured_inst_MDVs.keys())
                
                metabids = []
                atom_nos = []
                for emuid in self.model.target_EMUs:
                    metabid, atomNOs = emuid.split('_')
                    metabids.append(metabid)
                    atom_nos.append(atomNOs)
                
                EAMs = self.model._decompose_network(metabids, atom_nos, lump = False, n_jobs = n_jobs)
                for size, EAM in EAMs.items():
                    self.model.EAMs[size] = EAM
                
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_decomposition)
            
    
    def _set_timepoints(self):
        
        if not self.model.timepoints:
            self.calculator._set_timepoints()
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_timepoints)
        
        
    def _calculate_matrix_Ms_derivatives_p(self):
        
        if not self.model.matrix_Ms_der_p:
            self.calculator._calculate_matrix_Ms_derivatives_p()
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_matrix_Ms_derivatives_p)
            
        
    def _unset_matrix_Ms_derivatives_p(self):
        
        self.model.matrix_Ms_der_p.clear()    
        
    
    def _calculate_measured_inst_MDVs_inversed_covariance_matrix(self):
        
        if not self.model.measured_inst_MDVs_inv_cov:
            self.calculator._calculate_measured_inst_MDVs_inversed_covariance_matrix()
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_measured_inst_MDVs_inversed_covariance_matrix)
            
            
    def _unset_measured_inst_MDVs_inversed_covariance_matrix(self):
        
        self.model.measured_inst_MDVs_inv_cov = None    
        
    
    def _calculate_initial_matrix_Xs_derivatives_p(self):
        
        if not self.model.initial_matrix_Xs_der_p:
            self.calculator._calculate_initial_matrix_Xs_derivatives_p()
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_initial_matrix_Xs_derivatives_p)
    
    
    def _unset_initial_matrix_Xs_derivatives_p(self):
        
        self.model.initial_matrix_Xs_der_p.clear()
    
    
    def _calculate_initial_matrix_Ys_derivatives_p(self):
        
        if not self.model.initial_matrix_Ys_der_p:
            self.calculator._calculate_initial_matrix_Ys_derivatives_p()
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_initial_matrix_Ys_derivatives_p)
    
    
    def _unset_initial_matrix_Ys_derivatives_p(self):
        
        self.model.initial_matrix_Ys_der_p.clear()
        
        
    def _estimate_concentrations_range(self):
        
        if not self.model.concentrations_range:
            for metabid in self.model.concids:
                self.model.concentrations_range[metabid] = self.model.concentrations_bounds[metabid]
                
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_concentrations_range, self.model.concids))
            
        
    def _unset_concentrations_range(self, metabids):
        '''
        Parameters
        ----------
        metabids: str or list of str
            Metabolite ID(s).
        '''
        
        if not isinstance(metabids, Iterable):
            metabids = [metabids]
        
        for metabid in metabids:
            if metabid in self.model.concentrations_range:
                self.model.concentrations_range.pop(metabid)
        
        
    def prepare(self, dilution_from = None, n_jobs = 1):
        '''
        Parameters
        ----------
        dilution_from: str or list of str
            ID(s) of unlabeled (inactive) metabolite leading to dilution effect.
            These metabolites have zero stoichiometric coefficients in reaction network.
        n_jobs: int
            If n_jobs > 1, decomposition job will run in parallel.
        '''    
        
        self._decompose_network(n_jobs)
        self._set_timepoints()
        self._calculate_null_space()
        self._calculate_transform_matrix()
        self._lambdify_matrix_As_and_Bs()
        self._calculate_matrix_As_and_Bs_derivatives_p('inst', n_jobs)
        self._lambdify_matrix_Ms()
        self._calculate_matrix_Ms_derivatives_p()
        self._calculate_substrate_MDVs(dilution_from)
        self._calculate_substrate_MDV_derivatives_p('inst', dilution_from)
        self._calculate_measured_inst_MDVs_inversed_covariance_matrix()
        self._calculate_measured_fluxes_inversed_covariance_matrix()
        self._calculate_measured_fluxes_derivative_p('inst')
        self._calculate_initial_matrix_Xs()
        self._calculate_initial_matrix_Ys()
        self._calculate_initial_matrix_Xs_derivatives_p()
        self._calculate_initial_matrix_Ys_derivatives_p()
        self._set_default_flux_bounds()
        self._estimate_fluxes_range(self.model.unbalanced_metabolites)
        self._set_default_concentration_bounds()
        self._estimate_concentrations_range()
    
        
    def _check_dependencies(self, fit_measured_fluxes):
        '''
        Parameters
        ----------
        fit_measured_fluxes: bool
            Whether to fit measured fluxes.
        '''

        if not self.model.net_fluxes_bounds:
            raise ValueError('call set_flux_bounds first')
        if not self.model.concentrations_bounds:
            raise ValueError('call set_concentration_bounds first')
        if not self.model.measured_inst_MDVs:
            raise ValueError('call call set_measured_MDV or set_measured_MDVs_from_file first')
        if not self.model.measured_fluxes:
            raise ValueError('call set_measured_flux or set_measured_fluxes_from_file first')    
        if not self.model.labeling_strategy:
            raise ValueError('call labeling_strategy first')
            
        checklist = [not self.model.target_EMUs, 
                     self.model.transform_matrix is None, 
                     self.model.null_space is None,
                     self.model.measured_inst_MDVs_inv_cov is None, 
                     not self.model.matrix_As, 
                     not self.model.matrix_Bs,
                     not self.model.matrix_Ms, 
                     not self.model.matrix_As_der_p, 
                     not self.model.matrix_Bs_der_p, 
                     not self.model.matrix_Ms_der_p, 
                     not self.model.substrate_MDVs, 
                     not self.model.substrate_MDVs_der_p,
                     self.model.measured_fluxes_der_p is None, 
                     not self.model.initial_matrix_Xs, 
                     not self.model.initial_matrix_Ys,
                     not self.model.initial_matrix_Xs_der_p, 
                     not self.model.initial_matrix_Ys_der_p, 
                     not self.model.timepoints]
        if fit_measured_fluxes:
            checklist.append(self.model.measured_fluxes_inv_cov is None)
        
        if any(checklist):
            raise ValueError('call prepare first')


    def solve(self, fit_measured_fluxes = True, ini_fluxes = None, ini_concs = None, 
              solver = 'slsqp', tol = 1e-6, max_iters = 400, show_progress = True):
        '''
        Parameters
        ----------
        fit_measured_fluxes: bool
            Whether to fit measured fluxes.
        ini_fluxes: ser or file in .tsv or .xlsx
            Initial values of net fluxes.
        ini_concs: ser or file in .tsv or .xlsx
            Initial values of concentrations.
        solvor: {"slsqp", "ralg"}
            * If "slsqp", scipy.optimize.minimze will be used.
            * If "ralg", openopt NLP solver will be used.
        tol: float
            Tolerance for termination.
        max_iters: int
            Maximum # of iterations.
        show_progress: bool
            Whether to show the progress bar. 
        '''            
                
        self._check_dependencies(fit_measured_fluxes)

        if ini_fluxes is not None:
            iniFluxes = read_initial_values(ini_fluxes, self.model.netfluxids)
        else:
            iniFluxes = ini_fluxes

        if ini_concs is not None:
            iniConcs = read_initial_values(ini_concs, self.model.concids)
        else:
            iniConcs = ini_concs
        
        optModel = InstMFAModel(self.model, fit_measured_fluxes, solver)    
        optModel.build_objective()
        optModel.build_gradient()
        optModel.build_flux_and_conc_bound_constraints()
        optModel.build_initial_flux_and_conc_values(ini_netfluxes = iniFluxes, ini_concs = iniConcs)
        
        with Progress('INST fitting', silent = not show_progress):
            res = optModel.solve_flux(tol, max_iters)

        return InstFitResults(*res[:8], deepcopy(res[8]), res[9], deepcopy(res[10]), *res[11:])


    def _solve_with_confidence_intervals(self, fit_measured_fluxes, ini_fluxes, ini_concs, 
                                         solver, tol, max_iters, nruns):
        '''
        Parameters
        ----------
        fit_measured_fluxes: bool
            Whether to fit measured fluxes.
        ini_fluxes: ser or file in .tsv or .xlsx or None
            Initial values of net fluxes.
        ini_concs: ser or file in .tsv or .xlsx or None
            Initial values of concentrations.    
        solvor: {"slsqp", "ralg"}
            * If "slsqp", scipy.optimize.minimze will be used.
            * If "ralg", openopt NLP solver will be used.
        tol: float
            Tolerance for termination.
        max_iters: int
            Maximum # of iterations.
        nruns: int
            # of estimations in each worker.
        '''

        import platform
        if platform.system() == 'Linux':
            import os
            os.sched_setaffinity(os.getpid(), range(os.cpu_count()))

        self._lambdify_matrix_As_and_Bs()
        self._lambdify_matrix_Ms()
        
        if ini_fluxes is not None:
            iniFluxes = read_initial_values(ini_fluxes, self.model.netfluxids)
        else:
            iniFluxes = ini_fluxes

        # if ini_concs is not None:
        #     iniConcs = read_initial_values(ini_concs, self.model.concids)
        # else:
        #     iniConcs = ini_concs
        
        optTotalfluxesSet = []
        optNetfluxesSet = []
        optConcsSet = []
        for _ in range(nruns):
            self.calculator._generate_random_fluxes()
            self.calculator._generate_random_inst_MDVs()
            
            optModel = InstMFAModel(self.model, fit_measured_fluxes, solver)
            optModel.build_objective()
            optModel.build_gradient()
            optModel.build_flux_and_conc_bound_constraints()
            optModel.build_initial_flux_and_conc_values(ini_netfluxes = iniFluxes)
            
            while True:
                optTotalfluxes, optNetfluxes, optConcs, *_, isSuccess = optModel.solve_flux(tol, max_iters)
                if isSuccess:
                    break
            optTotalfluxesSet.append(optTotalfluxes)
            optNetfluxesSet.append(optNetfluxes)
            optConcsSet.append(optConcs)
            
            self.calculator._reset_measured_fluxes()
            self.calculator._reset_measured_inst_MDVs()
            
        return optTotalfluxesSet, optNetfluxesSet, optConcsSet


    def solve_with_confidence_intervals(self, fit_measured_fluxes = True, ini_fluxes = None, 
                                        ini_concs = None, solver = 'slsqp', tol = 1e-6, max_iters = 400, 
                                        n_runs = 100, n_jobs = 1, show_progress = True):
        '''
        Parameters
        ----------
        fit_measured_fluxes: bool
            Whether to fit measured fluxes.
        ini_fluxes: ser or file in .tsv or .xlsx
            Initial values of net fluxes.
        ini_concs: ser or file in .tsv or .xlsx
            Initial values of concentrations.   
        solvor: {"slsqp", "ralg"}
            * If "slsqp", scipy.optimize.minimze will be used.
            * If "ralg", openopt NLP solver will be used.
        tol: float
            Tolerance for termination.
        max_iters: int
            Max # of iterations.
        show_progress: bool
            Whether to show the progress bar.
        n_runs: int
            # of runs to estimate confidence intervals.
        n_jobs: int
            # of jobs to run in parallel.
        '''
        
        self._check_dependencies(fit_measured_fluxes)
        
        self._unset_matrix_As_and_Bs()
        self._unset_matrix_Ms()

        if n_runs <= n_jobs:
            nruns_worker = 1
        else:
            nruns_worker = ceil(n_runs/n_jobs)
        
        pool = Pool(processes = n_jobs)
        with Progress('INST fitting with CIs', silent = not show_progress):
            resSet = []
            for _ in range(n_jobs):
                res = pool.apply_async(func = self._solve_with_confidence_intervals, 
                                       args = (fit_measured_fluxes,
                                               ini_fluxes,
                                               ini_concs,
                                               solver,
                                               tol,
                                               max_iters,
                                               nruns_worker))
                resSet.append(res)
            
            pool.close()    
            pool.join()
        
        resSet = [res.get() for res in resSet]
        
        totalFluxesSet = []
        netFluxesSet = []
        concsSet = []
        for totalFluxesSubset, netFluxesSubset, concsSubset in resSet:
            totalFluxesSet.extend(totalFluxesSubset)
            netFluxesSet.extend(netFluxesSubset)
            concsSet.extend(concsSubset)
        
        return InstFitMCResults(totalFluxesSet, netFluxesSet, concsSet)     
        