#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '05/24/2022'




from collections.abc import Iterable
from functools import partial
from copy import deepcopy
from math import ceil
import numpy as np
from multiprocessing import Pool
from ..optim.optim import Optimizer
from .simulate import Simulator
from ..io.inputs import read_measurements_from_file, read_initial_values
from ..io.results import FitResults, FitMCResults
from ..utils.utils import Calculator
from ..utils.progress import Progress
from ..solver.nlpsolver import MFAModel 




class Fitter(Optimizer, Simulator):
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
        
        super().__init__(model)
        
        self.calculator = Calculator(self.model)
        
    
    def set_measured_MDV(self, fragmentid, mean, sd):
        '''
        Set measured MDV.
        
        Parameters
        ----------
        fragmentid: str
            Metabolite ID + "_" + atom NOs, e.g., "Glu_12345".
        mean: array
            Means of measured MDV vector.
        sd: array
            Standard deviations of measured MDV vector.
        '''
        
        self.model.measured_MDVs[fragmentid] = [np.array(mean), np.array(sd)]
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_measured_MDVs, fragmentid))
            
            
    def set_measured_MDVs_from_file(self, file):
        '''
        Read measured MDVs from file.
        
        Parameters
        ----------
        file: file path
            tsv or excel file with fields "fragment_ID", "mean" and "sd". 
            "fragment_ID" is metabolite ID + "_" + atom NOs, e.g., 'Glu_12345'; 
            "mean" and "sd" are the mean and standard deviation of MDV with 
            element seperated by ",".

            Header line starts with "#", and will be skiped.
        '''
        
        measMDVs = read_measurements_from_file(file)

        for emuid, [mean, sd] in measMDVs.iterrows():
            self.model.measured_MDVs[emuid] = [np.array(list(map(float, mean.split(',')))),
                                               np.array(list(map(float, sd.split(','))))]
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_measured_MDVs, measMDVs.index.tolist()))
        
        
    def _unset_measured_MDVs(self, fragmentids):
        '''
        Parameters
        ----------
        fragmentids: str or list of str
            Measured MDV ID(s).
        '''
        
        if not isinstance(fragmentids, Iterable):
            fragmentids = [fragmentids]
            
        for fragmentid in fragmentids:
            if fragmentid in self.model.measured_MDVs:
                self.model.measured_MDVs.pop(fragmentid)
        
        
    def set_measured_flux(self, fluxid, mean, sd):
        '''
        Set measured flux.
        
        Parameters
        ----------
        fluxid: str
            Flux ID, i.e., reaction ID. Typically, measured fluxes are substrate 
            consumption product formation or cell growth, etc. They should be irreversible.
        mean: float
            Mean of measured flux.
        sd: float
            Standard deviation of measured flux.
        '''
        
        self.model.measured_fluxes[fluxid] = [mean, sd]
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_measured_fluxes, fluxid))
            
            
    def set_measured_fluxes_from_file(self, file):
        '''
        Read measured fluxes from file.
        
        Parameters
        ----------
        file: file path
            tsv or excel file with fields "reaction_ID", "mean" and "sd".
            "reaction_ID" is reaction ID, typically measured fluxes are substrate 
            consumption, product formation or cell growth, etc. They should be irreversible;
            "mean" and "sd" are the mean and standard deviation of measured flux.

            Header line starts with "#", and will be skiped.
        '''
        
        measFluxes = read_measurements_from_file(file)
            
        for emuid, [mean, sd] in measFluxes.iterrows():    
            self.model.measured_fluxes[emuid] = [mean, sd]
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_measured_fluxes, measFluxes.index.tolist()))
    
    
    def _unset_measured_fluxes(self, fluxids):
        '''
        Parameters
        ----------
        fluxids: str or list of str
            Measured flux ID(s).
        '''
        
        if not isinstance(fluxids, Iterable):
            fluxids = [fluxids]
            
        for fluxid in fluxids: 
            if fluxid in self.model.measured_fluxes:
                self.model.measured_fluxes.pop(fluxid)
            
        
    def set_unbalanced_metabolites(self, metabids):
        '''
        Parameters
        ----------
        metabids: str or list of str
            Unbalanced metabolite ID(s).
        '''
        
        if not isinstance(metabids, Iterable):
            metabids = [metabids]
        
        self.model.unbalanced_metabolites.update(metabids)
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_unbalanced_metabolites, metabids))
        
        
    def _unset_unbalanced_metabolites(self, metabids):
        '''
        Parameters
        ----------
        metabids: list of str
            Unbalanced metabolite IDs.
        '''
        
        for metabid in metabids:
            if metabid in self.model.unbalanced_metabolites:
                self.model.unbalanced_metabolites.remove(metabid)
            
    
    def set_flux_bounds(self, fluxid, bounds):
        '''
        Set lower and upper bounds of flux.
        
        Parameters
        ----------
        fluxid: str or 'all'
            Flux ID, i.e., reaction ID. Since forward and backward fluxes of 
            reversible reaction are usually unknown, the method is used to set 
            the range of net fluxes.
            
            If 'all', all fluxes will be set to the range.
        bounds: 2-list
            [lower bound, upper bound]. Lower bound is not allow to equal upper bound.
            Use set_measured_flux (or set_measured_fluxes_from_file) to set fixed value
            of flux.

            For irreversible reaction, the lower bound will be set to zero ignorant 
            of bounds[0].
        '''
        
        fluxids = []
        if bounds[0] < bounds[1]:   # lower bound not allow to equal upper bound
            if fluxid == 'all':
                for rxnid in self.model.reactions:
                    self._set_flux_bounds(rxnid, bounds)
                    fluxids.append(rxnid)    
            elif fluxid in self.model.reactions:
                self._set_flux_bounds(fluxid, bounds)
                fluxids = [fluxid]
            else:
                raise ValueError('flux range set to nonexistent reaction %s' % fluxid)
        else:
            raise ValueError('lower bound of flux should be less than upper bound, \nuse set_measured_flux or set_measured_fluxes_from_file to set fixed flux value')
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_flux_bounds, fluxids))
    
    
    def _decompose_network(self, n_jobs):
        '''
        Parameters
        ----------
        n_jobs: int
            # of jobs to run in parallel.
        '''
        
        if not self.model.measured_MDVs:
            raise ValueError('call set_measured_MDV or set_measured_MDVs_from_file first')
        
        if not self.model.EAMs:
            if n_jobs <= 0:
                raise ValueError('n_jobs should be a positive value')    
            else:
                self.model.target_EMUs = list(self.model.measured_MDVs.keys())
                
                metabids = []
                atom_nos = []
                for emuid in self.model.target_EMUs:
                    metabid, atomNOs = emuid.split('_')
                    metabids.append(metabid)
                    atom_nos.append(atomNOs)
                
                EAMs = self.model._decompose_network(metabids, atom_nos, n_jobs = n_jobs)
                for size, EAM in EAMs.items():
                    self.model.EAMs[size] = EAM
                
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_decomposition)
            
            
    def _unset_decomposition(self):
        
        self.model.EAMs.clear()        
    
    
    def _calculate_matrix_As_and_Bs_derivatives_p(self, kind, n_jobs):
        '''
        Parameters
        ----------
        kind: {"ss", "inst"}
            * If "ss", variables are free fluxes only.
            * If "inst", variables include both free fluxes and concentrations.
        '''
        
        if not self.model.matrix_As_der_p or not self.model.matrix_Bs_der_p:
            self.calculator._calculate_matrix_As_and_Bs_derivatives_p(kind, n_jobs)
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_matrix_As_and_Bs_derivatives_p)
        
        
    def _unset_matrix_As_and_Bs_derivatives_p(self):
        
        self.model.matrix_As_der_p.clear()
        self.model.matrix_Bs_der_p.clear()
    
    
    def _calculate_substrate_MDV_derivatives_p(self, kind, extra_subs = None):
        '''
        Parameters
        ----------
        kind: {"ss", "inst"}
            * If "ss", variables are free fluxes only.
            * If "inst", variables include free fluxes and concentrations.
        extra_subs: str or list of str
            Metabolite ID(s), additional metabolites considered as substrates.    
        '''
        
        if extra_subs is not None and not isinstance(extra_subs, list):
            extra_subs = [extra_subs]

        if not self.model.substrate_MDVs_der_p:
            self.calculator._calculate_substrate_MDV_derivatives_p(kind, extra_subs)
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_substrate_MDV_derivatives_p)
            
    
    def _unset_substrate_MDV_derivatives_p(self):
        
        self.model.substrate_MDVs_der_p.clear()
        
    
    def _calculate_null_space(self):
        
        if self.model.null_space is None:
            self.calculator._calculate_null_space()
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_null_space)
            
            
    def _unset_null_space(self):
        
        self.model.null_space = None
            
    
    def _calculate_transform_matrix(self):
        
        if self.model.transform_matrix is None:
            self.calculator._calculate_transform_matrix()
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_transform_matrix)    
    
            
    def _unset_transform_matrix(self):
        
        self.model.transform_matrix = None
    
    
    def _calculate_measured_MDVs_inversed_covariance_matrix(self):
    
        if self.model.measured_MDVs_inv_cov is None:
            self.calculator._calculate_measured_MDVs_inversed_covariance_matrix()
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_measured_MDVs_inversed_covariance_matrix)    
    
    
    def _unset_measured_MDVs_inversed_covariance_matrix(self):
        
        self.model.measured_MDVs_inv_cov = None
        
    
    def _calculate_measured_fluxes_inversed_covariance_matrix(self):
        
        if self.model.measured_fluxes_inv_cov is None:
            self.calculator._calculate_measured_fluxes_inversed_covariance_matrix()
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_measured_fluxes_inversed_covariance_matrix)
    
            
    def _unset_measured_fluxes_inversed_covariance_matrix(self):
        
        self.model.measured_fluxes_inv_cov = None
        
    
    def _calculate_measured_fluxes_derivative_p(self, kind):
        '''
        Parameters
        ----------
        kind: {"ss", "inst"}
            * If "ss", variables are free fluxes only.
            * If "inst", variables include free fluxes and concentrations.
        '''
        
        if self.model.measured_fluxes_der_p is None:
            self.calculator._calculate_measured_fluxes_derivative_p(kind)
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_measured_fluxes_derivative_p)
        
            
    def _unset_measured_fluxes_derivative_p(self):
        
        self.model.measured_fluxes_der_p = None
    
    
    def _estimate_fluxes_range(self, exclude_metabs = None):
        '''
        Parameters
        ----------
        exclude_metabs: list
            Metabolite IDs, metabolites excluded from mass balance.
        '''
        
        fluxids = []
        if not self.model.net_fluxes_range:
            FVAres = self.estimate_fluxes_range(exclude_metabs = exclude_metabs, 
                                                show_progress = False)
            for fluxid, fluxRange in FVAres.flux_ranges.items():
                self.model.net_fluxes_range[fluxid] = fluxRange
                fluxids.append(fluxid)
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_net_fluxes_range, fluxids))
            
            
    def _unset_net_fluxes_range(self, fluxids):
        '''
        Parameters
        ----------
        fluxids: str or list of str
            Fluxe ID(s).
        '''
        
        if not isinstance(fluxids, Iterable):
            fluxids = [fluxids]
        
        for fluxid in fluxids:
            if fluxid in self.model.net_fluxes_range:
                self.model.net_fluxes_range.pop(fluxid)
        
    
    def prepare(self, dilution_from = None, n_jobs = 1):
        '''
        Parameters
        ----------
        dilution_from: str or list of str
            ID(s) of unlabeled (inactive) metabolite leading to dilution effect.
            These metabolites have zero stoichiometric coefficients in reaction network.
        n_jobs: int
            If n_jobs > 1, preparation will run in parallel.
        '''
        
        self._decompose_network(n_jobs)
        self._calculate_null_space()
        self._calculate_transform_matrix()
        self._lambdify_matrix_As_and_Bs()
        self._calculate_matrix_As_and_Bs_derivatives_p('ss', n_jobs)
        self._calculate_substrate_MDVs(dilution_from)
        self._calculate_substrate_MDV_derivatives_p('ss', dilution_from)
        self._calculate_measured_MDVs_inversed_covariance_matrix()
        self._calculate_measured_fluxes_inversed_covariance_matrix()
        self._calculate_measured_fluxes_derivative_p('ss')
        self._set_default_flux_bounds()
        self._estimate_fluxes_range(self.model.unbalanced_metabolites)
        
    
    def _check_dependencies(self, fit_measured_fluxes):
        '''
        Parameters
        ----------
        fit_measured_fluxes: bool
            Whether to fit measured fluxes.
        '''

        if not self.model.net_fluxes_bounds:
            raise ValueError('call set_flux_bounds first')
        if not self.model.measured_MDVs:
            raise ValueError('call set_measured_MDV or set_measured_MDVs_from_file first')
        if not self.model.measured_fluxes:
            raise ValueError('call set_measured_flux or set_measured_fluxes_from_file first')
        if not self.model.labeling_strategy:
            raise ValueError('call labeling_strategy first')
        
        checklist = [not self.model.target_EMUs, 
                     self.model.transform_matrix is None, 
                     self.model.null_space is None,
                     self.model.measured_MDVs_inv_cov is None, 
                     not self.model.matrix_As, 
                     not self.model.matrix_Bs,
                     not self.model.matrix_As_der_p, 
                     not self.model.matrix_Bs_der_p, 
                     not self.model.substrate_MDVs,
                     not self.model.substrate_MDVs_der_p, 
                     self.model.measured_fluxes_der_p is None]
        if fit_measured_fluxes:
            checklist.append(self.model.measured_fluxes_inv_cov is None)
        
        if any(checklist):
            raise ValueError('call prepare first')


    def solve(self, fit_measured_fluxes = True, ini_fluxes = None, solver = 'slsqp', 
              tol = 1e-6, max_iters = 400, show_progress = True):
        '''
        Parameters
        ----------
        fit_measured_fluxes: bool
            Whether to fit measured fluxes.
        ini_fluxes: ser or file in .tsv or .xlsx
            Initial values of net fluxes
        solver: {"slsqp", "ralg"}
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
            
        optModel = MFAModel(self.model, fit_measured_fluxes, solver)
        optModel.build_objective()
        optModel.build_gradient()
        optModel.build_flux_bound_constraints()
        optModel.build_initial_flux_values(ini_netfluxes = iniFluxes)
        
        with Progress('fitting', silent = not show_progress):
            res = optModel.solve_flux(tol, max_iters)

        return FitResults(*res[:7], deepcopy(res[7]), res[8], deepcopy(res[9]), *res[10:])
    
    
    def _solve_with_confidence_intervals(self, fit_measured_fluxes, ini_fluxes, solver, 
                                         tol, max_iters, nruns):
        '''
        Parameters
        ----------
        fit_measured_fluxes: bool
            Whether to fit measured fluxes.
        ini_fluxes: ser or file in .tsv or .xlsx or None
            Initial values of net fluxes.   
        solver: {"slsqp", "ralg"}
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
        
        if ini_fluxes is not None:
            iniFluxes = read_initial_values(ini_fluxes, self.model.netfluxids)
        else:
            iniFluxes = ini_fluxes

        optTotalfluxesSet = []
        optNetfluxesSet = []
        for _ in range(nruns):
            self.calculator._generate_random_fluxes()
            self.calculator._generate_random_MDVs()
            
            optModel = MFAModel(self.model, fit_measured_fluxes, solver)
            optModel.build_objective()
            optModel.build_gradient()
            optModel.build_flux_bound_constraints()
            optModel.build_initial_flux_values(ini_netfluxes = iniFluxes)

            while True:
                optTotalfluxes, optNetfluxes, *_, isSuccess = optModel.solve_flux(tol, max_iters)
                if isSuccess:
                    break
            optTotalfluxesSet.append(optTotalfluxes)
            optNetfluxesSet.append(optNetfluxes)
            
            self.calculator._reset_measured_fluxes()
            self.calculator._reset_measured_MDVs()
               
        return optTotalfluxesSet, optNetfluxesSet
        
    
    def solve_with_confidence_intervals(self, fit_measured_fluxes = True, ini_fluxes = None, 
                                        solver = 'slsqp', tol = 1e-6, max_iters = 400, 
                                        n_runs = 100, n_jobs = 1, show_progress = True):
        '''
        Parameters
        ----------
        fit_measured_fluxes: bool
            Whether to fit measured fluxes.
        ini_fluxes: ser or file in .tsv or .xlsx
            Initial values of net fluxes.  
        solver: {"slsqp", "ralg"}
            * If "slsqp", scipy.optimize.minimze will be used.
            * If "ralg", openopt NLP solver will be used.
        tol: float
            Tolerance for termination.
        max_iters: int
            Maximum # of iterations.
        show_progress: bool
            Whether to show the progress bar.
        n_runs: int
            # of runs to estimate confidence intervals.
        n_jobs: int
            # of jobs to run in parallel.
        '''
        
        self._check_dependencies(fit_measured_fluxes)
        
        self._unset_matrix_As_and_Bs()

        if n_runs <= n_jobs:
            nruns_worker = 1
        else:
            nruns_worker = ceil(n_runs/n_jobs)
        
        pool = Pool(processes = n_jobs)
        with Progress('fitting with CIs', silent = not show_progress):
            fluxesSet = []
            for _ in range(n_jobs):
                fluxes = pool.apply_async(func = self._solve_with_confidence_intervals, 
                                          args = (fit_measured_fluxes,
                                                  ini_fluxes,
                                                  solver,
                                                  tol,
                                                  max_iters,
                                                  nruns_worker))
                fluxesSet.append(fluxes)
            
            pool.close()    
            pool.join()
        
        fluxesSet = [fluxes.get() for fluxes in fluxesSet]
        
        totalFluxesSet = []
        netFluxesSet = []
        for totalFluxesSubset, netFluxesSubset in fluxesSet:
            totalFluxesSet.extend(totalFluxesSubset)
            netFluxesSet.extend(netFluxesSubset)
        
        return FitMCResults(totalFluxesSet, netFluxesSet)
