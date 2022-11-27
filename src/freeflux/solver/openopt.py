#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '05/19/2022'




import numpy as np
import pandas as pd
from scipy.linalg import pinv2
from openopt import NLP
from ..utils.utils import Calculator




class OoptModel():
    
    def __init__(self, model, fit_measured_fluxes):
        '''
        Parameters
        model: Model
        fit_measured_fluxes: bool
            whether to fit measured fluxes
        '''
        
        self.model = model
        self.calculator = Calculator(self.model)
        self.fit_measured_fluxes = fit_measured_fluxes
        
        self.N = self.model.null_space
        self.T = self.model.transform_matrix
        
    
    def _calculate_difference_sim_exp_MDVs(self):
        
        simMDVs = self.calculator._calculate_MDVs()
        expMDVs = self.model.measured_MDVs
        diff = np.concatenate([simMDVs[emuid] - expMDVs[emuid][0] for emuid in self.model.target_EMUs])
        
        return diff
            
            
    def _calculate_difference_sim_exp_fluxes(self):
        
        simFluxes = self.model.total_fluxes[self.model.measured_fluxes.keys()]
        expFluxes = np.array([mean for mean, _ in self.model.measured_fluxes.values()])
        diff = simFluxes - expFluxes
        
        return diff
        
        
    def _calculate_sim_MDVs_derivative(self):
        
        simMDVs, simMDVsDer = self.calculator._calculate_MDVs_and_derivatives()
        
        expMDVs = self.model.measured_MDVs
        diff = np.concatenate([simMDVs[emuid] - expMDVs[emuid][0] for emuid in self.model.target_EMUs])
        
        dxsimdv = np.vstack([simMDVsDer[emuid] for emuid in self.model.target_EMUs])
        dxsimdu = dxsimdv@self.N
        
        return diff, dxsimdu
        
        
    def _calculate_sim_fluxes_derivative(self):
        
        dvsimdv = self.model.measured_fluxes_der
        dvsimdu = dvsimdv@self.N
        
        return dvsimdu
        
        
    def build_objective(self):
        
        def _f(u):
            self.model.total_fluxes[:] = self.N@u
            
            MDV_diff = self._calculate_difference_sim_exp_MDVs()
            MDV_inv_cov = self.model.measured_MDVs_inv_cov
            obj = MDV_diff@MDV_inv_cov@MDV_diff
            
            return obj
        
        def f1(u):
            return _f(u)
            
        def f2(u):
            obj1 = _f(u)
            
            flux_diff = self._calculate_difference_sim_exp_fluxes()
            flux_inv_cov = self.model.measured_fluxes_inv_cov
            obj2 = flux_diff@flux_inv_cov@flux_diff
            
            return obj1 + obj2
            
        self.f = f2 if self.fit_measured_fluxes else f1
        
    
    def build_gradient(self):
        
        def _df(u):
            self.model.total_fluxes[:] = self.N@u
            
            MDV_diff, MDV_der = self._calculate_sim_MDVs_derivative()
            MDV_inv_cov = self.model.measured_MDVs_inv_cov
            grad = MDV_der.T@MDV_inv_cov@MDV_diff
            
            return grad
            
        def df1(u):
            return _df(u)
            
        def df2(u):
            grad1 = _df(u)
            
            flux_der = self._calculate_sim_fluxes_derivative()
            flux_diff = self._calculate_difference_sim_exp_fluxes()
            flux_inv_cov = self.model.measured_fluxes_inv_cov
            grad2 = flux_der.T@flux_inv_cov@flux_diff
            
            return grad1 + grad2
        
        self.df = df2 if self.fit_measured_fluxes else df1

    
    def build_hessian(self):
        
        def _ddf(u):
            self.model.total_fluxes[:] = self.N@u
            
            _, MDV_der = self._calculate_sim_MDVs_derivative()
            MDV_inv_cov = self.model.measured_MDVs_inv_cov
            hess = MDV_der.T@MDV_inv_cov@MDV_der
            
            return hess
            
        def ddf1(u):
            return _ddf(u)
        
        def ddf2(u):
            hess1 = _ddf(u)
            
            flux_der = self._calculate_sim_fluxes_derivative()
            flux_inv_cov = self.model.measured_fluxes_inv_cov
            hess2 = flux_der.T@flux_inv_cov@flux_der
            
            return hess1 + hess2
        
        self.ddf = ddf2 if self.fit_measured_fluxes else ddf1    
    
            
    def build_flux_bound_constraints(self):
        
        A1 = -self.N
        A2 = -self.T@self.N
        A3 = -A2
        
        b1 = np.zeros(len(self.calculator.totalfluxids))
        vnet_lb, vnet_ub = np.array(list(self.model.net_fluxes_range.values())).T
        b2 = -vnet_lb
        b3 = vnet_ub
        
        self.A = np.vstack((A1, A2, A3))
        self.b = np.concatenate((b1, b2, b3))
    
        
    def build_initial_flux_values(self):
        
        vnet_lb, vnet_ub = np.array(list(self.model.net_fluxes_range.values())).T
        vnet_ini = np.random.uniform(low = vnet_lb, high = vnet_ub)
        u_ini = pinv2(self.T@self.N)@vnet_ini
        
        self.x0 = u_ini
        
        
    def solve_flux(self, tol = 1e-6, max_iters = 400, iprint = -1):
        
        problem = NLP(f = self.f, 
                      x0 = self.x0, 
                      df = self.df, 
                      A = self.A, 
                      b = self.b, 
                      xtol = tol, 
                      ftol = tol, 
                      maxIter = max_iters, 
                      iprint = iprint)
        res = problem.solve('ralg')
        
        opt_obj = res.ff
        opt_u = res.xf
        
        opt_totalfluxes = pd.Series(self.N@opt_u, index = self.calculator.totalfluxids)
        opt_netfluxes = pd.Series(self.T@self.N@opt_u, index = self.calculator.netfluxids)
        
        opt_resids = self._calculate_residuals()
        
        sim_MDVs_all = self.calculator._calculate_MDVs()
        sim_MDVs = {emuid: sim_MDVs_all[emuid] for emuid in self.model.target_EMUs}
        exp_MDVs = self.model.measured_MDVs
        
        if self.fit_measured_fluxes:
            sim_fluxes = self.model.total_fluxes[self.model.measured_fluxes.keys()].to_dict()
            exp_fluxes = self.model.measured_fluxes
        else:
            sim_fluxes = {}
            exp_fluxes = {}
        
        nmeas = opt_resids.size - len(self.model.target_EMUs)
        nparams = opt_u.size
        
        self.build_hessian()
        hess = self.ddf(opt_u)
        
        return (opt_obj, opt_resids, opt_totalfluxes, opt_netfluxes, nmeas, nparams, sim_MDVs, exp_MDVs, 
                sim_fluxes, exp_fluxes, hess, self.N)
        
        
    def _calculate_residuals(self):
        
        MDV_diff = self._calculate_difference_sim_exp_MDVs()
        MDV_diag = np.diag(self.model.measured_MDVs_inv_cov)**0.5
        resids = MDV_diff*MDV_diag
        
        if self.fit_measured_fluxes:    
            flux_diff = self._calculate_difference_sim_exp_fluxes()
            flux_diag = np.diag(self.model.measured_fluxes_inv_cov)**0.5
            resids = np.concatenate((resids, flux_diff*flux_diag))
        
        return resids
        
        
        
        
        
        
        
        
        