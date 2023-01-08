#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '05/19/2022'




import numpy as np
import pandas as pd
from scipy.linalg import pinv2
from scipy.optimize import minimize
try:
    from openopt import NLP
except ModuleNotFoundError:
    OPENOPT_INSTALLED = False
else:
    OPENOPT_INSTALLED = True
from ..utils.utils import Calculator




class MFAModel():
    '''
    Parameters
    ----------
    model: Model
        Freeflux Model.
    fit_measured_fluxes: bool
        Whether to fit measured fluxes.
    solvor: {"slsqp", "ralg"}
        * If "slsqp", scipy.optimize.minimze will be used.
        * If "ralg", openopt NLP solver will be used.
    '''
    
    def __init__(self, model, fit_measured_fluxes, solver = 'slsqp'):
        '''
        Parameters
        ----------
        model: Model
            Freeflux Model.
        fit_measured_fluxes: bool
            Whether to fit measured fluxes.
        solvor: {"slsqp", "ralg"}
            * If "slsqp", scipy.optimize.minimze will be used.
            * If "ralg", openopt NLP solver will be used.
        '''
        
        self.model = model
        self.calculator = Calculator(self.model)
        self.fit_measured_fluxes = fit_measured_fluxes
        self.solver = solver
        
        self.N = self.model.null_space
        self.T = self.model.transform_matrix
        
        self.ntotalfluxes = len(self.model.totalfluxids)
        
    
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
        
        simMDVs, simMDVsDer = self.calculator._calculate_MDVs_and_derivatives_p()
        expMDVs = self.model.measured_MDVs
        diff = np.concatenate([simMDVs[emuid] - expMDVs[emuid][0] for emuid in self.model.target_EMUs])
        
        dxsim_dp = np.vstack([simMDVsDer[emuid] for emuid in self.model.target_EMUs])
        
        return diff, dxsim_dp
        
    
    def _calculate_sim_fluxes_derivative(self):
        
        dvsim_dp = self.model.measured_fluxes_der_p
        
        return dvsim_dp


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
        
        A1 = self.N
        A2 = self.T@self.N
        A3 = -A2
        
        b1 = np.zeros(self.ntotalfluxes)
        vnet_lb, vnet_ub = np.array(list(self.model.net_fluxes_range.values())).T
        b2 = vnet_lb
        b3 = -vnet_ub
        
        A = np.vstack((A1, A2, A3))
        b = np.concatenate((b1, b2, b3))

        if self.solver == 'slsqp':
            self.constrs = {'type': 'ineq', 'fun': lambda u: A@u - b}
        elif self.solver == 'ralg':
            self.A = -A
            self.b = -b
        
        
    def build_initial_flux_values(self, ini_netfluxes = None):
        '''
        Parameters
        ----------
        ini_netfluxes: array
            Initial guess of net fluxes.
        '''
        
        if ini_netfluxes is None:
            vnet_lb, vnet_ub = np.array(list(self.model.net_fluxes_range.values())).T
            vnet_ini = np.random.uniform(low = vnet_lb, high = vnet_ub)
        else:
            vnet_ini = ini_netfluxes
            
        u_ini = pinv2(self.T@self.N)@vnet_ini
        
        self.x0 = u_ini
        
        
    def _initialize_total_fluxes(self):
        
        for fluxid in self.model.totalfluxids:
            self.model.total_fluxes[fluxid] = 0.0
        

    def _solve_flux_slsqp(self, tol, max_iters, disp):
        
        res = minimize(fun = self.f,
                       x0 = self.x0,
                       method = 'SLSQP',
                       jac = self.df,
                       constraints = self.constrs,
                       options = {'ftol': tol,
                                  'maxiter': max_iters, 
                                  'disp': disp})
       
        return res.fun, res.x, res.status in [0, 2]


    def _solve_flux_ralg(self, tol, max_iters, disp):

        if OPENOPT_INSTALLED:
            res = res = NLP(f = self.f, 
                            x0 = self.x0, 
                            df = self.df, 
                            A = self.A, 
                            b = self.b, 
                            xtol = tol, 
                            ftol = tol, 
                            maxIter = max_iters, 
                            iprint = 1 if disp else -1).solve('ralg')

            return res.ff, res.xf, res.istop > 0 or res.istop == -7
        
        else:
            raise ModuleNotFoundError('install openopt first')


    def _calculate_residuals(self):
        
        MDV_diff = self._calculate_difference_sim_exp_MDVs()
        MDV_diag = np.diag(self.model.measured_MDVs_inv_cov)**0.5
        resids1 = MDV_diff*MDV_diag
        
        if self.fit_measured_fluxes:    
            flux_diff = self._calculate_difference_sim_exp_fluxes()
            flux_diag = np.diag(self.model.measured_fluxes_inv_cov)**0.5
            resids2 = flux_diff*flux_diag
            resids = np.concatenate((resids1, resids2))
        
        return resids         


    def _get_exp_and_sim_MDVs(self):

        exp_MDVs = self.model.measured_MDVs
        
        sim_MDVs_all = self.calculator._calculate_MDVs()
        sim_MDVs = {emuid: sim_MDVs_all[emuid] for emuid in self.model.target_EMUs}

        return exp_MDVs, sim_MDVs


    def _get_exp_and_sim_fluxes(self):

        if self.fit_measured_fluxes:
            exp_fluxes = self.model.measured_fluxes
            sim_fluxes = self.model.total_fluxes[self.model.measured_fluxes.keys()].to_dict()
        else:
            exp_fluxes = {}
            sim_fluxes = {}            

        return exp_fluxes, sim_fluxes


    def _get_nmeasurements(self, opt_resids):

        return opt_resids.size

         
    def _get_nparameters(self, opt_p):

        return opt_p.size


    def _get_hessian(self, opt_p):

        self.build_hessian()
        hess = self.ddf(opt_p)

        return hess


    def solve_flux(self, tol = 1e-6, max_iters = 400, disp = False):    
        
        self._initialize_total_fluxes()
        
        if self.solver == 'slsqp':
            opt_obj, opt_u, is_success = self._solve_flux_slsqp(tol, max_iters, disp)
        elif self.solver == 'ralg':
            opt_obj, opt_u, is_success = self._solve_flux_ralg(tol, max_iters, disp)
        else:
            raise ValueError('currently only "slsqp" and "ralg" are acceptable')
        
        opt_totalfluxes = pd.Series(self.N@opt_u, index = self.model.totalfluxids)
        opt_netfluxes = pd.Series(self.T@self.N@opt_u, index = self.model.netfluxids)
        
        opt_resids = self._calculate_residuals()
        
        exp_MDVs, sim_MDVs = self._get_exp_and_sim_MDVs()
        exp_fluxes, sim_fluxes = self._get_exp_and_sim_fluxes()
        
        nmeas = self._get_nmeasurements(opt_resids)
        nparams = self._get_nparameters(opt_u)
        
        hess = self._get_hessian(opt_u)

        _, dxsim_du = self._calculate_sim_MDVs_derivative()
        dvsim_du = self._calculate_sim_fluxes_derivative()
        
        return (opt_totalfluxes, opt_netfluxes, opt_obj, opt_resids, nmeas, nparams, sim_MDVs, exp_MDVs, 
                sim_fluxes, exp_fluxes, hess, self.N, self.T, dxsim_du, dvsim_du, 
                self.model.measured_MDVs_inv_cov, self.model.measured_fluxes_inv_cov, is_success)
        
        
    
        
class InstMFAModel(MFAModel):

    def __init__(self, *args):
        
        super().__init__(*args)
        
        self.nfreefluxes = self.N.shape[1]
        self.nconcs = len(self.model.concids)
        self.nnetfluxes = len(self.model.netfluxids)
        
    
    def _calculate_difference_sim_exp_MDVs(self):
        
        simMDVs = self.calculator._calculate_inst_MDVs()
        expMDVs = self.model.measured_inst_MDVs
        diff = np.concatenate([simMDVs[emuid][t] - expMDVs[emuid][t][0] for emuid in self.model.target_EMUs
                                                                        for t in expMDVs[emuid] if t != 0])
        
        return diff
        
        
    def _calculate_sim_MDVs_derivative(self):
        
        simMDVs, simMDVsDer = self.calculator._calculate_inst_MDVs_and_derivatives_p()
        expMDVs = self.model.measured_inst_MDVs
        diff = np.concatenate([simMDVs[emuid][t] - expMDVs[emuid][t][0] for emuid in self.model.target_EMUs
                                                                        for t in expMDVs[emuid] if t != 0])
        
        dxsim_dp = np.vstack([simMDVsDer[emuid][t] for emuid in self.model.target_EMUs
                                                   for t in expMDVs[emuid] if t != 0])
        
        return diff, dxsim_dp
        
    
    def build_objective(self):
        
        def _f(p):
            u, c = p[:self.nfreefluxes], p[self.nfreefluxes:]
            self.model.total_fluxes[:] = self.N@u
            self.model.concentrations[:] = c
            
            MDV_diff = self._calculate_difference_sim_exp_MDVs()
            MDV_inv_cov = self.model.measured_inst_MDVs_inv_cov
            
            obj = MDV_diff@MDV_inv_cov@MDV_diff
            
            return obj
        
        def f1(p):
            return _f(p)
            
        def f2(p):
            obj1 = _f(p)
            
            flux_diff = self._calculate_difference_sim_exp_fluxes()
            flux_inv_cov = self.model.measured_fluxes_inv_cov
            obj2 = flux_diff@flux_inv_cov@flux_diff
            
            return obj1 + obj2
            
        self.f = f2 if self.fit_measured_fluxes else f1
    
        
    def build_gradient(self):
        
        def _df(p):
            u, c = p[:self.nfreefluxes], p[self.nfreefluxes:]
            self.model.total_fluxes[:] = self.N@u
            self.model.concentrations[:] = c
            
            MDV_diff, MDV_der = self._calculate_sim_MDVs_derivative()
            MDV_inv_cov = self.model.measured_inst_MDVs_inv_cov
            grad = MDV_der.T@MDV_inv_cov@MDV_diff
            
            return grad
            
        def df1(p):
            return _df(p)
            
        def df2(p):
            grad1 = _df(p)
            
            flux_der = self._calculate_sim_fluxes_derivative()
            flux_diff = self._calculate_difference_sim_exp_fluxes()
            flux_inv_cov = self.model.measured_fluxes_inv_cov
            grad2 = flux_der.T@flux_inv_cov@flux_diff
            
            return grad1 + grad2
        
        self.df = df2 if self.fit_measured_fluxes else df1
        
    
    def build_hessian(self):
        
        def _ddf(p):
            u, c = p[:self.nfreefluxes], p[self.nfreefluxes:]
            self.model.total_fluxes[:] = self.N@u
            self.model.concentrations[:] = c
            
            _, MDV_der = self._calculate_sim_MDVs_derivative()
            MDV_inv_cov = self.model.measured_inst_MDVs_inv_cov
            hess = MDV_der.T@MDV_inv_cov@MDV_der
            
            return hess
            
        def ddf1(p):
            return _ddf(p)
        
        def ddf2(p):
            hess1 = _ddf(p)
            
            flux_der = self._calculate_sim_fluxes_derivative()
            flux_inv_cov = self.model.measured_fluxes_inv_cov
            hess2 = flux_der.T@flux_inv_cov@flux_der
            
            return hess1 + hess2
        
        self.ddf = ddf2 if self.fit_measured_fluxes else ddf1
    
    
    def build_flux_and_conc_bound_constraints(self):
        
        A1 = self.N
        A2 = self.T@self.N
        A3 = -A2
        A4 = np.eye(self.nconcs)
        
        b1 = np.zeros(self.ntotalfluxes)
        vnet_lb, vnet_ub = np.array(list(self.model.net_fluxes_range.values())).T
        b2 = vnet_lb
        b3 = -vnet_ub
        b4 = np.zeros(self.nconcs)
        
        A = np.zeros((self.ntotalfluxes+2*self.nnetfluxes+self.nconcs, self.nfreefluxes+self.nconcs))
        A[:(self.ntotalfluxes+2*self.nnetfluxes), :self.nfreefluxes] = np.vstack((A1, A2, A3))
        A[(self.ntotalfluxes+2*self.nnetfluxes):, self.nfreefluxes:] = A4
        b = np.concatenate((b1, b2, b3, b4))
        
        if self.solver == 'slsqp':
            self.constrs = {'type': 'ineq', 'fun': lambda p: A@p - b}
        elif self.solver == 'ralg':
            self.A = -A
            self.b = -b
    
    
    def build_initial_flux_and_conc_values(self, ini_netfluxes = None, ini_concs = None):
        '''
        Parameters
        ----------
        ini_netfluxes: array
            Initial guess of net fluxes.
        ini_concs: array
            Initial guess of concentrations.
        '''
        
        if ini_netfluxes is None:
            vnet_lb, vnet_ub = np.array(list(self.model.net_fluxes_range.values())).T
            vnet_ini = np.random.uniform(low = vnet_lb, high = vnet_ub)
        else:
            vnet_ini = ini_netfluxes
        u_ini = pinv2(self.T@self.N)@vnet_ini
        
        if ini_concs is None:
            c_lb, c_ub = np.array(list(self.model.concentrations_range.values())).T
            c_ini = np.random.uniform(low = c_lb, high = c_ub)
        else:
            c_ini = ini_concs

        self.x0 = np.concatenate((u_ini, c_ini))
    
    
    def _initialize_total_fluxes_and_concs(self):
        
        for fluxid in self.model.totalfluxids:
            self.model.total_fluxes[fluxid] = 0.0
    
        for metabid in self.model.concids:
            self.model.concentrations[metabid] = 0.0
            
    
    def _calculate_residuals(self):
        
        MDV_diff = self._calculate_difference_sim_exp_MDVs()
        MDV_diag = np.diag(self.model.measured_inst_MDVs_inv_cov)**0.5
        resids1 = MDV_diff*MDV_diag
        
        if self.fit_measured_fluxes:    
            flux_diff = self._calculate_difference_sim_exp_fluxes()
            flux_diag = np.diag(self.model.measured_fluxes_inv_cov)**0.5
            resids2 = flux_diff*flux_diag
            resids = np.concatenate((resids1, resids2))
        
        return resids
        

    def _get_exp_and_sim_inst_MDVs(self):

        exp_inst_MDVs = self.model.measured_inst_MDVs

        sim_inst_MDVs_all = self.calculator._calculate_inst_MDVs()
        self.calculator._build_initial_sim_MDVs()
        sim_inst_MDVs = {}
        for emuid in self.model.target_EMUs:
            instMDVs = self.model.initial_sim_MDVs[emuid].copy()
            instMDVs.update({t: sim_inst_MDVs_all[emuid][t] for t in sim_inst_MDVs_all[emuid]})
            sim_inst_MDVs[emuid] = instMDVs

        return exp_inst_MDVs, sim_inst_MDVs 

    
    def solve_flux(self, tol = 1e-6, max_iters = 400, disp = False):
        
        self._initialize_total_fluxes_and_concs()

        if self.solver == 'slsqp':
            opt_obj, opt_p, is_success = self._solve_flux_slsqp(tol, max_iters, disp)
        elif self.solver == 'ralg':
            opt_obj, opt_p, is_success = self._solve_flux_ralg(tol, max_iters, disp)
        else:
            raise ValueError('currently only "slsqp" and "ralg" are acceptable')    

        opt_u = opt_p[:self.nfreefluxes]
        opt_c = opt_p[self.nfreefluxes:]
        
        opt_totalfluxes = pd.Series(self.N@opt_u, index = self.model.totalfluxids)
        opt_netfluxes = pd.Series(self.T@self.N@opt_u, index = self.model.netfluxids)        
        opt_concs = pd.Series(opt_c, index = self.model.concids)
        
        opt_resids = self._calculate_residuals()
        
        exp_inst_MDVs, sim_inst_MDVs = self._get_exp_and_sim_inst_MDVs()
        exp_fluxes, sim_fluxes = self._get_exp_and_sim_fluxes()
        
        nmeas = self._get_nmeasurements(opt_resids)
        nparams = self._get_nparameters(opt_p)
        
        hess = self._get_hessian(opt_p)

        _, dxsim_dp = self._calculate_sim_MDVs_derivative()
        dxsim_du = dxsim_dp[:, :self.nfreefluxes]
        
        dvsim_dp = self._calculate_sim_fluxes_derivative()
        dvsim_du = dvsim_dp[:, :self.nfreefluxes]
        
        return (opt_totalfluxes, opt_netfluxes, opt_concs, opt_obj, opt_resids, nmeas, nparams, sim_inst_MDVs, 
                exp_inst_MDVs, sim_fluxes, exp_fluxes, hess, self.N, self.T, dxsim_du, dvsim_du, 
                self.model.measured_inst_MDVs_inv_cov, self.model.measured_fluxes_inv_cov, is_success)
        