#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '04/14/2022'




import numpy as np
import pandas as pd
from scipy.linalg import pinv2
from ..core.mdv import MDV
from ..analysis.stats import (_chi2_test, _normal_probability, _simulated_vs_measured_MDVs, 
                              _simulated_vs_measured_fluxes, _simulated_vs_measured_inst_MDVs,
                              _confidence_intervals_le, _confidence_intervals_mc, _MDV_kinetics,
                              _contribution_matrix, _sensitivity)




class pDict(dict):

    def __init__(self, *args, digits = 3, **kwargs):
        
        super().__init__(*args, **kwargs)
        self.digits = digits
    
        
    def __repr__(self):
        
        itemStrs = []
        for k, v in self.items():
            if isinstance(v, MDV):
                itemStr = '%s: %s' % (k, v)
            elif isinstance(v, list):
                itemStr = '%s: [%s]' % (k, ', '.join(map(str, np.round(v, self.digits))))
            else:
                itemStr = '%s: %s' % (k, np.round(v, self.digits))
            itemStrs.append(itemStr)
        
        return '\n'.join(itemStrs)

    
    
    
class FBAResults():
    '''
    Parameters
    ----------
    obj: dict
        Reaction ID => coefficient, i.e., the objective function.
    opt_obj: float
        Optimal objective.
    opt_fluxes: OrderedDict
        Optimal fluxes.
    '''
        
    def __init__(self, obj, opt_obj, opt_fluxes):
        '''
        Parameters
        ----------
        obj: dict
            Reaction ID => coefficient, i.e., the objective function.
        opt_obj: float
            Optimal objective.
        opt_fluxes: OrderedDict
            Optimal fluxes.
        '''
        
        self._obj = obj
        self._opt_obj = opt_obj
        self._opt_fluxes = opt_fluxes
        
    
    @property
    def objective(self):
        
        objStr = ''
        count = 0
        for fluxid, coe in self._obj.items():
            if coe > 0 and count > 0:
                objStr += '+%s*%s' % (coe, fluxid)
            else:
                objStr += '%s*%s' % (coe, fluxid)
            count += 1    
                
        return objStr
    
    
    @property
    def opt_objective(self):
        
        return round(self._opt_obj, 3)
        
        
    @property
    def opt_fluxes(self):
        
        return pDict(self._opt_fluxes)
        
    
    def __repr__(self):
    
        return 'objective: %s\noptimal objective: %s\noptimal fluxes\n%s' % (self.objective, 
                                                                             self.opt_objective, 
                                                                             self.opt_fluxes)
        



class FVAResults():
    '''
    Parameters
    ----------
    flux_ranges: dict
        Reaction ID => [lower bound, upper bound].
    '''
    
    def __init__(self, flux_ranges):
        '''
        Parameters
        ----------
        flux_ranges: dict
            Reaction ID => [lower bound, upper bound].
        '''
        
        self._flux_ranges = flux_ranges
    
        
    @property    
    def flux_ranges(self):
        
        return pDict(self._flux_ranges)
        
    
    def __repr__(self):
        
        return 'flux ranges\n%s' % self.flux_ranges
        
        
        
        
class SimResults():
    '''
    Parameters
    ----------
    simulated_MDVs: dict
        EMU ID => MDV.
    '''
    
    def __init__(self, simulated_MDVs):
        '''
        Parameters
        ----------
        simulated_MDVs: dict
            EMU ID => MDV.
        '''
        
        self._simulated_MDVs = simulated_MDVs
        self._simulated_EMUs = sorted(self._simulated_MDVs.keys())
    
    
    @property
    def simulated_EMUs(self):
        
        return self._simulated_EMUs
        
    
    def simulated_MDV(self, emuid):
        '''
        Parameters
        emuid: str
            EMU ID
        '''
        
        return self._simulated_MDVs[emuid]
        
        
    def __repr__(self):
        
        return 'simulated MDVs\n%s' % pDict(self._simulated_MDVs)
        
        
        
        
class InstSimResults():
    '''
    Parameters
    ----------
    simulated_inst_MDVs: dict
        EMU IDs => {timepoints => MDV}.
    '''
    
    def __init__(self, simulated_inst_MDVs):
        '''
        Parameters
        ----------
        simulated_inst_MDVs: dict
            EMU IDs => {timepoints => MDV}.
        '''
        
        self._simulated_inst_MDVs = simulated_inst_MDVs
        self._simulated_EMUs = sorted(self._simulated_inst_MDVs)
        self._timepoints = sorted(self._simulated_inst_MDVs[self._simulated_EMUs[0]])
        
    
    @property
    def simulated_EMUs(self):
        
        return self._simulated_EMUs
        
        
    @property
    def timepoints(self):
    
        return self._timepoints
        
    
    def simulated_MDV(self, emuid):
        '''
        Parameters
        ----------
        emuid: str
            EMU ID.
        '''
        
        return pDict(self._simulated_inst_MDVs[emuid])        
        
        
    def plot_MDV_kinetics(self, emuid, show_fig = True, output_dir = None):
        '''
        Parameters
        ----------
        emuid: str
            EMU ID.
        show_fig: bool
            Whether to show figure.
        output_dir: str
            Output directory.
        '''    
        
        _MDV_kinetics(emuid, self._simulated_inst_MDVs[emuid], show_fig, output_dir)
        
        
    def __repr__(self):
        
        mdvstrs = []
        for emuid in self.simulated_EMUs:
            mdvstr = '%s\n%s' % (emuid, self.simulated_MDV(emuid))
            mdvstrs.append(mdvstr)
            
        return 'simulated_MDVs\n' + '\n\n'.join(mdvstrs)
        
        
        
        
class FitResults():
    '''
    Parameters
    ----------
    opt_total_fluxes: ser
        Total fluxes at optimal objective.
    opt_net_fluxes: ser
        Net fluxes at optimal objective.
    opt_obj: float
        Optimal value of objective.
    opt_resids: array
        Optimal weighted residuals.
    n_meas: int
        # of measurements.
    n_params: int
        # of parameters.
    sim_MDVs: dict
        EMU ID => simulated MDV.
    exp_MDVs: dict
        EMU ID => [means, sds].
    sim_fluxes: dict
        Flux ID => simulated flux.
    exp_fluxes: dict
        Flux ID => [mean, sd].
    hessian: array
        Hessian matrix at convergence.
    null_space: array
        Null space of stoichiometric matrix.
    transform_matrix: array
        Transform matrix from total fluxes to net fluxes.
    sim_MDVs_der_u: array
        Derivative of simulated MDVs w.r.t. free fluxes.
    sim_fluxes_der_u: array
        Derivative of simualted fluxes w.r.t. free fluxes.
    exp_MDVs_inv_cov: array
        Inversed covariance matrix of measured MDVs.
    exp_fluxes_inv_cov: array
        Inversed covariance matrix of measured fluxes.
    is_success: bool
        Whether the optimization is successful.
    '''
    
    def __init__(self, opt_total_fluxes, opt_net_fluxes, opt_obj, opt_resids, n_meas, n_params, sim_MDVs, 
                 exp_MDVs, sim_fluxes, exp_fluxes, hessian, null_space, transform_matrix, sim_MDVs_der_u, 
                 sim_fluxes_der_u, exp_MDVs_inv_cov, exp_fluxes_inv_cov, is_success):
        '''
        Parameters
        ----------
        opt_total_fluxes: ser
            Total fluxes at optimal objective.
        opt_net_fluxes: ser
            Net fluxes at optimal objective.
        opt_obj: float
            Optimal value of objective.
        opt_resids: array
            Optimal weighted residuals.
        n_meas: int
            # of measurements.
        n_params: int
            # of parameters.
        sim_MDVs: dict
            EMU ID => simulated MDV.
        exp_MDVs: dict
            EMU ID => [means, sds].
        sim_fluxes: dict
            Flux ID => simulated flux.
        exp_fluxes: dict
            Flux ID => [mean, sd].
        hessian: array
            Hessian matrix at convergence.
        null_space: array
            Null space of stoichiometric matrix.
        transform_matrix: array
            Transform matrix from total fluxes to net fluxes.
        sim_MDVs_der_u: array
            Derivative of simulated MDVs w.r.t. free fluxes.
        sim_fluxes_der_u: array
            Derivative of simualted fluxes w.r.t. free fluxes.
        exp_MDVs_inv_cov: array
            Inversed covariance matrix of measured MDVs.
        exp_fluxes_inv_cov: array
            Inversed covariance matrix of measured fluxes.
        is_success: bool
            Whether the optimization is successful.
        '''
        
        self._opt_total_fluxes = opt_total_fluxes
        self._opt_net_fluxes = opt_net_fluxes
        self._opt_obj = opt_obj
        self.opt_resids = opt_resids
        self.simulated_MDVs = sim_MDVs
        self.measured_MDVs = exp_MDVs
        self.simulated_fluxes = sim_fluxes
        self.measured_fluxes = exp_fluxes
        self.hessian = hessian
        self.null_space = null_space
        self.transform_matrix = transform_matrix
        self.n_meas = n_meas
        self.n_params = n_params
        self.dof = self.n_meas - self.n_params
        self.sim_MDVs_der_u = sim_MDVs_der_u
        self.sim_fluxes_der_u = sim_fluxes_der_u
        self.exp_MDVs_inv_cov = exp_MDVs_inv_cov
        self.exp_fluxes_inv_cov = exp_fluxes_inv_cov
        self.is_success = is_success
        
    
    @property
    def n_total_fluxes(self):
        
        return len(self._opt_total_fluxes)
        
        
    @property
    def n_net_fluxes(self):
        
        return len(self._opt_net_fluxes)
        
    
    @property
    def n_free_fluxes(self):
        
        return self.null_space.shape[1]
        
    
    @property
    def opt_objective(self):
        
        return round(self._opt_obj, 3)
        
    
    @property    
    def opt_total_fluxes(self):
        
        return pDict(self._opt_total_fluxes)
    
    
    @property
    def opt_net_fluxes(self):
        
        return pDict(self._opt_net_fluxes)


    @property
    def optimization_successful(self):

        return self.is_success    
    
    
    def chi2_test(self, confidence_level = 0.999):
        '''
        This method performs chi square test of the optimal objective.
        Actually, SSR < LB of chi square interval can be also considered as successful.
        
        Parameters
        ----------
        confidence_level: float
            Confidence level, e.g. 0.95 as 95% confidence level.
        '''
        
        _chi2_test(self._opt_obj, self.dof, confidence_level)
        
        
    def plot_normal_probability(self, show_fig = True, output_dir = None):
        '''
        This method performs normal probability plot for residuals.
        
        Parameters
        ----------
        show_fig: bool
            Whether to show figure.
        output_dir: str
            Output directory. 
        '''
        
        _normal_probability(self.opt_resids, show_fig, output_dir)
        
    
    def plot_simulated_vs_measured_MDVs(self, show_fig = True, output_dir = None):
        '''
        This method plots simulated and measured MDVs.
        
        Parameters
        ----------
        show_fig: bool
            Whether to show figure.
        output_dir: str
            Output directory.
        '''
        
        _simulated_vs_measured_MDVs(self.simulated_MDVs, self.measured_MDVs, show_fig, output_dir)
        
        
    def plot_simulated_vs_measured_fluxes(self, show_fig = True, output_dir = None):
        '''
        This method plots simulated and measured fluxes.
        
        Parameters
        ----------
        show_fig: bool
            Whether to show figure.
        output_dir: str
            Output directory.
        '''
        
        _simulated_vs_measured_fluxes(self.simulated_fluxes, self.measured_fluxes, show_fig, output_dir)
        
    
    def estimate_confidence_intervals(self, which = 'net', confidence_level = 0.95):
        '''
        This method calculates CI of net (total) fluxes using local estimation.
        
        Parameters
        ----------
        which: {"net", "total"}
            * "net" if net fluxes.
            * "total" if total fluxes.
        confidence_level: float
            Confidence level, e.g. 0.95 as 95% confidence level.    
        '''
        
        if which == 'net':
            totalFluxesCov = self.null_space@pinv2(self.hessian)@self.null_space.T
            netFluxesCov = self.transform_matrix@totalFluxesCov@self.transform_matrix.T
            irrRxns = (self._opt_net_fluxes.index&self._opt_total_fluxes.index).tolist()
            netFluxesRange = _confidence_intervals_le(self._opt_net_fluxes, irrRxns, 
                                                      netFluxesCov, self.dof, confidence_level)
            
            return pDict(netFluxesRange)       
        
        elif which == 'total':
            totalFluxesCov = self.null_space@pinv2(self.hessian)@self.null_space.T
            irrRxns = self._opt_total_fluxes.index.tolist()
            totalFluxesRange = _confidence_intervals_le(self._opt_total_fluxes, irrRxns, 
                                                        totalFluxesCov, self.dof, confidence_level)
            
            return pDict(totalFluxesRange)
        
        else:
            raise ValueError('only "net" and "total" are acceptable for which argument')


    def estimate_contribution_matrix(self, which = 'net'):
        '''
        This method calculates contribution matrix of measurement variance to net (total) flux variance.

        Parameters
        ----------
        which: {"net", "total"}
            * "net" if net fluxes.
            * "total" if total fluxes.
        '''

        freeFluxesCov = pinv2(self.hessian)
        expMDVsCov = pinv2(self.exp_MDVs_inv_cov)
        expFluxesCov = pinv2(self.exp_fluxes_inv_cov)

        if which == 'net':
            transMat = self.transform_matrix@self.null_space
            fluxIdx = self._opt_net_fluxes.index
        elif which == 'total':
            transMat = self.null_space
            fluxIdx = self._opt_total_fluxes.index
        else:
            raise ValueError('only "net" and "total" are acceptable for which argument')    

        contribMat_MDV = _contribution_matrix(freeFluxesCov, transMat, self.sim_MDVs_der_u, expMDVsCov)
        contribMat_MDV = pd.DataFrame(contribMat_MDV, index = fluxIdx, 
                                      columns = self._get_name_of_measurements(self.measured_MDVs))

        contribMat_flux = _contribution_matrix(freeFluxesCov, transMat, self.sim_fluxes_der_u, expFluxesCov)
        contribMat_flux = pd.DataFrame(contribMat_flux, index = fluxIdx,
                                       columns = self.measured_fluxes.keys())
        
        contribMat = pd.concat((contribMat_MDV, contribMat_flux), axis = 1)

        return contribMat


    def estimate_sensitivity(self, which = 'net'):
        '''
        This method calculates sensitivity matrix of estimated net (total) flux w.r.t. measurement changes
        
        Parameters
        ----------
        which: {"net", "total"}
            * "net" if net fluxes.
            * "total" if total fluxes.
        '''

        freeFluxesCov = pinv2(self.hessian)

        if which == 'net':
            transMat = self.transform_matrix@self.null_space
            fluxIdx = self._opt_net_fluxes.index
        elif which == 'total':
            transMat = self.null_space
            fluxIdx = self._opt_total_fluxes.index
        else:
            raise ValueError('only "net" and "total" are acceptable for which argument')    

        senMat_MDV = _sensitivity(freeFluxesCov, transMat, self.sim_MDVs_der_u, self.exp_MDVs_inv_cov)
        senMat_MDV = pd.DataFrame(senMat_MDV, index = fluxIdx, 
                                  columns = self._get_name_of_measurements(self.measured_MDVs))
        
        senMat_flux = _sensitivity(freeFluxesCov, transMat, self.sim_fluxes_der_u, self.exp_fluxes_inv_cov)
        senMat_flux = pd.DataFrame(senMat_flux, index = fluxIdx,
                                   columns = self.measured_fluxes.keys())
        
        senMat = pd.concat((senMat_MDV, senMat_flux), axis = 1)

        return senMat


    @staticmethod
    def _get_name_of_measurements(measured_MDVs):
        '''
        Parameters
        ----------
        measured_MDVs: dict
            EMU ID => [means, sds].
        '''

        names = []
        for emuid in measured_MDVs:
            _, atoms = emuid.split('_')
            for atom in '0'+atoms:
                names.append('%s_m%s' %(emuid, atom))

        return names        

    
    def __repr__(self):
    
        return 'optimal objective: %s at\n%s' % (self.opt_objective, self.opt_net_fluxes)
    
        
    
    
class FitMCResults():
    '''
    Parameters
    ----------
    total_fluxes_set: list of ser
        Set of optimal total fluxes.
    net_fluxes_set: list of ser
        Set of optimal net fluxes.    
    '''
    
    def __init__(self, total_fluxes_set, net_fluxes_set):
        '''
        Parameters
        ----------
        total_fluxes_set: list of ser
            Set of optimal total fluxes.
        net_fluxes_set: list of ser
            Set of optimal net fluxes.    
        '''
        
        self.total_fluxes_set = total_fluxes_set
        self.net_fluxes_set = net_fluxes_set
        
        
    def estimate_confidence_intervals(self, which = 'net', confidence_level = 0.95):
        '''
        This method estimates CI from a set of fluxes.
        
        Parameters
        ----------
        which: {"net", "total"}
            * "net" if net fluxes.
            * "total" if total fluxes.
        confidence_level: float
            Confidence level, e.g. 0.95 as 95% confidence level.    
        '''
        
        if which == 'net':
            irrRxns = (self.net_fluxes_set[0].index&self.total_fluxes_set[0].index).tolist()
            netFluxesRange = _confidence_intervals_mc(self.net_fluxes_set, irrRxns, confidence_level)
            return pDict(netFluxesRange)
        
        elif which == 'total':
            irrRxns = self.total_fluxes_set[0].index.tolist()
            totalFluxesRange = _confidence_intervals_mc(self.total_fluxes_set, irrRxns, confidence_level)
            return pDict(totalFluxesRange)
        
        else:
            raise ValueError('only "net" and "total" are acceptable for which argument')
        
        
    def __repr__(self):
        
        return '%s' % self.estimate_confidence_intervals('net', 0.95)
    
        
        
        
class InstFitResults(FitResults):
    '''
    Parameters
    ----------
    opt_total_fluxes: ser
        Total fluxes at optimal objective.
    opt_net_fluxes: ser
        Net fluxes at optimal objective.
    opt_concs: ser
        Concentrations at optimal objective.
    opt_obj: float
        Optimal value of objective.
    opt_resids: array
        Optimal weighted residuals.
    n_meas: int
        # of measurements.
    n_params: int
        # of parameters.
    sim_inst_MDVs: dict
        EMU ID => {t => simulated MDV}.
    exp_inst_MDVs: dict
        EMU ID => {t => [means, sds]}.
    sim_fluxes: dict
        Flux ID => simulated flux.
    exp_fluxes: dict
        Flux ID => [mean, sd].
    hessian: array
        Hessian matrix at convergence.
    null_space: array
        Null space of stoichiometric matrix.
    transform_matrix: array
        Transform matrix from total fluxes to net fluxes.
    sim_inst_MDVs_der_u: array
        Derivative of simulated MDVs w.r.t. free fluxes.
    sim_fluxes_der_u: array
        Derivative of simualted fluxes w.r.t. free fluxes.
    exp_inst_MDVs_inv_cov: array
        Inversed covariance matrix of measured MDVs.
    exp_fluxes_inv_cov: array
        Inversed covariance matrix of measured fluxes.
    is_success: bool
        Whether the optimization is successful.
    '''
    
    def __init__(self, opt_total_fluxes, opt_net_fluxes, opt_concs, opt_obj, opt_resids, n_meas, n_params, 
                 sim_inst_MDVs, exp_inst_MDVs, sim_fluxes, exp_fluxes, hessian, null_space, transform_matrix,
                 sim_inst_MDVs_der_u, sim_fluxes_der_u, exp_inst_MDVs_inv_cov, exp_fluxes_inv_cov, is_success):
        '''
        Parameters
        ----------
        opt_total_fluxes: ser
            Total fluxes at optimal objective.
        opt_net_fluxes: ser
            Net fluxes at optimal objective.
        opt_concs: ser
            Concentrations at optimal objective.
        opt_obj: float
            Optimal value of objective.
        opt_resids: array
            Optimal weighted residuals.
        n_meas: int
            # of measurements.
        n_params: int
            # of parameters.
        sim_inst_MDVs: dict
            EMU ID => {t => simulated MDV}.
        exp_inst_MDVs: dict
            EMU ID => {t => [means, sds]}.
        sim_fluxes: dict
            Flux ID => simulated flux.
        exp_fluxes: dict
            Flux ID => [mean, sd].
        hessian: array
            Hessian matrix at convergence.
        null_space: array
            Null space of stoichiometric matrix.
        transform_matrix: array
            Transform matrix from total fluxes to net fluxes.
        sim_inst_MDVs_der_u: array
            Derivative of simulated MDVs w.r.t. free fluxes.
        sim_fluxes_der_u: array
            Derivative of simualted fluxes w.r.t. free fluxes.
        exp_inst_MDVs_inv_cov: array
            Inversed covariance matrix of measured MDVs.
        exp_fluxes_inv_cov: array
            Inversed covariance matrix of measured fluxes.
        is_success: bool
            Whether the optimization is successful.
        '''
        
        self._opt_total_fluxes = opt_total_fluxes
        self._opt_net_fluxes = opt_net_fluxes
        self._opt_concs = opt_concs
        self._opt_obj = opt_obj
        self.opt_resids = opt_resids
        self.simulated_inst_MDVs = sim_inst_MDVs
        self.measured_inst_MDVs = exp_inst_MDVs
        self.simulated_fluxes = sim_fluxes
        self.measured_fluxes = exp_fluxes
        self.hessian = hessian
        self.null_space = null_space
        self.transform_matrix = transform_matrix
        self.n_meas = n_meas
        self.n_params = n_params
        self.dof = self.n_meas - self.n_params
        self.sim_inst_MDVs_der_u = sim_inst_MDVs_der_u
        self.sim_fluxes_der_u = sim_fluxes_der_u
        self.exp_inst_MDVs_inv_cov = exp_inst_MDVs_inv_cov
        self.exp_fluxes_inv_cov = exp_fluxes_inv_cov
        self.is_success = is_success
        
    
    @property
    def n_concentrations(self):
        
        return len(self._opt_concs)
        
    
    @property    
    def opt_concentrations(self):
        
        return pDict(self._opt_concs)

        
    def plot_simulated_vs_measured_MDVs(self, show_fig = True, output_dir = None):
        '''
        This method plots simulated and measured MDVs.
        
        Parameters
        ----------
        show_fig: bool
            Whether to show figure.
        output_dir: str
            Output directory.
        '''
        
        _simulated_vs_measured_inst_MDVs(self.simulated_inst_MDVs, self.measured_inst_MDVs, 
                                         show_fig, output_dir)    
    
    
    def estimate_confidence_intervals(self, which = 'net', confidence_level = 0.95):
        '''
        This method calculates CI of fluxes and concentrations using local estimation.
        
        Parameters
        ----------
        which: {"net", "total", "conc"}
            * "net" if net fluxes.
            * "total" if total fluxes.
            * "conc" if concentrations.
        confidence_level: float
            Confidence level, e.g. 0.95 as 95% confidence level.    
        '''
        
        if which == 'net':
            fluxHessian = self.hessian[:self.n_free_fluxes,:self.n_free_fluxes]
            totalFluxesCov = self.null_space@pinv2(fluxHessian)@self.null_space.T
            netFluxesCov = self.transform_matrix@totalFluxesCov@self.transform_matrix.T
            irrItems = (self._opt_net_fluxes.index&self._opt_total_fluxes.index).tolist()
            netFluxesRange = _confidence_intervals_le(self._opt_net_fluxes, irrItems, 
                                                      netFluxesCov, self.dof, confidence_level)
            
            return pDict(netFluxesRange)                                                   
        
        elif which == 'total':
            fluxHessian = self.hessian[:self.n_free_fluxes,:self.n_free_fluxes]
            totalFluxesCov = self.null_space@pinv2(fluxHessian)@self.null_space.T
            irrItems = self._opt_total_fluxes.index.tolist()
            totalFluxesRange = _confidence_intervals_le(self._opt_total_fluxes, irrItems, 
                                                        totalFluxesCov, self.dof, confidence_level)
            
            return pDict(totalFluxesRange)
        
        elif which == 'conc':
            concHessian = self.hessian[self.n_free_fluxes:,self.n_free_fluxes:]
            concsCov = pinv2(concHessian)
            irrItems = self._opt_concs.index.tolist()
            concsRange = _confidence_intervals_le(self._opt_concs, irrItems, concsCov, 
                                                  self.dof, confidence_level)
            
            return pDict(concsRange)
            
        else:
            raise ValueError('only "net", "total" and "conc" are acceptable for which argument')
            

    def estimate_contribution_matrix(self, which = 'net'):
        '''
        This method calculates contribution matrix of measurement variance to net (total) flux variance.

        Parameters
        ----------
        which: {"net", "total"}
            * "net" if net fluxes.
            * "total" if total fluxes.
        '''

        fluxHessian = self.hessian[:self.n_free_fluxes,:self.n_free_fluxes]
        freeFluxesCov = pinv2(fluxHessian)
        expMDVsCov = pinv2(self.exp_inst_MDVs_inv_cov)
        expFluxesCov = pinv2(self.exp_fluxes_inv_cov)

        if which == 'net':
            transMat = self.transform_matrix@self.null_space
            fluxIdx = self._opt_net_fluxes.index
        elif which == 'total':
            transMat = self.null_space
            fluxIdx = self._opt_total_fluxes.index
        else:
            raise ValueError('only "net" and "total" are acceptable for which argument')    

        contribMat_MDV = _contribution_matrix(freeFluxesCov, transMat, self.sim_inst_MDVs_der_u, expMDVsCov)
        contribMat_MDV = pd.DataFrame(contribMat_MDV, index = fluxIdx, 
                                      columns = self._get_name_of_measurements(self.measured_inst_MDVs))

        contribMat_flux = _contribution_matrix(freeFluxesCov, transMat, self.sim_fluxes_der_u, expFluxesCov)
        contribMat_flux = pd.DataFrame(contribMat_flux, index = fluxIdx,
                                       columns = self.measured_fluxes.keys())
        
        contribMat = pd.concat((contribMat_MDV, contribMat_flux), axis = 1)

        return contribMat


    def estimate_sensitivity(self, which = 'net'):
        '''
        This method calculates sensitivity matrix of estimated net (total) flux w.r.t. measurement changes.
        
        Parameters
        ----------
        which: {"net", "total"}
            * "net" if net fluxes.
            * "total" if total fluxes.
        '''

        fluxHessian = self.hessian[:self.n_free_fluxes,:self.n_free_fluxes]
        freeFluxesCov = pinv2(fluxHessian)

        if which == 'net':
            transMat = self.transform_matrix@self.null_space
            fluxIdx = self._opt_net_fluxes.index
        elif which == 'total':
            transMat = self.null_space
            fluxIdx = self._opt_total_fluxes.index
        else:
            raise ValueError('only "net" and "total" are acceptable for which argument')    

        senMat_MDV = _sensitivity(freeFluxesCov, transMat, self.sim_inst_MDVs_der_u, self.exp_inst_MDVs_inv_cov)
        senMat_MDV = pd.DataFrame(senMat_MDV, index = fluxIdx, 
                                  columns = self._get_name_of_measurements(self.measured_inst_MDVs))
        
        senMat_flux = _sensitivity(freeFluxesCov, transMat, self.sim_fluxes_der_u, self.exp_fluxes_inv_cov)
        senMat_flux = pd.DataFrame(senMat_flux, index = fluxIdx,
                                   columns = self.measured_fluxes.keys())
        
        senMat = pd.concat((senMat_MDV, senMat_flux), axis = 1)

        return senMat
        

    @staticmethod
    def _get_name_of_measurements(measured_MDVs):
        '''
        Parameters
        ----------
        measured_MDVs: dict
            EMU ID => {t => [means, sds]}.
        '''

        names = []
        for emuid in measured_MDVs:
            for t in measured_MDVs[emuid]:
                if t != 0:
                    _, atoms = emuid.split('_')
                    for atom in '0'+atoms:
                        names.append('%s_m%s_%s' %(emuid, atom, t))

        return names




class InstFitMCResults(FitMCResults):
    '''
    Parameters
    ----------
    total_fluxes_set: list of ser
        Set of optimal total fluxes.
    net_fluxes_set: list of ser
        Set of optimal net fluxes.
    concs_set: list of ser
        Set of optimal concentrations.
    '''
    
    def __init__(self, total_fluxes_set, net_fluxes_set, concs_set):
        '''
        Parameters
        ----------
        total_fluxes_set: list of ser
            Set of optimal total fluxes.
        net_fluxes_set: list of ser
            Set of optimal net fluxes.
        concs_set: list of ser
            Set of optimal concentrations.
        '''
        
        super().__init__(total_fluxes_set, net_fluxes_set)
        
        self.concs_set = concs_set


    def estimate_confidence_intervals(self, which = 'net', confidence_level = 0.95):
        '''
        This method estimates CI from a set of fluxes.
        
        Parameters
        ----------
        which: {"net", "total", "conc"}
            * "net" if net fluxes.
            * "total" if total fluxes.
            * "conc" if concentrations.
        confidence_level: float
            Confidence level, e.g. 0.95 as 95% confidence level.    
        '''

        if which == 'net':
            irrItems = (self.net_fluxes_set[0].index&self.total_fluxes_set[0].index).tolist()
            netFluxesRange = _confidence_intervals_mc(self.net_fluxes_set, irrItems, confidence_level)
            return pDict(netFluxesRange)
        
        elif which == 'total':
            irrItems = self.total_fluxes_set[0].index.tolist()
            totalFluxesRange = _confidence_intervals_mc(self.total_fluxes_set, irrItems, confidence_level)
            return pDict(totalFluxesRange)

        elif which == 'conc':
            irrItems = self.concs_set[0].index.tolist()
            concsRange = _confidence_intervals_mc(self.concs_set, irrItems, confidence_level)
            return pDict(concsRange)
        
        else:
            raise ValueError('only "net", "total" and "conc" are acceptable for which argument')
