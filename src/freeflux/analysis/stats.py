#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '05/08/2022'




from os import makedirs
import numpy as np
import pandas as pd
from scipy.stats import t, chi2, probplot, zscore
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt




def _chi2_test(obj_value, dof, confidence_level):
    '''
    Parameters
    ----------
    obj_value: float
        Objective value.
    dof: int
        Degree of freedom.
    confidence_level: float
        Confidence level, e.g., 0.95 as 95% confidence level.
    '''
    
    chi2Lb = chi2.ppf((1-confidence_level)/2, dof)
    chi2Ub = chi2.ppf((1+confidence_level)/2, dof)
    
    flag = '' if chi2Lb <= obj_value <= chi2Ub else 'not '
    
    print('SSR %.2f %sacceptable, which is %sin [%.2f, %.2f] of %s%% confidence level by chi2 test' % (obj_value, 
                                                                                                       flag,
                                                                                                       flag,
                                                                                                       chi2Lb, 
                                                                                                       chi2Ub, 
                                                                                                       confidence_level*100))


def _normal_probability(resids, show_fig, output_dir):
    '''
    Parameters
    ----------
    resids: array
        Residuals.
    show_fig: bool
        Whether to show figure.
    output_dir: str
        Output directory.
    '''
    
    probplot(resids, plot = plt)
    ax = plt.gca()
    for item in [ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels():
        item.set_fontsize(15)
    ax.set_ylabel('Weighted residuals')
    ax.set_title('')
    
    if output_dir:
        makedirs(output_dir, exist_ok = True)
        plt.savefig('%s/normal_probability_plot.jpg' % output_dir, dpi = 300, bbox_inches = 'tight')
    
    if show_fig:
        plt.show()
    
    plt.close()


def _simulated_vs_measured_measurements(sim_meas, exp_meas, exp_meas_err, xlabel, ylabel, 
                                        xticklabels, filename, show_fig, output_dir):
    '''
    Parameters
    ----------
    sim_meas: array or list
        Simulated measurements.
    exp_meas: array or list
        Measured measurements.
    exp_meas_err: array or list
        Errors of measured measurements.
    xlabel: str
        Xlabel.
    ylabel: str
        Ylabel.
    xticklabels: array or list
        Xticklabels.
    filename: str
        File name.    
    show_fig: bool
        Whether to show figure.
    output_dir: str
        Output directory.
    '''
    
    sim_meas = np.array(sim_meas)
    exp_meas = np.array(exp_meas)
    exp_meas_err = np.array(exp_meas_err)
    
    barWidth = 0.4
    xPos = np.arange(1, sim_meas.size+1)
    
    plt.bar(x = xPos - 0.5 * barWidth, height = sim_meas, width = barWidth, label = 'simulated')
    plt.bar(x = xPos + 0.5 * barWidth, height = exp_meas, yerr = exp_meas_err, width = barWidth, 
            label = 'measured', capsize = 3)
    plt.xticks(ticks = xPos, labels = xticklabels, fontsize = 15)
    plt.tick_params(axis = 'y', labelsize = 15)
    plt.xlabel(xlabel, fontsize = 15)
    plt.ylabel(ylabel, fontsize = 15)
    
    plt.legend(loc = 'center', bbox_to_anchor = (1.2, 0.5), fontsize = 15, frameon = False)
    
    if output_dir:
        makedirs(output_dir, exist_ok = True)
        #filename = xlabel if xlabel else ylabel
        plt.savefig('%s/%s.jpg' % (output_dir, filename), dpi = 300, bbox_inches = 'tight')
    
    if show_fig:
        plt.show()
    
    plt.close()
    
    
def _simulated_vs_measured_MDVs(simulated_MDVs, measured_MDVs, show_fig, output_dir):
    '''
    Parameters
    ----------
    simulated_MDVs: dict
        EMU ID => simulated MDV.
    measured_MDVs: dict
        EMU ID => [means, sds].
    show_fig: bool
        Whether to show figure.
    output_dir: str
        Output directory.    
    '''
    
    for emuid in measured_MDVs:
        simMDV = simulated_MDVs[emuid]
        expMDV, expMDVerr = measured_MDVs[emuid]
        xticklabels = ['M%s' % str(i) for i in range(simMDV.size)]
        _simulated_vs_measured_measurements(simMDV, expMDV, expMDVerr, emuid, 'Isotopomer fraction', 
                                            xticklabels, 'sim_vs_exp_MDVs-'+emuid, show_fig, output_dir)
        
    
def _simulated_vs_measured_fluxes(simulated_fluxes, measured_fluxes, show_fig, output_dir):
    '''
    Parameters
    ----------
    simulated_fluxes: dict
        Flux ID => simulated flux.
    measured_fluxes: dict
        Flux ID => [mean, sd].
    show_fig: bool
        Whether to show figure.
    output_dir: str
        Output directory.
    '''
    
    simFluxes = []
    expFluxes = []
    expFluxesErr = []
    xticklabels = []
    for fluxid, simFlux in simulated_fluxes.items():
        expFlux, expFluxErr = measured_fluxes[fluxid]
        simFluxes.append(simFlux)
        expFluxes.append(expFlux)
        expFluxesErr.append(expFluxErr)
        xticklabels.append(fluxid)
    
    _simulated_vs_measured_measurements(simFluxes, expFluxes, expFluxesErr, '', 'Flux', 
                                        xticklabels, 'sim_vs_exp_fluxes', show_fig, output_dir)
    


def _confidence_intervals_le(res, irr, cov, dof, confidence_level):
    '''
    Parameters
    ----------
    res: ser
        Optimal results, e.g., net fluxes, total fluxes or concentrations.
    irr: list
        Irreversible items. total fluxes and concentrations are all considered
        irreversible.
    cov: array
        Corvariance matrix of free fluxes obtained from hessian at convergence.
    dof: int
        Degree of freedom.
    confidence_level: float
        Confidence level, e.g., 0.95 as 95% confidence level.
    '''
    
    errors = t.ppf((1+confidence_level)/2, dof)/(dof+1)**0.5*np.diag(cov)**0.5
    lbs = res - errors
    ubs = res + errors
    
    ranges = {}
    for varid in res.index:
        if varid in irr:
            ranges[varid] = [max(0, lbs[varid]), ubs[varid]]
        else:
            ranges[varid] = [lbs[varid], ubs[varid]]
   
    return ranges
    

def _confidence_intervals_mc(res_set, irr, confidence_level):
    '''
    Parameters
    ----------
    res_set: list of ser
        Set of optimal net fluxes. Total fluxes or concentrations.
    irr: list
        Irreversible items. Total fluxes and concentrations are all considered
        irreversible.
    confidence_level: float
        Confidence level, e.g., 0.95 as 95% confidence level.
    '''
    
    resSet = pd.DataFrame(res_set)
    resSet = resSet[(np.abs(zscore(resSet)) < 3).all(axis = 1)]
    
    quantile = resSet.quantile([(1-confidence_level)/2, (1+confidence_level)/2])
    mask = resSet.apply(lambda col: col.between(*quantile[col.name]), axis = 0)
    selected = resSet[mask.all(axis = 1)]

    if selected.empty:
        raise ValueError('can not estimate CIs, need more Monte Carlo runs')
    else:
        ranges = {}
        for varid, col in selected.iteritems():
            if varid in irr:
                ranges[varid] = [max(0, col.min()), col.max()]
            else:
                ranges[varid] = [col.min(), col.max()]
      
        return ranges


def _MDV_kinetics(emuid, simulated_inst_MDVs, show_fig, output_dir):
    '''
    Parameters
    ----------
    emuid: str
        EMU ID.
    simulated_inst_MDVs: dict
        Timepoint => MDV.
    show_fig: bool
        Whether to show figure.
    output_dir: str
        Output directory.
    '''
    
    tpoints = sorted(simulated_inst_MDVs)
    ts = np.linspace(tpoints[0], tpoints[-1], 100)
    mdvs = np.array(list(simulated_inst_MDVs.values()))
    interpMDVs = interp1d(tpoints, mdvs, axis = 0, kind = 1)(ts)
    
    for i in range(interpMDVs.shape[1]):
        plt.plot(ts, interpMDVs[:, i], label = 'M%s' % i, linewidth = 2)
        
    plt.xlabel('Time (s)', fontsize = 15)
    plt.ylabel('EMU %s' % emuid, fontsize = 15)
    plt.tick_params(labelsize = 15)
    plt.legend(loc = 'center', bbox_to_anchor = (1.15, 0.5), fontsize = 15, frameon = False)
    
    if output_dir:
        makedirs(output_dir, exist_ok = True)
        plt.savefig('%s/%s_kinetics.jpg' % (output_dir, emuid), dpi = 300, bbox_inches = 'tight')
    
    if show_fig:
        plt.show()
    
    plt.close()
    
    
def _simulated_vs_measured_inst_MDVs(simulated_inst_MDVs, measured_inst_MDVs, show_fig, output_dir):
    '''
    Parameters
    ----------
    simulated_inst_MDVs: dict
        EMU ID => {t => simulated MDV}.
    measured_inst_MDVs: dict
        EMU ID => {t => [means, sds]}.
    show_fig: bool
        Whether to show figure.
    output_dir: str
        Output directory.
    '''
    
    for emuid in measured_inst_MDVs:
        simInstMDVs = simulated_inst_MDVs[emuid]
        expInstMDVs = measured_inst_MDVs[emuid]
    
        tpoints_sim = sorted(simInstMDVs)
        ts = np.linspace(tpoints_sim[0], tpoints_sim[-1], 100)
        mdvs = np.array(list(simInstMDVs.values()))
        interpMDVs = interp1d(tpoints_sim, mdvs, axis = 0, kind = 1)(ts)
        
        tpoint_exp = sorted(expInstMDVs)
        expMDVs = np.array([expInstMDVs[t][0] for t in tpoint_exp])
        expMDVerrs = np.array([expInstMDVs[t][1] for t in tpoint_exp])
        
        colors = []
        for i, simMDV in enumerate(interpMDVs.T):
            p = plt.plot(ts, simMDV, label = 'M%s' % i, linewidth = 2)
            colors.append(p[0].get_color())
        legend = plt.legend(loc = 'center', bbox_to_anchor = (1.15, 0.5), title = 'simulated', 
                            title_fontsize = 15, fontsize = 15, frameon = False)
        plt.gca().add_artist(legend)    
        
        handles = []
        labels = []
        for i, (expMDV, expMDVerr) in enumerate(zip(expMDVs.T, expMDVerrs.T)):
            e = plt.errorbar(tpoint_exp, expMDV, expMDVerr, color = colors[i], linestyle = '',
                             marker = '.', markersize = 7, elinewidth = 2, capsize = 3)
            handles.append(e[0])
            labels.append('M%s' % i)
        plt.legend(handles = handles, labels = labels, loc = 'center', bbox_to_anchor = (1.4, 0.5), 
                   title = 'measured', title_fontsize = 15, fontsize = 15, frameon = False)
            
        plt.xlabel('Time (s)', fontsize = 15)
        plt.ylabel('EMU %s' % emuid, fontsize = 15)
        plt.tick_params(labelsize = 15)
        
        if output_dir:
            makedirs(output_dir, exist_ok = True)
            plt.savefig('%s/%s_inst.jpg' % (output_dir, emuid), dpi = 300, bbox_inches = 'tight')
        
        if show_fig:
            plt.show()
        
        plt.close()
            

def _contribution_matrix(cov, trans_mat, simulated_der, measured_cov):
    '''
    Parameters
    ----------
    cov: array
        Corvariance matrix of free fluxes obtained from hessian at convergence.
    trans_mat: array
        Transformation matrix from free flux to total flux or net flux, i.e.,
        N to total flux, T@N to net flux.
    simulated_der: array
        Derivative of simulated measurements w.r.t. free fluxes.
    measured_cov: array
        Corvariance matrix of measurements.
    '''

    fluxesCov = trans_mat@cov@trans_mat.T
    fluxes_measCov = trans_mat@cov@simulated_der.T

    contribMat = np.empty_like(fluxes_measCov)
    for i in range(contribMat.shape[0]):
        for j in range(contribMat.shape[1]):
            contribMat[i,j] = fluxes_measCov[i,j]**2/(fluxesCov[i,i]*measured_cov[j,j])

    return contribMat


def _sensitivity(cov, trans_mat, simulated_der, measured_inv_cov):
    '''
    Parameters
    ----------
    cov: array
        Corvariance matrix of free fluxes obtained from hessian at convergence.
    trans_mat: array
        Transformation matrix from free flux to total flux or net flux, i.e.,
        N to total flux, T@N to net flux.
    simulated_der: array
        Derivative of simulated measurements w.r.t. free fluxes.
    measured_inv_cov: array
        Inversed corvariance matrix of measurements.
    '''

    return trans_mat@cov@simulated_der.T@measured_inv_cov
