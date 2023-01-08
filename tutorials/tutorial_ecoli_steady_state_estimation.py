#!/usr/bin/env python
# -*- coding: UTF-8 -*-



MODEL_FILE = '../models/ecoli/reactions.xlsx' 
MEASURED_MDVS = '../models/ecoli/measured_MDVs.xlsx'
MEASURED_FLUXES = '../models/ecoli/measured_fluxes.xlsx'
OUT_DIR = '..results/ecoli/steady_state_estimation'


from os import makedirs
import pandas as pd
from freeflux import Model


# estimate fluxes at steadt state
def ecoli_steady_state_fitting():

    ecoli = Model('ecoli')
    ecoli.read_from_file(MODEL_FILE)
    
    with ecoli.fitter('ss') as fit:
        # specify the lableing strategy
        fit.set_labeling_strategy('Glc.ex', 
                                  labeling_pattern = ['100000', '111111'], 
                                  percentage = [0.754, 0.246], 
                                  purity = [0.997, 0.994])   # call this method for each labeled substrate
        
        # read measurements
        fit.set_measured_MDVs_from_file(MEASURED_MDVS)
        fit.set_measured_fluxes_from_file(MEASURED_FLUXES)
        
        # set upper and lower bounds for fluxes
        fit.set_flux_bounds('all', bounds = [-100, 100]) 
        
        # solve the fluxes
        fit.prepare(n_jobs = 3)
        while True:
            res = fit.solve(solver = 'ralg')
            if res.optimization_successful:
                break
    
    # save the results
    pd.Series(res.opt_net_fluxes).to_excel(OUT_DIR+'/estimated_net_fluxes.xlsx')
    pd.Series(res.opt_total_fluxes).to_excel(OUT_DIR+'/estimated_total_fluxes.xlsx')

    net_cis = res.estimate_confidence_intervals(which = 'net', confidence_level = 0.95)
    pd.DataFrame(net_cis, index = ['LB', 'UB']).T.to_excel(OUT_DIR+'/netflux_le_CIs.xlsx')
    
    total_cis = res.estimate_confidence_intervals(which = 'total', confidence_level = 0.95)
    pd.DataFrame(total_cis, index = ['LB', 'UB']).T.to_excel(OUT_DIR+'/totalflux_le_CIs.xlsx')

    # normal probability plot of residuals
    res.plot_normal_probability(show_fig = False, output_dir = OUT_DIR)
    
    # compare simulations and measurements
    res.plot_simulated_vs_measured_MDVs(show_fig = False, output_dir = OUT_DIR)
    res.plot_simulated_vs_measured_fluxes(show_fig = False, output_dir = OUT_DIR)
    
    # export the contribution matrix of measurements to the uncertainties of estimated fluxes
    res.estimate_contribution_matrix(which = 'net').to_excel(OUT_DIR+'/netflux_contribMat.xlsx')
    res.estimate_contribution_matrix(which = 'total').to_excel(OUT_DIR+'/totalflux_contribMat.xlsx')
    
    # export the sensitivity matrix of estimated fluxes w.r.t. measurement
    res.estimate_sensitivity(which = 'net').to_excel(OUT_DIR+'/netflux_senMat.xlsx')
    res.estimate_sensitivity(which = 'total').to_excel(OUT_DIR+'/totalflux_senMat.xlsx')


# estimate the confidence intervals of fluxes
def ecoli_steady_state_fitting_CIs():

    ecoli = Model('ecoli')
    ecoli.read_from_file(MODEL_FILE)
    
    with ecoli.fitter('ss') as fit:
        # specify the lableing strategy
        fit.set_labeling_strategy('Glc.ex', 
                                  labeling_pattern = ['100000', '111111'], 
                                  percentage = [0.754, 0.246], 
                                  purity = [0.997, 0.994])   # call this method for each labeled substrate
        
        # read measurements
        fit.set_measured_MDVs_from_file(MEASURED_MDVS)
        fit.set_measured_fluxes_from_file(MEASURED_FLUXES)
        
        # set upper and lower bounds for fluxes
        fit.set_flux_bounds('all', bounds = [-100, 100]) 
        
        # estimate confidence intervals, highly recommended to run with parallel jobs
        fit.prepare(n_jobs = 3)
        res = fit.solve_with_confidence_intervals(solver = 'ralg', n_runs = 100, n_jobs = 30)
    
    # save the CIs
    net_cis = res.estimate_confidence_intervals(which = 'net', confidence_level = 0.95)
    pd.DataFrame(net_cis, index = ['LB', 'UB']).T.to_excel(OUT_DIR+'/netflux_mc_CIs.xlsx')

    total_cis = res.estimate_confidence_intervals(which = 'total', confidence_level = 0.95)
    pd.DataFrame(total_cis, index = ['LB', 'UB']).T.to_excel(OUT_DIR+'/totalflux_mc_CIs.xlsx')




if __name__ == '__main__':

    makedirs(OUT_DIR, exist_ok = True)
    ecoli_steady_state_fitting()
    ecoli_steady_state_fitting_CIs()
