'''Example of Flux estimation at steady state with a E. coli model using experimental data.
'''


from os import makedirs
import pandas as pd
from freeflux import Model


MODEL_FILE = '../models/ecoli/experimental_data/reactions.xlsx' 
MEASURED_MDVS = '../models/ecoli/experimental_data/measured_MDVs.xlsx'
MEASURED_FLUXES = '../models/ecoli/experimental_data/measured_fluxes.xlsx'
OUT_DIR = '../results/ecoli/experimental_data/steady_state_estimation'

DILUTION_FROM = [
    'CO2u', 
    'Alau', 
    'Glyu', 
    'Valu', 
    'Leuu', 
    'Ileu', 
    'Seru', 
    'Pheu', 
    'Aspu', 
    'Gluu', 
    'Tyru'
]


# estimate fluxes at steady state
def steady_state_fitting():

    ecoli = Model('ecoli')
    ecoli.read_from_file(MODEL_FILE)
    
    with ecoli.fitter('ss') as fit:
        # specify the lableing strategy, 
        # use this method for every labeled substrate
        fit.set_labeling_strategy(
            'Glc.ex', 
            labeling_pattern = ['100000', '111111'], 
            percentage = [0.77, 0.205], 
            purity = [0.99, 0.985]
        )
        
        # read measurements
        fit.set_measured_MDVs_from_file(MEASURED_MDVS)
        fit.set_measured_fluxes_from_file(MEASURED_FLUXES)
        
        # set upper and lower bounds for fluxes
        fit.set_flux_bounds('all', bounds = [-100, 100]) 
        
        # solve the fluxes
        fit.prepare(
            dilution_from = DILUTION_FROM, 
            n_jobs = 30
        )
        while True:
            res = fit.solve(solver = 'ralg', max_iters = 1000)
            if res.optimization_successful:
                break
    
    # save the results
    pd.Series(res.opt_net_fluxes).to_excel(
        OUT_DIR+'/estimated_net_fluxes.xlsx'
    )
    pd.Series(res.opt_total_fluxes).to_excel(
        OUT_DIR+'/estimated_total_fluxes.xlsx'
    )

    net_cis = res.estimate_confidence_intervals(
        which = 'net', 
        confidence_level = 0.95
    )
    total_cis = res.estimate_confidence_intervals(
        which = 'total', 
        confidence_level = 0.95
    )
    pd.DataFrame(net_cis, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/netflux_le_CIs.xlsx'
    )
    pd.DataFrame(total_cis, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/totalflux_le_CIs.xlsx'
    )

    # normal probability plot of residuals
    res.plot_normal_probability(show_fig = False, output_dir = OUT_DIR)
    
    # compare simulations and measurements
    res.plot_simulated_vs_measured_MDVs(show_fig = False, output_dir = OUT_DIR)
    res.plot_simulated_vs_measured_fluxes(show_fig = False, output_dir = OUT_DIR)
    
    # export the contribution matrix
    res.estimate_contribution_matrix(which = 'net').to_excel(
        OUT_DIR+'/netflux_contribMat.xlsx'
    )
    res.estimate_contribution_matrix(which = 'total').to_excel(
        OUT_DIR+'/totalflux_contribMat.xlsx'
    )
    
    # export the sensitivity matrix
    res.estimate_sensitivity(which = 'net').to_excel(
        OUT_DIR+'/netflux_senMat.xlsx'
    )
    res.estimate_sensitivity(which = 'total').to_excel(
        OUT_DIR+'/totalflux_senMat.xlsx'
    )


# estimate with confidence intervals
def steady_state_fitting_CIs():

    ecoli = Model('ecoli')
    ecoli.read_from_file(MODEL_FILE)
    
    with ecoli.fitter('ss') as fit:
        # specify the lableing strategy, 
        # use this method for every labeled substrate
        fit.set_labeling_strategy(
            'Glc.ex', 
            labeling_pattern = ['100000', '111111'], 
            percentage = [0.77, 0.205], 
            purity = [0.99, 0.985]
        )
        
        # read measurements
        fit.set_measured_MDVs_from_file(MEASURED_MDVS)
        fit.set_measured_fluxes_from_file(MEASURED_FLUXES)
        
        # set upper and lower bounds for fluxes
        fit.set_flux_bounds('all', bounds = [-100, 100]) 
        
        # estimate confidence intervals, 
        # highly recommended to run with parallel jobs
        fit.prepare(
            dilution_from = DILUTION_FROM, 
            n_jobs = 30
        )
        res = fit.solve_with_confidence_intervals(
            solver = 'ralg', 
            max_iters = 1000,
            n_runs = 500, 
            n_jobs = 30
        )

    # save the CIs
    net_cis = res.estimate_confidence_intervals(
        which = 'net', 
        confidence_level = 0.95
    )
    pd.DataFrame(net_cis, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/netflux_mc_CIs.xlsx'
    )

    total_cis = res.estimate_confidence_intervals(
        which = 'total', 
        confidence_level = 0.95
    )
    pd.DataFrame(total_cis, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/totalflux_mc_CIs.xlsx'
    )




if __name__ == '__main__':

    makedirs(OUT_DIR, exist_ok = True)
    steady_state_fitting()
    steady_state_fitting_CIs()
