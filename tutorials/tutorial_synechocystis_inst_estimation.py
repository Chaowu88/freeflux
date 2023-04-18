'''Example of flux estimation at isotopically nonstationary (INST) state with a E. coli model.
'''


from os import makedirs
import pandas as pd
from freeflux import Model


MODEL_FILE = '../models/synechocystis/reactions.xlsx' 
MEASURED_MDVS = '../models/synechocystis/measured_inst_MDVs.xlsx'
MEASURED_FLUXES = '../models/synechocystis/measured_fluxes.xlsx'
OUT_DIR = '..results/synechocystis/inst_estimation'


# estimate fluxes and concentrations at INST state
def syn_inst_fitting():

    model = Model('cyano')
    model.read_from_file(MODEL_FILE)
    
    with model.fitter('inst') as ifit:
        # specify the lableing strategy, use this method for every labeled substrate
        ifit.set_labeling_strategy(
            'CO2.ex', 
            labeling_pattern = ['1'], 
            percentage = [0.5], 
            purity = [0.997]
        )
        
        # read measurements
        ifit.set_measured_MDVs_from_file(MEASURED_MDVS)
        ifit.set_measured_fluxes_from_file(MEASURED_FLUXES)
        
        # set bounds for fluxes and concentrations
        ifit.set_flux_bounds('all', bounds = [-100, 100])
        ifit.set_concentration_bounds('all', bounds = [0.1, 100])
        
        # solve fluxes and concentrations
        ifit.prepare(n_jobs = 3)
        while True:
            res = ifit.solve(solver = 'ralg')
            if res.optimization_successful:
                break

    # save the results
    pd.Series(res.opt_net_fluxes).to_excel(
        OUT_DIR+'/estimated_net_fluxes.xlsx'
    )
    pd.Series(res.opt_total_fluxes).to_excel(
        OUT_DIR+'/estimated_total_fluxes.xlsx'
    )
    pd.Series(res.opt_concentrations).to_excel(
        OUT_DIR+'/estimated_concentrations.xlsx'
    )

    net_cis = res.estimate_confidence_intervals(
        which = 'net', 
        confidence_level = 0.95
    )
    pd.DataFrame(net_cis, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/netflux_le_CIs.xlsx'
    )

    total_cis = res.estimate_confidence_intervals(
        which = 'total', 
        confidence_level = 0.95
    )
    pd.DataFrame(total_cis, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/totalflux_le_CIs.xlsx'
    )

    conc_cis = res.estimate_confidence_intervals(
        which = 'conc', 
        confidence_level = 0.95
    )
    pd.DataFrame(conc_cis, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/conc_le_CIs.xlsx'
    )

    # normal probability plot of residuals
    res.plot_normal_probability(show_fig = False, output_dir = OUT_DIR)
    
    # compare simulations and measurements
    res.plot_simulated_vs_measured_MDVs(show_fig = False, output_dir = OUT_DIR)
    res.plot_simulated_vs_measured_fluxes(show_fig = False, output_dir = OUT_DIR)
    
    # export the contribution matrix of measurements to the uncertainties of estimated fluxes
    res.estimate_contribution_matrix(which = 'net').to_excel(
        OUT_DIR+'/netflux_contribMat.xlsx'
    )
    res.estimate_contribution_matrix(which = 'total').to_excel(
        OUT_DIR+'/totalflux_contribMat.xlsx'
    )
    
    # export the sensitivity matrix of estimated fluxes w.r.t. measurement
    res.estimate_sensitivity(which = 'net').to_excel(
        OUT_DIR+'/netflux_senMat.xlsx'
    )
    res.estimate_sensitivity(which = 'total').to_excel(
        OUT_DIR+'/totalflux_senMat.xlsx'
    )


def syn_inst_fitting_CIs():

    model = Model('cyano')
    model.read_from_file(MODEL_FILE)
    
    with model.fitter('inst') as ifit:
        # specify the lableing strategy, use this method for every labeled substrate
        ifit.set_labeling_strategy(
            'CO2.ex', 
            labeling_pattern = ['1'], 
            percentage = [0.5], 
            purity = [0.997]
        )
        
        # read measurements
        ifit.set_measured_MDVs_from_file(MEASURED_MDVS)
        ifit.set_measured_fluxes_from_file(MEASURED_FLUXES)
        
        # set bounds for fluxes and concentrations
        ifit.set_flux_bounds('all', bounds = [-100, 100])
        ifit.set_concentration_bounds('all', bounds = [0.1, 100])

        # estimate the confidence intervals, highly recommended to run with parallel jobs
        ifit.prepare(n_jobs = 3)
        res = ifit.solve_with_confidence_intervals(
            solver = 'ralg', 
            max_iters = 800, 
            n_runs = 100, 
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

    conc_cis = res.estimate_confidence_intervals(
        which = 'conc', 
        confidence_level = 0.95
    )
    pd.DataFrame(conc_cis, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/concentration_mc_CIs.xlsx'
    )




if __name__ == '__main__':

    makedirs(OUT_DIR, exist_ok = True)
    syn_inst_fitting()
    syn_inst_fitting_CIs()
