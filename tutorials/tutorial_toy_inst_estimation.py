'''Example of flux estimation at isotopically nonstationary (INST) state with a toy model.
'''


from os import makedirs
import pandas as pd
from freeflux import Model


MODEL_FILE = '../models/toy/reactions.tsv' 
MEASURED_MDVS = '../models/toy/measured_inst_MDVs.tsv'
MEASURED_FLUXES = '../models/toy/measured_fluxes.tsv'
OUT_DIR = '../results/toy/inst_estimation'


# estimate fluxes and concentrations at INST state
def toy_model_inst_fitting():

    model = Model('demo')
    model.read_from_file(MODEL_FILE)

    with model.fitter('inst') as ifit:
        # specify the lableing strategy, 
        # use this method for every labeled substrate
        ifit.set_labeling_strategy(
            'AcCoA', 
            labeling_pattern = ['01', '11'], 
            percentage = [0.25, 0.25], 
            purity = [1, 1]
        )
        
        # read measurements
        ifit.set_measured_MDVs_from_file(MEASURED_MDVS)
        ifit.set_measured_fluxes_from_file(MEASURED_FLUXES)
        
        # set bounds for fluxes and concentrations
        ifit.set_flux_bounds('all', bounds = [-100, 100])
        ifit.set_concentration_bounds('all', bounds = [0.001, 10])
        
        # solve fluxes and concentrations
        ifit.prepare(n_jobs = 3)
        res = ifit.solve(solver = 'ralg')
    
    # print(res.optimization_successful)

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
def toy_model_inst_fitting_CIs():

    model = Model('demo')
    model.read_from_file(MODEL_FILE)

    with model.fitter('inst') as ifit:
        # specify the lableing strategy, 
        # use this method for every labeled substrate
        ifit.set_labeling_strategy(
            'AcCoA', 
            labeling_pattern = ['01', '11'], 
            percentage = [0.25, 0.25], 
            purity = [1, 1]
        )
        
        # read measurements
        ifit.set_measured_MDVs_from_file(MEASURED_MDVS)
        ifit.set_measured_fluxes_from_file(MEASURED_FLUXES)
        
        # set bounds for fluxes and concentrations
        ifit.set_flux_bounds('all', bounds = [-100, 100])
        ifit.set_concentration_bounds('all', bounds = [0.001, 10])

        # estimate the confidence intervals, 
        # highly recommended to run with parallel jobs
        ifit.prepare(n_jobs = 3)
        res = ifit.solve_with_confidence_intervals(
            solver = 'ralg', 
            n_runs = 30, 
            n_jobs = 3
        )

    # save the CIs
    net_cis = res.estimate_confidence_intervals(
        which = 'net', 
        confidence_level = 0.95
    )
    pd.DataFrame(net_cis, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/netflux_CIs.xlsx'
    )

    total_cis = res.estimate_confidence_intervals(
        which = 'total', 
        confidence_level = 0.95
    )
    pd.DataFrame(total_cis, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/totalflux_CIs.xlsx'
    )

    conc_cis = res.estimate_confidence_intervals(
        which = 'conc', 
        confidence_level = 0.95
    )
    pd.DataFrame(conc_cis, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/concentration_CIs.xlsx'
    )




if __name__ == '__main__':

    makedirs(OUT_DIR, exist_ok = True)
    toy_model_inst_fitting()
    toy_model_inst_fitting_CIs()
