'''Example of flux balabce analysis with a E. coli model.
'''


from os import makedirs
import pandas as pd
from freeflux import Model


MODEL_FILE = '../models/ecoli/synthetic_data/reactions.xlsx'
OUT_DIR = '../results/ecoli/fba'


# Flux estimation by Flux Balance Analysis
def ecoli_fba():
    
    model = Model('ecoli')
    model.read_from_file(MODEL_FILE)
    
    with model.optimizer() as opt:
        # set bounds for fluxes
        opt.set_flux_bounds('all', bounds = [-100, 100]) 
        opt.set_flux_bounds('glk', bounds = [10, 10])
        
        # FBA
        opt.prepare()
        res = opt.optimize(objective = {'biom': 1})
        
    # save the results
    pd.Series(res.opt_fluxes).to_excel(OUT_DIR+'/estimated_fluxes.xlsx')


# Flux range estimation by Flux Variability Analysis
def ecoli_fva():

    model = Model('ecoli')
    model.read_from_file(MODEL_FILE)
    
    with model.optimizer() as opt:
        # set bounds for fluxes
        opt.set_flux_bounds('all', bounds = [-100, 100]) 
        opt.set_flux_bounds('glk', bounds = [10, 10])
        
        # FVA
        opt.prepare()
        res = opt.estimate_fluxes_range(
            objective = {'biom': 1}, 
            gamma = 0
        )
    
    # save the results
    pd.DataFrame(res.flux_ranges, index = ['LB', 'UB']).T.to_excel(
        OUT_DIR+'/estimated_flux_ranges.xlsx'
    )




if __name__ == '__main__':
    
    makedirs(OUT_DIR, exist_ok = True)    
    ecoli_fba()
    ecoli_fva()
