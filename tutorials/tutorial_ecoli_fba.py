#!/usr/bin/env python
# -*- coding: UTF-8 -*-


r'''
python C:\Users\cwu\Desktop\Software\FreeFlux\tutorials\tutorial_ecoli_fba.py   # TODO delete
'''


MODEL_FILE = r'C:\Users\cwu\Desktop\Software\FreeFlux/models/ecoli/reactions.xlsx'
OUT_DIR = r'C:\Users\cwu\Desktop\Software\FreeFlux/results/ecoli/fba'
# TODO
# MODEL_FILE = '../models/ecoli/reactions.xlsx'
# OUT_DIR = '..results/ecoli/fba'


from os import makedirs
import pandas as pd
from freeflux import Model


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
        
    # print(res)
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
        res = opt.estimate_fluxes_range(objective = {'biom': 1}, gamma = 0)
    
    # print(res)
    pd.DataFrame(res.flux_ranges, index = ['LB', 'UB']).T.to_excel(OUT_DIR+'/estimated_flux_ranges.xlsx')




if __name__ == '__main__':
    
    makedirs(OUT_DIR, exist_ok = True)    
    ecoli_fba()
    ecoli_fva()
