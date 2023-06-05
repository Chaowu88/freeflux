'''Example of MDV simulation at steady state with a E. coli model.
'''


from freeflux import Model


MODEL_FILE = '../models/ecoli/synthetic_data/reactions.xlsx' 
FLUXES = '../models/ecoli/synthetic_data/fluxes.xlsx'


def ecoli_steady_state_simulation():

    model = Model('ecoli')
    model.read_from_file(MODEL_FILE)
    
    with model.simulator('ss') as sim:
        sim.set_target_EMUs({
            'Ala': ['123', [2,3]], 
            'Glu': [1,2,3,4,5], 
            'Phe': [[1,2,3,4,5], '123456789']
        })
        
        # specify the lableing strategy
        sim.set_labeling_strategy(
            'Glc.ex', 
            labeling_pattern = ['100000', '111111'], 
            percentage = [0.754, 0.246], 
            purity = [0.997, 0.994]
        )   # call this method for each labeled substrate
        
        # read the flux distribution
        sim.set_fluxes_from_file(FLUXES)
        
        # simulate MDVs
        sim.prepare(n_jobs = 3)
        res = sim.simulate()
        
    print(res)




if __name__ == '__main__':

    ecoli_steady_state_simulation()    
