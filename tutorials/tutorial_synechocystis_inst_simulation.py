'''Example of MDV simulation at isotopically nonstationary (INST) state with a Synechocystis 
model.
'''


from os import makedirs
from freeflux import Model


MODEL_FILE = '../models/synechocystis/synthetic_data/reactions.xlsx' 
FLUXES = '../models/synechocystis/synthetic_data/fluxes.xlsx'
CONCS = '../models/synechocystis/synthetic_data/concentrations.xlsx'
OUT_DIR = '../results/synechocystis/synthetic_data/inst_simulation'


def syn_inst_simulation():
    
    model = Model('syn')
    model.read_from_file(MODEL_FILE)
    
    with model.simulator('inst') as isim:
        isim.set_target_EMUs({
            'G3P': ['23', '123'], 
            'DHAP': '123', 
            'PEP': '123', 
            'Fum': '1234', 
            'R5P': '12345', 
            'RuBP': '12345',
            'Cit': ['12345', '123456'],
            'S7P': '1234567'
        })
        
        # enter the timepoints when MDVs will be simulated
        isim.set_timepoints([0, 10, 30, 60, 120, 240])
        
        # specify the lableing strategy
        isim.set_labeling_strategy(
            'CO2.ex', 
            labeling_pattern = ['1'], 
            percentage = [0.5], 
            purity = [0.997]
        )   # call this method for each labeled substrate
        
        # read the flux distribution and concentrations
        isim.set_fluxes_from_file(FLUXES)
        isim.set_concentrations_from_file(CONCS)
        
        # simulate MDVs
        isim.prepare(n_jobs = 3)
        res = isim.simulate()
    
    for emuid in ['G3P_23', 'G3P_123', 'DHAP_123', 'PEP_123', 'Fum_1234', 'R5P_12345', 
                  'RuBP_12345', 'Cit_12345', 'Cit_123456', 'S7P_1234567']:
        res.plot_MDV_kinetics(emuid, show_fig = False, output_dir = OUT_DIR)




if __name__ == '__main__':

    makedirs(OUT_DIR, exist_ok = True)
    syn_inst_simulation()    
