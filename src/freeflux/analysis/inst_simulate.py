#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '06/02/2022'




from collections.abc import Iterable
from functools import partial
from ..core.mdv import MDV
from ..analysis.simulate import Simulator
from ..io.inputs import read_preset_values_from_file 
from ..io.results import InstSimResults




class InstSimulator(Simulator):
    '''
    Provided fluxes should be in the unit of umol/gCDW/s if concentrations in the unit of 
    umol/gCDW and timepoints in the unit of s.
    '''
    
    def set_concentration(self, metabid, value):
        '''
        Set metabolite concentration in unit of umol/gCDW.
        
        Parameters
        ----------
        metabid: str
            Metabolite ID.
        value: float
            Intracellular concentration.
        '''
        
        self.model.concentrations[metabid] = value
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_concentrations, metabid))
        
        
    def set_concentrations_from_file(self, file):
        '''
        Read metabolite concentrations (umol/gCDW) from file.
        
        Parameters
        ----------
        file: file path
            tsv or excel file with fields "metabolite_ID" and "value".
        '''
        
        concs = read_preset_values_from_file(file)
                
        for metabid, value in concs.iteritems():
            self.model.concentrations[metabid] = value
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(partial(self._unset_concentrations, concs.index.tolist()))    

    
    def _unset_concentrations(self, metabids):
        '''
        Parameters
        ----------
        metabids: str or list of str
            Metabolite ID(s).
        '''
        
        if not isinstance(metabids, Iterable):
            metabids = [metabids]
        
        for metabid in metabids: 
            if metabid in self.model.concentrations:
                self.model.concentrations.drop(metabids, inplace = True)
                
                
    def set_timepoints(self, timepoints):
        '''
        Set timepoints in unit of s.
        
        Parameters
        timepoints: list of float
            Timepoints when MDVs are simulated.
        '''
        
        self.model.timepoints = sorted(set(self.model.timepoints + timepoints))
        if 0 not in self.model.timepoints:
            self.model.timepoints = [0] + self.model.timepoints
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_timepoints)
        
    
    def _unset_timepoints(self):
        
        self.model.timepoints.clear()
    
    
    def _lambdify_matrix_Ms(self):
        
        if not self.model.matrix_Ms:
            self.calculator._lambdify_matrix_Ms()
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_matrix_Ms)
        
        
    def _unset_matrix_Ms(self):
        
        self.model.matrix_Ms.clear()
            
            
    def _calculate_initial_matrix_Xs(self):
        
        if not self.model.initial_matrix_Xs:
            self.calculator._calculate_initial_matrix_Xs()
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_initial_matrix_Xs)
        
        
    def _unset_initial_matrix_Xs(self):
        
        self.model.initial_matrix_Xs.clear()
        
        
    def _calculate_initial_matrix_Ys(self):
        
        if not self.model.initial_matrix_Ys:
            self.calculator._calculate_initial_matrix_Ys()
        
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_initial_matrix_Ys)
            
            
    def _unset_initial_matrix_Ys(self):
        
        self.model.initial_matrix_Ys.clear()
        
        
    def _build_initial_sim_MDVs(self):
        
        if not self.model.initial_sim_MDVs:
            self.calculator._build_initial_sim_MDVs()
            
        if self.contexts:
            context = self.contexts[-1]
            context.add_undo(self._unset_initial_sim_MDVs)
            
            
    def _unset_initial_sim_MDVs(self):
        
        self.model.initial_sim_MDVs.clear()
        
    
    def prepare(self, n_jobs = 1):
        '''
        Parameters
        ----------
        n_jobs: int
            If n_jobs > 1, decomposition job will run in parallel.
        '''
        
        self._decompose_network(n_jobs, lump = False)
        self._lambdify_matrix_As_and_Bs()
        self._lambdify_matrix_Ms()        
        self._calculate_substrate_MDVs()
        self._calculate_initial_matrix_Xs()
        self._calculate_initial_matrix_Ys()
        self._build_initial_sim_MDVs()
        

    def _check_dependencies(self):

        if self.model.total_fluxes.empty:
            raise ValueError('call set_flux or set_fluxes_from_file first')
        if self.model.concentrations.empty:
            raise ValueError('call set_concentration or set_concentrations_from_file first')
        if not self.model.timepoints:
            raise ValueError('call set_timepoints first')
        if not self.model.target_EMUs:
            raise ValueError('call set_target_EMUs first')    
        if not self.model.labeling_strategy:
            raise ValueError('call labeling_strategy first')
        
        if any([not self.model.EAMs, 
                not self.model.substrate_MDVs, 
                not self.model.initial_matrix_Xs, 
                not self.model.matrix_As, 
                not self.model.matrix_Bs, 
                not self.model.matrix_Ms,
                not self.model.timepoints]):
            raise ValueError('call prepare first')


    def simulate(self):
        '''
        Returns
        -------
        instTargetMDVs: dict of dict
           EMU IDs => {measurement timepoints => simulated MDVs}
        '''
        
        simInstMDVs = self.calculator._calculate_inst_MDVs()
        
        targetInstMDVs = {}
        for emuid in self.model.target_EMUs:
            instMDVs = self.model.initial_sim_MDVs[emuid].copy()
            instMDVs.update({t: MDV(simInstMDVs[emuid][t]) for t in simInstMDVs[emuid]})
            targetInstMDVs[emuid] = instMDVs
        
        return InstSimResults(targetInstMDVs)
