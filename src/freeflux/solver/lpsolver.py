'''Define the FBAModel class.'''


__author__ = 'Chao Wu'


from pyomo.environ import (ConcreteModel, Var, Objective, Constraint, 
                           SolverFactory, maximize, minimize, value)
                                                   
                           
class FBAModel():
    
    def __init__(self):
        
        self.model = ConcreteModel()
        
        
    def build_flux_variables(self, fluxids, flux_bounds):
        '''
        Parameters
        ----------
        fluxids: list
            Flux IDs.
        flux_bounds: dict
            Reaction ID => [lower bound, upper bound].
        '''
        
        if set(fluxids) != flux_bounds.keys():
            raise ValueError('flux IDs not consistent with those with bounds')

        self.fluxids = fluxids
        
        def flux_bounds_rule(model, fluxid):
            return flux_bounds[fluxid]
            
        self.model.fluxes = Var(self.fluxids, bounds = flux_bounds_rule)

    
    def build_objective(self, objective, direction):
        '''
        Parameters
        ----------
        objective: dict
            Reaction ID => coefficient, i.e., objective function.
        direction: {"max", "min"}
            Optimization direction.
        '''
        
        if direction == 'max':
            sense = maximize
        elif direction == 'min':
            sense = minimize
        
        def obj_rule(model):
            objExpr = [coe*model.fluxes[rxnid] for rxnid, coe in objective.items()]
            return sum(objExpr)
        
        self.model.obj = Objective(rule = obj_rule, sense = sense)
        
        
    def build_mass_balance_constraints(self, stoy_mat):
        '''
        Parameters
        ----------
        stoy_mat: df
            Stoichiometric matrix.
        '''
        
        def mb_rule(model, metabid):
            fluxesExpr = [stoy_mat.loc[metabid, rxnid]*model.fluxes[rxnid] 
                          for rxnid in self.fluxids]
            return sum(fluxesExpr) == 0
            
        self.model.MBcstrs = Constraint(stoy_mat.index.tolist(), rule = mb_rule)
        
        
    def build_objective_constraint(self, objective, max_obj, gamma):
        '''
        Parameters
        ----------
        objective: dict
            Reaction ID => coefficient, i.e., objective function.
        max_obj: float
            Optimal objective.
        gamma: float
            A value in [0, 1].
        '''
        
        def objcstr_rule(model):
            objcstrExpr = [coe*model.fluxes[rxnid] for rxnid, coe in objective.items()]
            return sum(objcstrExpr) >= gamma*max_obj
            
        self.model.OBJcstr = Constraint(rule = objcstr_rule)
        
    
    def remove_objective(self):
        
        self.model.del_component(self.model.obj)
        
        
    def solve_flux(self):
        
        solver = SolverFactory('glpk')
        solver.solve(self.model)
        
        optObj = value(self.model.obj)
        
        optFluxes = {}
        for fluxid in self.fluxids:
            try:
                flux = value(self.model.fluxes[fluxid])
            except ValueError:
                print(f'flux value of {fluxid} not available')
                continue
            optFluxes[fluxid] = flux
        
        return optObj, optFluxes
        