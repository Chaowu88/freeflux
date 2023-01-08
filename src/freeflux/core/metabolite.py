#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '02/16/2022'




import re




class Metabolite():
    '''
    Define Metabolite object which constitutes a Reaction.

    Metabolites are considered identical when they have the same IDs and atoms.
    
    Parameters
    ----------
    id: str
        Metabolite ID.
    atoms: str or list of str or None
        Carbons in metabolite, e.g., 'abcd' for metabolite without equivalents, 
        and ['abcd', 'dcba'] for metabolite with equivalents.

    Attributes
    ----------
    id: str
        Metabolite ID.
    atoms_info: dict
        equivalent => coefficient, e.g., {'abcd': 1.0} or {'abcd': 0.5, 'dcba': 0.5}.
    n_carbons: int
        # of carbons in metabolite. 
    host_reactions: set of Reaction or None
        Reactions hosting the metabolite.
    '''
    
    def __init__(self, id, atoms = None):
        '''
        Parameters
        ----------
        id: str
            Metabolite ID.
        atoms: str or list of str or None
            Carbons in metabolite, e.g., 'abcd' for metabolite without equivalents, 
            and ['abcd', 'dcba'] for metabolite with equivalents.
        '''
        
        self.id = id
        self.atoms = atoms
        if self.atoms and isinstance(self.atoms, str):
            self.atoms = [self.atoms]
        self.host_reactions = None
        
    
    def __hash__(self):
        
        if self.atoms:
            return hash(self.id) + sum([hash(atoms) for atoms in self.atoms])
        
        else:
            return hash(self.id)
            
        
    def __eq__(self, other):
        '''
        Parameters
        ----------
        other: Metabolite
        '''
        
        if self.atoms:
            return self.id == other.id and set(self.atoms) == set(other.atoms)
        
        else:
            if other.atoms:
                return False
            else:
                return self.id == other.id
        
    
    @property
    def atoms_info(self):
    
        if self.atoms:
            return {atoms: 1/len(self.atoms) for atoms in self.atoms}
        
        else:
            return None
    
    
    @property
    def n_carbons(self):
        
        if self.atoms:
            return len(re.search(r'[a-z]+', self.atoms).group())
        
        else:
            return None
            
            
    def __repr__(self):
        
        return '%s %s%s' % (self.__class__.__name__, self.id, '(' + ','.join(self.atoms) + ')' if self.atoms else '')
