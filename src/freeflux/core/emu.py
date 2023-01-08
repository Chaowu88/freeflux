#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '02/26/2022'




from functools import lru_cache
from collections.abc import Iterable
from .metabolite import Metabolite




class EMU():
    '''
    Define EMU (i.e., elementary metabolite unit) object and its operations.

    EMUs in the same metabolite and with the same atom NOs are considered as identical, 
    while metabolites which they derived from could be different.
    
    EMUs can be compared based self.metabolite_id and self.atom_nos.
    EMU and iterable object of EMUs can also be compared. In this case EMU will be put into
    the same iterable object with single item, and comparison between two iterables are performed.
    
    Currently only binary equivalents are considered.
    
    Parameters
    ----------
    id: str
        EMU ID
    metabolite: Metabolite or str
        Which metabolite the EMU comes from.
    atom_nos: list of int or str
        Atom NOs, sorted by number.

    Attributes
    ----------
    id: str
        EMU ID
    metabolite: Metabolite
        Which metabolite the EMU comes from.    
    metabolite_id: str
        Metabolite ID.
    atom_nos: list of int
        Atom NOs, sorted by number.
    size: int
        Size of EMU.
    equivalent_atom_nos: None or list of int
        Equivalent atom NOs, sorted by number.   
    equivalent: EMU
        Equivalent of EMU.
    '''
    
    def __init__(self, id, metabolite, atom_nos):
        '''
        Parameters
        ----------
        id: str
            EMU ID
        metabolite: Metabolite or str
            Which metabolite the EMU comes from.
        atom_nos: list of int or str
            Atom NOs, sorted by number.
        '''
        
        self.id = id
        if isinstance(metabolite, Metabolite):
            self.metabolite = metabolite
            self.metabolite_id = self.metabolite.id
        elif isinstance(metabolite, str):
            self.metabolite = Metabolite(metabolite)
            self.metabolite_id = metabolite
        if isinstance(atom_nos, list):
            self.atom_nos = sorted(atom_nos)
        elif isinstance(atom_nos, str):
            self.atom_nos = sorted(list(map(int, atom_nos)))
        self.size = len(self.atom_nos)
        
    
    def __hash__(self):
        
        return hash(self.metabolite_id) + hash(sum(self.atom_nos))
        
        
    def __eq__(self, other):
        '''
        Parameters
        ----------
        other: EMU or iterable
        '''
        
        if isinstance(other, Iterable):
            return type(other)([self]) == other
        else:
            return self.metabolite_id == other.metabolite_id and self.atom_nos == other.atom_nos
        
        
    def __lt__(self, other):
        '''
        Parameters
        ----------
        other: EMU or iterable
        '''
        
        if isinstance(other, Iterable):
            return type(other)([self]) < other
        else:
            if self.metabolite_id != other.metabolite_id:
                return self.metabolite_id < other.metabolite_id
            else:
                return self.atom_nos < other.atom_nos
                
    
    def __gt__(self, other):
        '''
        Parameters
        ----------
        other: EMU or iterable
        '''
        
        if isinstance(other, Iterable):
            return type(other)([self]) > other
        else:
            if self.metabolite_id != other.metabolite_id:
                return self.metabolite_id > other.metabolite_id
            else:
                return self.atom_nos > other.atom_nos
        
    
    @property
    @lru_cache()
    def equivalent_atom_nos(self):
        '''
        Returns
        -------
        equivAtomNOs: list of int or None
            Equivalent atom NOs, sorted by number.
        '''
        
        if len(self.metabolite.atoms_info) == 1:
            return None

        else:
            refAtoms, equivAtoms = self.metabolite.atoms_info   # only binary equivalents considered
            
            mapping = dict(zip(refAtoms, range(1, len(refAtoms)+1)))
            
            equivAtomNOs = sorted([mapping[equivAtoms[no-1]] for no in self.atom_nos])
            
            if equivAtomNOs == self.atom_nos:
                return None
            else:
                return equivAtomNOs
                
                
    @property
    @lru_cache()
    def equivalent(self):
        '''
        Returns
        -------
        EMU: EMU
            Equivalent of current EMU.
        '''
        
        equivAtomNOs = self.equivalent_atom_nos
        
        if equivAtomNOs:
            id = self.metabolite_id + '_' + ''.join(map(str, self.equivalent_atom_nos))
            metab = self.metabolite
            return EMU(id, metab, equivAtomNOs)
        else:
            return None
            
        
    def __repr__(self):
        
        return '%s %s_%s' % (self.__class__.__name__, self.metabolite_id, ''.join(map(str, self.atom_nos)))
        