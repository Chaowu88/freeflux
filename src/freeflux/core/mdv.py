#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '02/18/2022'




import re
from functools import lru_cache
from collections import Iterable, OrderedDict, Counter
from itertools import combinations_with_replacement
from numbers import Real
import numpy as np
from scipy.special import comb
from scipy.linalg import inv


natAbuns = {'H': [0.999885, 0.000115],
            'C': [0.9893, 0.0107],
            'N': [0.99632, 0.00368],
            'O': [0.99757, 0.00038, 0.00205],
            'Si': [0.922297, 0.046832, 0.030872],
            'S': [0.9493, 0.0076, 0.0429, 0.0002]}
# ref. https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf            




class MDV():
    '''
    Define MDV (i.e., mass isotopomer distribution vector) object and its operations.

    Convolution between mdv1 and mdv2 can be performed in three ways: 
    1. mdv1.conv(mdv2);
    2. mdv1*mdv2;
    3. mdv.conv(mdv1, mdv2).
    
    The zero element for convolution is MDV([0]), and the identity element is MDV([1]).
    
    Scalar multiplication (a*mdv) and MDV addition (mdv1 + mdv2) are also supported.
    In these cases, the resulting MDV are not automatically normalized.
    
    In addition to C, MDV can be built based on H, N, O, Si, and S.

    Parameters
    ----------
    fractions: list or array
        MDV vector.
    nonnegative: bool
        Whether to keep the elements >= 0.
    normalize: bool
        Whether to normalize MDV vector to ensure the sum equals one.
    base_atom: str
        Base atom for MDV.
            
    Attributes
    ----------
    value: array
        MDV vector.
    n_atoms: int
        # of atom.
    base_atom: str
        Base atom for MDV.
    fl: float
        Fractional labeling.
    '''
    
    def __init__(self, fractions, nonnegative = True, normalize = True, base_atom = 'C'):
        '''
        Parameters
        ----------
        fractions: list or array
            MDV vector.
        nonnegative: bool
            Whether to keep the elements >= 0.
        normalize: bool
            Whether to normalize MDV vector to ensure the sum == 1.
        base_atom: str
            Base atom for MDV.
        '''
        
        self.value = np.array(fractions)
        
        if nonnegative:
            self.value[self.value < 0] = 0
        
        if normalize and self.value.sum() != 0:
            self.value = self.value / self.value.sum()
            
        self.base_atom = base_atom
    

    def __iter__(self):
        
        return self.value
        
    
    def __getitem__(self, key):
        '''
        Parameters
        ----------
        key: int or slice
            Index.
        '''
        
        return self.value[key]
        
        
    def __array__(self):
        
        return self.value
    
    
    def conv(self, mdv):
        '''
        Parameters
        ----------
        mdv: list or array or MDV 
        '''
        
        if not isinstance(mdv, self.__class__):
            mdv = self.__class__(mdv)
        
        mdvConv = gen_conv(self.value, mdv.value)
        
        '''
        if mdv.n_atoms > self.n_atoms:
            mdv1, mdv2 = mdv, self
        else:
            mdv1, mdv2 = self, mdv
        
        mdvConv = []
        for i in range(mdv1.n_atoms+mdv2.n_atoms+1):
            
            if i <= mdv2.n_atoms:
                mdvConv.append(sum([mdv1[i-j] * mdv2[j] for j in range(i+1)]))
            
            elif mdv2.n_atoms < i <= mdv1.n_atoms:
                mdvConv.append(sum([mdv1[i-j] * mdv2[j] for j in range(mdv2.n_atoms+1)]))
            
            else:
                mdvConv.append(sum([mdv1[i-j] * mdv2[j] for j in range(i-mdv1.n_atoms, mdv2.n_atoms+1)]))
        '''    
        return self.__class__(mdvConv)    
    
        
    def __mul__(self, other):
        '''
        Parameters
        ----------
        other: scalar, list, array or MDV
        '''
        
        if isinstance(other, Real):
            return self.__class__(other * self.value, nonnegative = False, normalize = False)                
        elif isinstance(other, Iterable):
            return self.conv(other)    
        else:
            return NotImplemented
            
    
    def __rmul__(self, other):
        '''
        Parameters
        ----------
        other: MDV
        '''
        
        return self * other
        
        
    def __add__(self, other):
        '''
        Parameters
        ----------
        other: MDV
        '''
        
        if isinstance(other, type(self)):
            return self.__class__(self.value + other.value, nonnegative = False, normalize = False)
        else:
            return NotImplemented
        
        
    def __radd__(self, other):
        '''
        Parameters
        ----------
        other: MDV
        '''
        
        return self + other
    
    
    def _correction_matrix(self, X, n_Xs):
        '''
        Parameters
        ----------
        X: str
            Which element the correction matrix will be generated for.
        n_Xs: int
            # of X atoms in metabolite or metabolite fragment, these atoms will be corrected.
        '''
        
        corrMat = np.zeros((self.n_atoms+1, self.n_atoms+1))
        
        natMDV = get_natural_MDV(n_Xs, base_atom = X) 
        
        for i in range(self.n_atoms+1):
            for j in range(i+1):    
                if i - j <= natMDV.n_atoms:
                    corrMat[i, j] = natMDV.value[i-j] 
                else:
                    corrMat[i, j] = 0
                    
        return corrMat
                
    
    def correct_for_natural_abundance(self, atom_dict):
        '''
        Parameters
        ----------
        atom_dict: dict
            element needs to be corrected => # of corresponding atoms in metabolite (fragment).
        '''
        
        corrMat = np.eye(self.n_atoms+1)
        for X, n_Xs in atom_dict.items():
            corrMatX = self._correction_matrix(X, n_Xs)
            corrMat = np.dot(corrMat, corrMatX)
        
        corrMDV = np.dot(inv(corrMat), self.value)

        return self.__class__(corrMDV)
        
        
    def correct_for_inoculum(self, fraction):
        '''
        Parameters
        ----------
        fraction: float in [0, 1]
            Fraction of inoculum in biomass measured.
        '''
        
        natMDV = get_natural_MDV(self.n_atoms, base_atom = self.base_atom)
        corrMDV =  (self.value - fraction * natMDV.value) / (1 - fraction)

        return self.__class__(corrMDV)
    
    
    @property
    @lru_cache()
    def n_atoms(self):
        
        return self.value.size - 1
        
        
    @property
    @lru_cache()
    def fl(self):
    
        if self.n_atoms == 0:
            raise ValueError('no atom in MDV')
        else:
            return round((self.value * np.arange(self.n_atoms+1)).sum() / self.n_atoms, 4)
        
        
    def __repr__(self):
    
        return '%s([%s])' % (self.__class__.__name__, ', '.join(map(str, self.value.round(3))))

    
def _isotopomer_combination(n_atoms, n_natural_isotops):
    '''
    Parameters
    ----------    
    n_atoms: int
        # of atoms.
    n_natural_isotops: int
        # of natural isotopmers.
    
    Returns
    -------
    combs2: dict
    
    Example
    -------
    >>> _isotopomer_combination(2, 3)
    OrderedDict([(0, [Counter({0: 2})]), (1, [Counter({0: 1, 1: 1})]), (2, [Counter({0: 1, 2: 1}), 
    Counter({1: 2})]), (3, [Counter({1: 1, 2: 1})]), (4, [Counter({2: 2})])])
    '''
    
    allCombos = combinations_with_replacement(range(n_natural_isotops), n_atoms)   
    
    combos1 = {}
    for combo in allCombos:
        combos1.setdefault(sum(combo), []).append(combo)
    
    combos2 = OrderedDict()
    for mass in sorted(combos1):
        for combo in combos1[mass]:        
            combos2.setdefault(mass, []).append(Counter(combo))
    
    return combos2
    
    
def get_natural_MDV(n_atoms, base_atom = 'C'):
    '''
    Calculate the MDV of an unlabeled fragment.

    Parameters
    ----------
    n_atoms: int
        # of atoms.
    base_atom: str
        Base atom for MDV.
    '''
    
    natAbun = natAbuns[base_atom]
    
    combos = _isotopomer_combination(n_atoms, len(natAbun))
    
    mdv = []
    for _, counters in combos.items():    
        ele = 0
        for counter in counters:
            item = 1
            left = n_atoms
            for isotop, count in counter.items():
                item *= comb(left, count) * natAbun[isotop]**count
                left -= count
            ele += item
        mdv.append(ele)
    
    return MDV(mdv)


def get_substrate_MDV(atom_nos, labeling_pattern, percentage, purity):
    '''
    Calculate the MDV of a fragment from labeled substrate. 
    Currently this function only supports C-labeled substrate.
    
    Parameters
    ----------
    atom_nos: list of int
        Atom NOs. 
    labeling_pattern: str or list of str
        Labeling pattern of substrate, '0' for unlabeled atom, '1' for labeled atom, 
        e.g., '100000' for 1-13C glucose. 
        
        List if tracer with multiple labeling patterns are used. 
        
        Natural substrate (with all '0's) don't need to be explicitly set.
        
        If str, labeling_pattern should not be natural substrate.
    percentage: float or list of float
        Molar percentage (in range of [0,1]) of corresponding tracer. 
        Sum of percentage should be <= 1, and the rest will be considered as natural substrate.
        
        List if tracer with multiple labeling patterns are used. 
        
        * If list, len(percentage) should be equal to len(labeling_pattern).
        * If float, labeling_pattern should not be natural substrate.
    purity: float or list of float
        Labeled atom purity (in range of [0,1]) of corresponding tracer.
            
        List if tracer with multiple labeling patterns are used.

        * If list, len(purity) should be equal to len(labeling_pattern).
        * If float, labeling_pattern should not be natural substrate.
    '''
    
    if not isinstance(labeling_pattern, list):
        
        if re.match(r'^0+$', labeling_pattern):
            raise ValueError('use mdv.get_natural_MDV for natural substrate')

        labeling_pattern = [labeling_pattern]
        
    if not isinstance(percentage, list):
        percentage = [percentage]
        
    if sum(percentage) > 1:
        raise ValueError('sum of tacer percentage should be <= 1')

    if not isinstance(purity, list):
        purity = [purity]

    
    baseAtom = 'C'
    
    nAtoms = len(atom_nos)
    
    singleAtomMdv = {}
    for pat, pur, in zip(labeling_pattern, purity):
        singleAtomMdv[pat] = {'1': MDV([1-pur, pur]), '0': get_natural_MDV(1, base_atom = baseAtom)}
    
    
    mdvAllTracer = MDV(np.zeros(nAtoms+1))
    for pat, per, pur in zip(labeling_pattern, percentage, purity):
        
        mdvPerTracer = MDV([1])
        for atomNO in atom_nos:
            mdvPerTracer *= singleAtomMdv[pat][pat[atomNO-1]]
        
        mdvAllTracer += per * mdvPerTracer     
    
    return mdvAllTracer + (1 - sum(percentage)) * get_natural_MDV(nAtoms, base_atom = baseAtom)
    
    
def gen_conv(arr1, arr2):
    '''
    Parameters
    ----------
    arr1, arr2: list or array
        If one of the arrs is 2-D, the function performs 2-D convolution of MDV (array) and MDV derivative.
    '''
    
    nAtoms1 = len(arr1) - 1
    nAtoms2 = len(arr2) - 1
    
    if nAtoms2 > nAtoms1:
        arr1, arr2 = arr2, arr1
        nAtoms1, nAtoms2 = nAtoms2, nAtoms1
    
    arrConv = []
    for i in range(nAtoms1+nAtoms2+1):
        if i <= nAtoms2:
            arrConv.append(sum([arr1[i-j] * arr2[j] for j in range(i+1)]))    
        elif nAtoms2 < i <= nAtoms1:
            arrConv.append(sum([arr1[i-j] * arr2[j] for j in range(nAtoms2+1)]))
        else:
            arrConv.append(sum([arr1[i-j] * arr2[j] for j in range(i-nAtoms1, nAtoms2+1)]))
        
    return np.array(arrConv)    
    
    
def conv(mdv1, mdv2):
    '''
    Perform convolution between two MDVs.

    Parameters
    ----------
    mdv1: list or array or MDV
    mdv2: list or array or MDV
    
    Returns
    -------
    mdv: MDV
    '''
    
    if not isinstance(mdv1, MDV):
        mdv1 = MDV(mdv1)
    
    if not isinstance(mdv2, MDV):
        mdv2 = MDV(mdv2)
    
    return mdv1 * mdv2
        

def diff_conv(mdv_mdvder1, mdv_mdvder2):
    '''
    Parameters
    ----------
    mdv_mdvder1: 2-list of arrays
        namely [MDV, MDVder], MDVder in shape of (len(MDV), len(v))
    mdv_mdvder2:: 2-list of arrays
        namely [MDV, MDVder], MDVder in shape of (len(MDV), len(v))
        
    Returns
    -------
    mdv: array
    mdvder: 2-D array
        In shape of (len(MDV), len(v)).
    '''
    
    mdv1, mdvder1 = mdv_mdvder1
    mdv2, mdvder2 = mdv_mdvder2
    
    if isinstance(mdv1, MDV):
        mdv1 = mdv1.value
        
    if isinstance(mdv2, MDV):
        mdv2 = mdv2.value
        
    mdv = gen_conv(mdv1, mdv2)
    mdvder = gen_conv(mdv1, mdvder2) + gen_conv(mdv2, mdvder1)
    
    return mdv, mdvder
