#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '05/23/2022'




from os.path import splitext, isfile
import re
import pandas as pd




def read_model_from_file(file):
    '''
    Parameters
    ----------
    file: file path
        tsv or excel file with reactions.
    '''
    
    ext = splitext(file)[1]
    if re.search(r'tsv', ext):
        data = pd.read_csv(file, sep = '\t', comment = '#', header = None, names = ['subs', 'pros', 'rev'], 
                           index_col = 0)
    elif re.search(r'xls', ext):
        data = pd.read_excel(file, comment = '#', header = None, names = ['subs', 'pros', 'rev'], 
                             index_col = 0)
    else:
        raise TypeError('can only read from .tsv or excel file')
    
    return data


def read_preset_values_from_file(file):
    '''
    Parameters
    ----------
    file: file path
        tsv or excel file.
    '''
    
    ext = splitext(file)[1]
    if re.search(r'tsv', ext):
        data = pd.read_csv(file, sep = '\t', comment = '#', header = None, names = ['value'], 
                           index_col = 0, squeeze = True)
    elif re.search(r'xls', ext):
        data = pd.read_excel(file, comment = '#', header = None, names = ['value'], 
                             index_col = 0, squeeze = True)
    else:
        raise TypeError('can only read from .tsv or excel file')
    
    return data
        

def read_initial_values(ini, ids):
    '''
    Parameters
    ----------
    ini: ser of file in .tsv or .xlsx
        Initial values.
    ids: list
        IDs of fluxes or concentrations in correct order.
    '''

    if isfile(ini):
        ini = read_preset_values_from_file(ini)
        ini = ini[ids]
    elif isinstance(ini, pd.Series):
        ini = ini[ids]
    else:
        raise ValueError('initial values should be in pd.Series or file')

    return ini    

        
def read_measurements_from_file(file, inst_data = False):
    '''
    Parameters
    ----------
    file: file path
        tsv or excel file.
    inst_data: bool
        If True, additional column of timepoints is used as index (multiindex).
    '''
    
    if inst_data:
        indexCols = [0, 1]
    else:
        indexCols = 0
    
    ext = splitext(file)[1]
    if re.search(r'tsv', ext):
        data = pd.read_csv(file, sep = '\t', comment = '#', header = None, names = ['mean', 'sd'], 
                           index_col = indexCols)                        
    elif re.search(r'xls', ext):
        data = pd.read_excel(file, header = 0, index_col = indexCols)
        data.columns = ['mean', 'sd']
    else:
        raise TypeError('can only read from .tsv or excel file')
        
    return data
