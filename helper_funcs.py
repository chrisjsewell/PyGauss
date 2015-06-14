# -*- coding: utf-8 -*-
"""
Created on Tue May 19 03:08:21 2015

@author: chris sewell
"""
import numpy as np
import pygauss as pg

def check_opts_set(file_pattern, folder='', 
                   server=None, username=None, passwrd=None):
    analysis = pg.Analysis(folderpath=folder, 
                   server=server, username=username, passwrd=passwrd)
    files = analysis.folder.list_files(file_pattern)
    if not files:
        return None
        
    for file_name in files:
        analysis.add_run({'File':file_name}, 
                         opt_fname=file_name, add_if_error=True)
    
    analysis.add_basic_properties(['optimised'])
    analysis.add_basic_properties(['opt_error'])
    analysis.add_mol_property('E(au)', 'get_optimisation_E',
                              units='hartree')
    analysis.add_mol_property('allE', 'get_optimisation_E',
                              units='hartree', final=False)
    
    
    df = analysis.get_table(row_index='File')
    df['Steps'] = df.allE.map(lambda x: x.shape[0])
    df['Ediff'] = df.allE.map(lambda x: x[-1] - x[-2] if x.shape[0]>1 else np.nan)
    return df.drop('allE', 1)

def check_freqs_set(file_pattern, folder='', 
                   server=None, username=None, passwrd=None):
    analysis = pg.Analysis(folderpath=folder, 
                   server=server, username=username, passwrd=passwrd)
    files = analysis.folder.list_files(file_pattern)
    if not files:
        return None

    for file_name in files:
        analysis.add_run({'File':file_name}, 
                         freq_fname=file_name, add_if_error=True)
        analysis.add_basic_properties(['conformer'])
    
    df = analysis.get_table(row_index='File')
    return df
    
if __name__=='__main__':
    df = check_freqs_set('CJS1_coro-emim*P*+g*freq*log', '/work/cjs14/CJS1', 
    'login.cx1.hpc.ic.ac.uk', 'cjs14')
    print df