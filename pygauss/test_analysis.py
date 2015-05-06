# -*- coding: utf-8 -*-
"""
Created on Fri May 01 04:11:16 2015

@author: chris
"""
import pygauss as pg
folder = pg.get_test_folder()

analysis = pg.analysis.Analysis(folder)
analysis.add_run({'Cation':'emim'},
                init_fname='CJS1_emim_-_init.com', 
                opt_fname='CJS1_emim_-_6-311+g-d-p-_gd3bj_opt_.log',
                freq_fname='CJS1_emim_-_6-311+g-d-p-_gd3bj_freq_.log',
                nbo_fname='CJS1_emim_-_6-311+g-d-p-_gd3bj_pop-nbo-full-_.log')

df, errors = analysis.add_runs(headers=['Cation', 'Anion', 'Initial'], 
                               values=[['emim'], ['cl'],
                                       ['B', 'BE', 'BM', 'F', 'FE', 'FM', 'T']],
            init_pattern='CJS1_{0}-{1}_{2}_init.com',
            opt_pattern='CJS1_{0}-{1}_{2}_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log',
            freq_pattern='CJS1_{0}-{1}_{2}_6-311+g-d-p-_gd3bj_freq_unfrz.log',
            nbo_pattern='CJS1_{0}-{1}_{2}_6-311+g-d-p-_gd3bj_pop-nbo-full-_unfrz.log')

print 'Read Errors:', errors

analysis.add_basic_properties()

analysis.add_mol_property('Energy (au)', 'get_optimisation_E', units='hartree')
analysis.add_mol_property('Cation chain, $\\psi$', 'calc_dihedral_angle', [1, 4, 9, 10])
analysis.add_mol_property('Cation Charge', 'calc_nbo_charge', range(1, 20))
analysis.add_mol_property(['Cation Charge center, $r$', 'Cation Charge center, $\\theta$', 
                           'Cation Charge center, $\\phi$'], 
                               'calc_nbo_charge_center', 3, 2, 1, atoms=range(1, 20))
analysis.add_mol_property('Anion Charge', 'calc_nbo_charge', [20])
analysis.add_mol_property(['Anion-Cation, $r$', 'Anion-Cation, $\\theta$', 'Anion-Cation, $\\phi$'], 
                               'calc_polar_coords_from_plane', 3, 2, 1, 20)
analysis.add_mol_property('Anion-Cation, $d_{min}$', 'calc_min_dist', range(1, 20), [20])

analysis.add_mol_property_subset('Energy (eV)', 'get_optimisation_E', rows=2, kwargs={'units':'eV'})

#print analysis.get_table() 
from IPython.display import display
mols = analysis.yield_mol_images(mtype='optimised',
                                      align_to=[3,2,1], 
                                      rotations=[[0, 0, 90], [-90, 90, 0]],
                                      axis_length=0.5)
for mol in mols:
    display(mol)
    
analysis.plot_radviz_comparison('Anion', columns=range(7, 16), rows=range(1, 6))

print analysis.calc_kmean_groups('Anion', 'cl', 4, columns=range(7, 16), rows=range(1, 6))                              

