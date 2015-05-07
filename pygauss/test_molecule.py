# -*- coding: utf-8 -*-
"""
Created on Tue May 05 22:56:25 2015

@author: chris sewell
"""
from IPython.display import display

import pygauss as pg
folder = pg.get_test_folder()

mol = pg.molecule.Molecule(folder,
                init_fname='CJS1_emim-cl_B_init.com', 
                opt_fname=['CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_difrz.log',
                           'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_difrz_err.log',
                           'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log'],
                freq_fname='CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_freq_unfrz.log',
                nbo_fname='CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_pop-nbo-full-_unfrz.log')

display(mol.show_initial())
display(mol.show_initial(ball_stick=True))

print(mol.is_optimised())
print(mol.is_conformer())
print(mol.get_optimisation_E(units='hartree'))
mol.plot_optimisation_E(units='hartree')
mol.plot_IRfreqs()

display(mol.show_optimisation(ball_stick=True, axis_length=0.4, 
                              rotations=[[0,0,90], [-90, 90, 0]]))
mol.set_alignment_atoms(3,2,1)
display(mol.show_optimisation(ball_stick=True, axis_length=0.4, 
                              rotations=[[0,0,90], [-90, 90, 0]]))
display(mol.show_highlight_atoms([[3, 2, 1], [20]], ball_stick=True, axis_length=0.4, 
                              rotations=[[0,0,90], [-90, 90, 0]]))
display(mol.show_highlight_atoms([[3, 2, 1], [20]], ball_stick=False, axis_length=0.4, 
                              rotations=[[0,0,90], [-90, 90, 0]]))

mol.plot_opt_trajectory(20, [3,2,1])

mol2 = pg.molecule.Molecule(folder, 
                           pes_fname=['CJS_emim_6311_plus_d3_scan.log', 
                                      'CJS_emim_6311_plus_d3_scan_bck.log'], 
                           alignto=[3,2,1])   

mol2.plot_pes_scans([1,4,9,10], rotation=[0,0,90], img_pos='local_maxs', zoom=0.5)