# -*- coding: utf-8 -*-
"""
Created on Tue May 05 22:56:25 2015

@author: chris
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

display(mol.show_optimisation(ball_stick=True, axis_length=0.4, 
                              rotations=[[0,0,90], [-90, 90, 0]]))
mol.set_alignment_atoms(3,2,1)
display(mol.show_optimisation(ball_stick=True, axis_length=0.4, 
                              rotations=[[0,0,90], [-90, 90, 0]]))

mol.plot_opt_trajectory(20, [3,2,1])