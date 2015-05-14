# -*- coding: utf-8 -*-
"""
Created on Thu May 14 18:55:52 2015

@author: cjs14
"""
from nose import tools as that
from nose_parameterized import parameterized

import pygauss as pg

class Test_Start(object):
    def setUp(self):
        self.folder = pg.get_test_folder()
        
    def test_no_args(self):
        that.assert_raises(TypeError, pg.Molecule)
        
    @parameterized(['init_fname', 'opt_fname', 'freq_fname', 'nbo_fname'])
    def test_bad_fname(self, ftype):
        that.assert_raises(IOError, pg.Molecule, self.folder, **{ftype:'abcxyz'})

    def test_good_fname(self):
        fnames={'init_fname':'CJS1_emim-cl_B_init.com',
                 'opt_fname':'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log',
                 'freq_fname':'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_freq_unfrz.log',
                 'nbo_fname':'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_pop-nbo-full-_unfrz.log'}        
        pg.Molecule(self.folder, **fnames)

    def test_wildcard_fname(self):
        fnames={'init_fname':'*_emim-cl_B_init.com',
                 'opt_fname':'*_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log',
                 'freq_fname':'*_emim-cl_B_6-311+g-d-p-_gd3bj_freq_unfrz.log',
                 'nbo_fname':'*_emim-cl_B_6-311+g-d-p-_gd3bj_pop-nbo-full-_unfrz.log'}        
        pg.Molecule(self.folder, **fnames)

class Test_Init(object):
    pass

class Test_Opt(object):
    """
File Extracts:

Standard basis: 6-311+G(d,p) (5D, 7F)
   272 basis functions,   429 primitive gaussians,   281 cartesian basis functions

SCF Done:  E(RB3LYP) =  -805.105260336     A.U. after   11 cycles

Item               Value     Threshold  Converged?
Maximum Force            0.000070     0.000450     YES
RMS     Force            0.000016     0.000300     YES
Maximum Displacement     0.001604     0.001800     YES
RMS     Displacement     0.000326     0.001200     YES  
    """
    def setUp(self):
        self.mol = pg.Molecule(pg.get_test_folder(),
        opt_fname='CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log') 
        
    def test_basis_descript(self):
       
        that.assert_equal(self.mol.get_basis_descript(), '6-311+G(d,p) (5D, 7F)')
        

    def test_basis_funcs(self):
        
         that.assert_equal(self.mol.get_basis_funcs(), 272)
        
    def test_is_optimised(self):
        
        that.assert_true(self.mol.is_optimised())

    def test_optimisation_E(self):
        
        that.assert_equal(self.mol.get_optimisation_E(units='hartree'), -805.105260336)       

class Test_Freq(object):
    pass

class Test_NBO(object):
    pass

    
if __name__=='__main__':
    import nose
    nose.run()
