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
                
    @parameterized(['init_fname', 'opt_fname', 'freq_fname', 'nbo_fname'])
    def test_bad_fname(self, ftype):
        that.assert_raises(IOError, pg.Molecule, folder_obj=self.folder, **{ftype:'abcxyz'})

    def test_good_fname(self):
        fnames={'init_fname':'CJS1_emim-cl_B_init.com',
                 'opt_fname':'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log',
                 'freq_fname':'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_freq_unfrz.log',
                 'nbo_fname':'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_pop-nbo-full-_unfrz.log'}        
        pg.Molecule(folder_obj=self.folder, **fnames)

    def test_wildcard_fname(self):
        fnames={'init_fname':'*_emim-cl_B_init.com',
                 'opt_fname':'*_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log',
                 'freq_fname':'*_emim-cl_B_6-311+g-d-p-_gd3bj_freq_unfrz.log',
                 'nbo_fname':'*_emim-cl_B_6-311+g-d-p-_gd3bj_pop-nbo-full-_unfrz.log'}        
        pg.Molecule(folder_obj=self.folder, **fnames)

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
        self.mol = pg.Molecule(folder_obj=pg.get_test_folder(),
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
    """
File Extracts:

 Second Order Perturbation Theory Analysis of Fock Matrix in NBO Basis

     Threshold for printing:   0.50 kcal/mol
    (Intermolecular threshold: 0.05 kcal/mol)
                                                                              E(2)  E(j)-E(i) F(i,j)
         Donor NBO (i)                     Acceptor NBO (j)                 kcal/mol   a.u.    a.u. 
 ===================================================================================================

 within unit  1
   7. CR (   1) O   1                / 20. RY*(   1) H   2                    0.75   19.70    0.108
  11. LP (   2) O   1                / 21. RY*(   1) H   3                    0.62    1.41    0.026

 from unit  1 to unit  2
   1. BD (   1) O   1 - H   2        / 23. RY*(   2) O   4                    0.05    2.35    0.010
   1. BD (   1) O   1 - H   2        / 27. RY*(   1) H   6                    0.07    1.43    0.009
   1. BD (   1) O   1 - H   2        / 36. BD*(   1) O   4 - H   5            0.52    0.95    0.020
   1. BD (   1) O   1 - H   2        / 37. BD*(   1) O   4 - H   6            0.05    0.95    0.006
   2. BD (   1) O   1 - H   3        / 27. RY*(   1) H   6                    0.07    1.41    0.009
   2. BD (   1) O   1 - H   3        / 37. BD*(   1) O   4 - H   6            0.08    0.93    0.008
  11. LP (   2) O   1                / 27. RY*(   1) H   6                    0.07    1.36    0.009
  11. LP (   2) O   1                / 36. BD*(   1) O   4 - H   5            0.20    0.88    0.012
  11. LP (   2) O   1                / 37. BD*(   1) O   4 - H   6            0.10    0.88    0.008

 from unit  1 to unit  3
   1. BD (   1) O   1 - H   2        / 39. BD*(   1) O   7 - H   9            0.11    1.04    0.010

 from unit  2 to unit  1
   3. BD (   1) O   4 - H   5        / 20. RY*(   1) H   2                    0.06    1.66    0.009
   3. BD (   1) O   4 - H   5        / 34. BD*(   1) O   1 - H   2            0.36    1.11    0.018
   4. BD (   1) O   4 - H   6        / 16. RY*(   1) O   1                    0.08    3.32    0.015
   4. BD (   1) O   4 - H   6        / 20. RY*(   1) H   2                    0.27    1.66    0.019
   4. BD (   1) O   4 - H   6        / 34. BD*(   1) O   1 - H   2            0.91    1.11    0.029
   8. CR (   1) O   4                / 20. RY*(   1) H   2                    0.06   19.82    0.031
   8. CR (   1) O   4                / 34. BD*(   1) O   1 - H   2            0.26   19.27    0.064
  13. LP (   2) O   4                / 34. BD*(   1) O   1 - H   2            9.59    1.06    0.091

 within unit  2
  13. LP (   2) O   4                / 26. RY*(   1) H   5                    0.63    1.45    0.027
  13. LP (   2) O   4                / 27. RY*(   1) H   6                    0.63    1.45    0.027

 from unit  2 to unit  3
   3. BD (   1) O   4 - H   5        / 28. RY*(   1) O   7                    0.08    3.32    0.015
   3. BD (   1) O   4 - H   5        / 33. RY*(   1) H   9                    0.27    1.66    0.019
   3. BD (   1) O   4 - H   5        / 39. BD*(   1) O   7 - H   9            0.91    1.11    0.029
   4. BD (   1) O   4 - H   6        / 33. RY*(   1) H   9                    0.06    1.66    0.009
   4. BD (   1) O   4 - H   6        / 39. BD*(   1) O   7 - H   9            0.36    1.11    0.018
   8. CR (   1) O   4                / 33. RY*(   1) H   9                    0.06   19.82    0.031
   8. CR (   1) O   4                / 39. BD*(   1) O   7 - H   9            0.26   19.27    0.063
  13. LP (   2) O   4                / 39. BD*(   1) O   7 - H   9            9.49    1.06    0.090

 from unit  3 to unit  1
   6. BD (   1) O   7 - H   9        / 34. BD*(   1) O   1 - H   2            0.11    1.04    0.010

 from unit  3 to unit  2
   5. BD (   1) O   7 - H   8        / 26. RY*(   1) H   5                    0.08    1.41    0.009
   5. BD (   1) O   7 - H   8        / 36. BD*(   1) O   4 - H   5            0.08    0.93    0.008
   6. BD (   1) O   7 - H   9        / 23. RY*(   2) O   4                    0.05    2.35    0.010
   6. BD (   1) O   7 - H   9        / 26. RY*(   1) H   5                    0.07    1.43    0.009
   6. BD (   1) O   7 - H   9        / 36. BD*(   1) O   4 - H   5            0.06    0.95    0.006
   6. BD (   1) O   7 - H   9        / 37. BD*(   1) O   4 - H   6            0.52    0.95    0.020
  15. LP (   2) O   7                / 26. RY*(   1) H   5                    0.07    1.36    0.009
  15. LP (   2) O   7                / 36. BD*(   1) O   4 - H   5            0.10    0.88    0.009
  15. LP (   2) O   7                / 37. BD*(   1) O   4 - H   6            0.20    0.88    0.012

 within unit  3
   9. CR (   1) O   7                / 33. RY*(   1) H   9                    0.75   19.70    0.108
  15. LP (   2) O   7                / 32. RY*(   1) H   8                    0.62    1.41    0.026

    
    """
    def setUp(self):
        self.mol = pg.Molecule(folder_obj=pg.get_test_folder(),
                               nbo_fname='test_3h2o_pop.log') 
    
    def test_total_energy(self):
        total_e = sum([0.75, 0.62, 0.05, 0.07, 0.52, 0.05, 0.07, 0.08, 0.07, 
                       0.2, 0.1, 0.11, 0.06, 0.36, 0.08, 0.27, 0.91, 0.06, 0.26, 
                       9.59, 0.63, 0.63, 0.08, 0.27, 0.91, 0.06, 0.36, 0.06, 0.26, 
                       9.49, 0.11, 0.08, 0.08, 0.05, 0.07, 0.06, 0.52, 0.07, 0.1, 
                       0.2, 0.75, 0.62])
        that.assert_almost_equal(self.mol.calc_sopt_energy(eunits='kcal'), total_e)
    
    def test_group_energy(self):
        """
   9. CR (   1) O   7                / 33. RY*(   1) H   9                    0.75   19.70    0.108
  15. LP (   2) O   7                / 32. RY*(   1) H   8                    0.62    1.41    0.026       
        """
        that.assert_almost_equal(
            self.mol.calc_sopt_energy(atom_groups=[[7], [8,9]], eunits='kcal'),
            0.75+0.62)

    def test_hbond_energy(self):
        """
  11. LP (   2) O   1                / 36. BD*(   1) O   4 - H   5            0.20    0.88    0.012
  11. LP (   2) O   1                / 37. BD*(   1) O   4 - H   6            0.10    0.88    0.008
  13. LP (   2) O   4                / 34. BD*(   1) O   1 - H   2            9.59    1.06    0.091
  13. LP (   2) O   4                / 39. BD*(   1) O   7 - H   9            9.49    1.06    0.090
  15. LP (   2) O   7                / 36. BD*(   1) O   4 - H   5            0.10    0.88    0.009
  15. LP (   2) O   7                / 37. BD*(   1) O   4 - H   6            0.20    0.88    0.012
        """
        hbond_e = sum([0.20, 0.10, 9.59, 9.49, 0.10, 0.20])
        that.assert_almost_equal(self.mol.calc_hbond_energy(eunits='kcal'), hbond_e)
        
    def test_no_hbond_energy(self):
        total_e = sum([0.75, 0.62, 0.05, 0.07, 0.52, 0.05, 0.07, 0.08, 0.07, 
                       0.2, 0.1, 0.11, 0.06, 0.36, 0.08, 0.27, 0.91, 0.06, 0.26, 
                       9.59, 0.63, 0.63, 0.08, 0.27, 0.91, 0.06, 0.36, 0.06, 0.26, 
                       9.49, 0.11, 0.08, 0.08, 0.05, 0.07, 0.06, 0.52, 0.07, 0.1, 
                       0.2, 0.75, 0.62])
        hbond_e = sum([0.20, 0.10, 9.59, 9.49, 0.10, 0.20])
        that.assert_almost_equal(self.mol.calc_sopt_energy(eunits='kcal', no_hbonds=True), 
                                 total_e - hbond_e)

class Test_Images(object):
    def setUp(self):
        self.mol = pg.Molecule(folder_obj=pg.get_test_folder(),
                init_fname='CJS1_emim-cl_B_init.com',
                opt_fname='CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log',
                nbo_fname='CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_pop-nbo-full-_unfrz.log') 

    def test_returns_init_img(self):
        self.mol.show_initial()
    
    def test_returns_opt_img(self):
        self.mol.show_optimisation()

    def test_returns_nbo_charges_img(self):
        self.mol.show_nbo_charges()
    
if __name__=='__main__':
    import nose
    nose.run()
