# -*- coding: utf-8 -*-
"""
Created on Thu May 14 18:55:52 2015

@author: cjs14

"""
from nose import tools as that
from nose_parameterized import parameterized, param

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

    def test_optimisation_energy(self):
        
        that.assert_equal(self.mol.get_opt_energy(units='hartree'), -805.105260336)

    def test_plot_optimisation(self):

        ax = self.mol.plot_opt_energy(units='hartree')

        data_line = ax.get_lines()[0]
        xd = data_line.get_xdata()
        yd = data_line.get_ydata()

        that.assert_equal(yd[-1], -805.105260336)

class Test_Freq(object):
    def setUp(self):
        self.mol = pg.Molecule(folder_obj=pg.get_test_folder(),
        freq_fname='CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_freq_unfrz.log')

    def test_is_conformer(self):
    
        that.assert_true(self.mol.is_conformer())

    def test_plot_freq(self):
    
        df = self.mol.get_freq_analysis()
        
        that.assert_greater_equal(0., df['Frequency ($cm^{-1}$)'].min())
        that.assert_greater_equal(0., df['IR Intensity ($km/mol$)'].min())

    def test_plot_freq(self):

        ax = self.mol.plot_freq_analysis()
        
    def test_zpe_correction(self):
        """zero-point error correction
        
        File extract:
        Zero-point correction=                           0.168038 (Hartree/Particle)
        """
        
        that.assert_equal(self.mol.get_zeropt_energy(units='hartree'), 0.168038)


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
        """test_group_energy
        File lines to filter:        
   9. CR (   1) O   7                / 33. RY*(   1) H   9                    0.75   19.70    0.108
  15. LP (   2) O   7                / 32. RY*(   1) H   8                    0.62    1.41    0.026       
        """
        that.assert_almost_equal(
            self.mol.calc_sopt_energy(atom_groups=[[7], [8,9]], eunits='kcal'),
            0.75+0.62)

    def test_hbond_energy(self):
        """test_hbond_energy
         File lines to filter:       
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

    def test_returns_init(self):
        self.mol.show_initial()
    
    def test_returns_init_with_ball_stick(self):
        self.mol.show_initial(represent='ball_stick')

    def test_returns_init_with_wire(self):
        self.mol.show_initial(represent='wire')

    def test_returns_init_with_rotations(self):
        self.mol.show_initial(rotations=[[0,0,0],[45,45,45]])

    def test_returns_init_with_axes(self):
        self.mol.show_initial(axis_length=1)

    def test_returns_opt(self):
        self.mol.show_optimisation()

    def test_returns_nbo_charges(self):
        self.mol.show_nbo_charges()

class Test_PES(object):
    
    @parameterized(['local_maxs', 'local_mins', 'global_min', 'global_max'])
    def test_plot_pes(self, img_pos):
        mol2 = pg.molecule.Molecule(folder_obj=pg.get_test_folder(),
                                    alignto=[3,2,1],
                                    pes_fname=['CJS_emim_6311_plus_d3_scan.log',
                                               'CJS_emim_6311_plus_d3_scan_bck.log'])
        ax = mol2.plot_pes_scans([1,4,9,10], rotation=[0,0,90], img_pos=img_pos, zoom=0.5)


def mock_isosurface(*args, **kwargs):
    return None

class Test_Orbitals(object):
    
    def setUp(self):
        self.mol = pg.Molecule(folder_obj=pg.get_test_folder(),
                          opt_fname='CJS1_emim-cl_F_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log',
                          nbo_fname='CJS1_emim-cl_F_6-311+g-d-p-_gd3bj_pop-nbo-full-_unfrz.log',
                          atom_groups={'emim':range(1,20), 'cl':[20]},
                          alignto=[3,2,1])

    def test_orbital_count(self):
        that.assert_equal(self.mol.get_orbital_count(), 272)

    def test_orbital_homo_lumo(self):
        that.assert_equal(self.mol.get_orbital_homo_lumo(), (39,40))

    @parameterized([param(1,-101.39532), 
                    param(10, -9.30938),
                    param(39, -0.19600),
                    param(40, -0.04310),
                    param(272, 215.85458),
                    ])
    def test_orbital_energies(self, orbital, energy):        
        that.assert_almost_equal(self.mol.get_orbital_energies(orbital, eunits='hartree')[0], 
                          energy)
        
    def test_plot_dos(self):
        ax = self.mol.plot_dos(lbound=-20, ubound=10)
        
    def test_yield_orbital_images_no_isos(self):
        pg.isosurface.get_isosurface = mock_isosurface
        self.mol.yield_orbital_images(1)

from io import BytesIO
class MockFolder(object):
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        return
    def write_file(self, *arg, **kwargs):
        return BytesIO()

class Test_Combine_Molecules(object):
    
    def setUp(self):
        self.mol1 = pg.Molecule(folder_obj=pg.get_test_folder(),
                          opt_fname='CJS1_emim-cl_F_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log',
                          atom_groups={'emim':range(1,20), 'cl':[20]},
                          alignto=[3,2,1])
        self.mol2 = pg.Molecule(folder_obj=pg.get_test_folder(),
                          opt_fname='CJS1_emim-cl_F_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log',
                          atom_groups={'emim':range(1,20), 'cl':[20]},
                          alignto=[3,2,1])

    def test_combine_all(self):
        self.mol1.combine_molecules(self.mol2)
    
    def test_combine_part(self):
        self.mol1.combine_molecules(self.mol2, self_atoms='cl', other_atoms='emim')

    def test_combine_transform(self):
        self.mol1.combine_molecules(self.mol2, 
                        self_rotation=[10, 10, 10], other_rotation=[20, 20, 20], 
                        self_transpose=[1, 1, 1], other_transpose=[0.5, 0.5, 0.5])

    def test_combine_to_file(self):
        
        self.mol1.combine_molecules(self.mol2, charge=0, multiplicity=1,
                                    out_name='out', folder_obj=MockFolder())
        
    
if __name__=='__main__':
    import nose
    nose.run()
