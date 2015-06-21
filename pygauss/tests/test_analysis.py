# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 10:55:57 2015

@author: chris
"""

from nose import tools as that
from nose_parameterized import parameterized, param

import pygauss as pg

class Test_Analysis(object):
    def setUp(self):
        self.folder = pg.get_test_folder()
    
    def test_initialises(self):
        pg.Analysis(folder_obj=self.folder)

    def test_add_run(self):

        analysis = pg.Analysis(folder_obj=self.folder)
                              
        errors = analysis.add_run(identifiers={'Test':'test'},
                    init_fname='*emim-cl_B_init.com',
                    opt_fname='*emim-cl_B_6-311+g-d-p-_gd3bj_opt*unfrz.log',
                    freq_fname='*emim-cl_B_6-311+g-d-p-_gd3bj_freq*.log',
                    nbo_fname='*emim-cl_B_6-311+g-d-p-_gd3bj_pop-nbo-full-*.log',
                    alignto=[3,2,1], atom_groups={'emim':range(20), 'cl':[20]})
        that.assert_equals(len(errors), 0)
        that.assert_equals(analysis.count_runs(), 1)
        

    def test_add_run_fails(self):

        analysis = pg.Analysis(folder_obj=self.folder)
                              
        errors = analysis.add_run(identifiers={'Test':'test'},
                    init_fname='none.com',
                    opt_fname='*emim-cl_B_6-311+g-d-p-_gd3bj_opt*unfrz.log',
                    freq_fname='*emim-cl_B_6-311+g-d-p-_gd3bj_freq*.log',
                    nbo_fname='*emim-cl_B_6-311+g-d-p-_gd3bj_pop-nbo-full-*.log',
                    alignto=[3,2,1], atom_groups={'emim':range(20), 'cl':[20]})
        that.assert_equals(len(errors), 1)
        that.assert_equals(analysis.count_runs(), 0)

    def test_add_run_on_error(self):

        analysis = pg.Analysis(folder_obj=self.folder)
                              
        errors = analysis.add_run(identifiers={'Test':'test'},
                    init_fname='none.com',
                    opt_fname='*emim-cl_B_6-311+g-d-p-_gd3bj_opt*unfrz.log',
                    freq_fname='*emim-cl_B_6-311+g-d-p-_gd3bj_freq*.log',
                    nbo_fname='*emim-cl_B_6-311+g-d-p-_gd3bj_pop-nbo-full-*.log',
                    alignto=[3,2,1], atom_groups={'emim':range(20), 'cl':[20]},
                    add_if_error=True)
        that.assert_equals(len(errors), 1)
        that.assert_equals(analysis.count_runs(), 1)

    def test_add_runs(self):

        analysis = pg.Analysis(folder_obj=self.folder)
        analysis.add_runs(headers=['Cation', 'Anion', 'Initial'], 
                                       values=[['emim'], ['cl'],
                                               ['B', 'F', 'FE']],
                    init_pattern='*{0}-{1}_{2}_init.com',
                    opt_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_opt*unfrz.log',
                    freq_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_freq*.log',
                    nbo_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_pop-nbo-full-*.log',
                    alignto=[3,2,1], atom_groups={'emim':range(20), 'cl':[20]})
        
        that.assert_equals(analysis.count_runs(), 3)


    @parameterized([param(False, 1),
                    param(True, 2)])
    def test_add_runs_with_errors(self, add_if_error, count):

        analysis = pg.Analysis(folder_obj=self.folder)
        errors = analysis.add_runs(headers=['Cation', 'Anion', 'Initial'], 
                                       values=[['emim'], ['cl'],
                                               ['B', 'Error']],
                    init_pattern='*{0}-{1}_{2}_init.com',
                    opt_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_opt*unfrz.log',
                    freq_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_freq*.log',
                    nbo_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_pop-nbo-full-*.log',
                    add_if_error=add_if_error)
        
        that.assert_equals(analysis.count_runs(), count)

    @parameterized([param('name','get_opt_energy'),
                    param('name','calc_bond_angle', [1, 4, 9]),
                    param('name','calc_dihedral_angle', [1, 4, 9, 10]),
                    param('name','calc_min_dist', 'emim', 'cl'),
                    param('name','calc_2plane_angle', [1,2,3], [8,9,10]),
                    param('name','calc_nbo_charge', 'emim'),
                    param('name', 'calc_hbond_energy'),
                    param(['name', 'name', 'name'], 'calc_polar_coords_from_plane', 3, 2, 1, 20)
                  ])
    def test_add_mol_property(self, name, prop, *args, **kwargs):
        analysis = pg.Analysis(folder_obj=self.folder)
        analysis.add_runs(headers=['Cation', 'Anion', 'Initial'],
                          values=[['emim'], ['cl'],
                                  ['B', 'F', 'FE']],
                          init_pattern='*{0}-{1}_{2}_init.com',
                          opt_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_opt*unfrz.log',
                          freq_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_freq*.log',
                          nbo_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_pop-nbo-full-*.log',
                          alignto=[3,2,1], atom_groups={'emim':range(20), 'cl':[20]})

        analysis.add_mol_property(name, prop, *args, **kwargs)

    @parameterized(['initial', 'optimised', 'highlight', 'nbo', 'sopt', 'hbond'])
    def test_tbl_images(self, mtype):

        analysis = pg.Analysis(folder_obj=self.folder)
        analysis.add_runs(headers=['Cation', 'Anion', 'Initial'], 
                                       values=[['emim'], ['cl'],
                                               ['B', 'F', 'FE']],
                    init_pattern='*{0}-{1}_{2}_init.com',
                    opt_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_opt*unfrz.log',
                    freq_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_freq*.log',
                    nbo_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_pop-nbo-full-*.log',
                    alignto=[3,2,1], atom_groups={'emim':range(20), 'cl':[20]})

        fig, caption = analysis.plot_mol_images(mtype=mtype, max_cols=2,
                                info_columns=['Cation', 'Anion', 'Initial'],
                                rotations=[[0,0,90]])
    
    @parameterized([param('energy', True),
                    param('energy', False),
                    param('freq', True),
                    param('freq', False),
                    param('dos', False, atom_groups=['cl'], group_colors=['blue'], 
                          group_labels=['Cl'], group_fill=True, lbound=-20, ubound=10),
                  ])
    def test_tbl_graphs(self, gtype, share_plot, **kwargs):

        analysis = pg.Analysis(folder_obj=pg.get_test_folder())
        analysis.add_runs(headers=['Cation', 'Anion', 'Initial'], 
                               values=[['emim'], ['cl'],
                                       ['B', 'F', 'FE']],
            init_pattern='*{0}-{1}_{2}_init.com',
            opt_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_opt*unfrz.log',
            freq_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_freq*.log',
            nbo_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_pop-nbo-full-*.log',
            alignto=[3,2,1], atom_groups={'emim':range(20), 'cl':[20]})
            
        fig, caption = analysis.plot_mol_graphs(gtype=gtype, share_plot=share_plot, 
                                                info_columns=['Cation', 'Anion', 'Initial'], 
                                                max_cols=2, info_incl_id=True, tick_rotation=45, 
                                                **kwargs)
        
        that.assert_equal(type(caption), str)
