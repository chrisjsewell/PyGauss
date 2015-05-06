# -*- coding: utf-8 -*-
import os
from itertools import product

import pandas as pd
from pandas.tools.plotting import radviz

from sklearn.cluster import KMeans

from IPython.core.display import clear_output

from .molecule import Molecule

class Analysis(object):
    def __init__(self, folderpath='', 
                 headers=[]):
        """analysis class """    
        self._folderpath = False
        if folderpath: self.set_folderpath(folderpath)
            
        heads = headers[:].append('Molecule')
        
        self._df = pd.DataFrame(columns=heads)
        self._next_index = 0
        
    def __repr__(self):
        return self.get_table().to_string()
        
    def get_folderpath(self):
        if not self._folderpath:
            raise IOError('folder path not set')
        return self._folderpath
    def set_folderpath(self, folderpath): 
        if not os.path.exists(str(folderpath)):
            raise ValueError("{0} does not exist".format(folderpath))
        self._folderpath = folderpath
        
    folderpath = property(get_folderpath, set_folderpath, 
                        doc="The folderpath for gaussian runs")        

    def add_run(self, identifiers={}, 
                      init_fname=None, opt_fname=None, 
                      freq_fname=None, nbo_fname=None,
                      alignto=[]):
        """add single Gaussian run input/outputs """             
        molecule = Molecule(self.folderpath, init_fname, opt_fname, 
                                            freq_fname, nbo_fname,
                                            alignto=alignto)
        identifiers['Molecule'] = molecule
        series = pd.DataFrame(identifiers, 
                              index=[self._next_index])
        self._df = self._df.copy().append(series)#, ignore_index=True)
        self._next_index += 1
        
        return self.get_table()
     
    def add_runs(self, headers=[], values=[], 
                 init_pattern=None, opt_pattern=None, 
                 freq_pattern=None, nbo_pattern=None,
                 alignto=[], ipython_print=False):
        """add multiple Gaussian run inputs/outputs """             
        read_errors=[]
        for idents in product(*values):
            identifiers = dict(zip(headers, idents))
            if ipython_print: print identifiers
            init = init_pattern.format(*idents) if init_pattern else None
            opt = opt_pattern.format(*idents) if opt_pattern else None
            freq = freq_pattern.format(*idents) if freq_pattern else None
            nbo = nbo_pattern.format(*idents) if nbo_pattern else None
            
            try:
                self.add_run(identifiers, init, opt, freq, nbo,
                             alignto=alignto)
            except IOError:
                read_errors.append(identifiers)

            if ipython_print: clear_output(wait=True)            
                        
        return self.get_table(), read_errors
                
    def get_table(self, rows=[], columns=[],  filters={},
                  precision=4, head=False, mol=False, 
                  row_index=[], column_index=[]):
        """return pandas table of requested data in requested format

        rows : integer or list of integers
            select row ids
        columns : string/integer or list of strings/integers
            select column names/positions
        filters : dict
            filter for rows with certain value(s) in specific columns
        precision : int
            decimal precision of displayed values
        head : int
            return only first n rows
        mol : bool
            include column containing the molecule objects
        row_index : string or list of strings
            columns to use as new index        
        column_index : list of strings
            srings to place in to higher order column indexs        
        """
        pd.set_option('precision', precision)
        
        if mol:
            df = self._df.copy()
        else:
            df = self._df.drop('Molecule', axis=1)
        
        for key, val in filters.iteritems():
            if type(val) is list or type(val) is tuple:
                 df = df[getattr(df, key).isin(val)]
            else:
                df = df[getattr(df, key)==val]
            
        if type(rows) is not list and type(rows) is not tuple:
            rows = [rows]
        if type(columns) is not list and type(columns) is not tuple:
            columns = [columns]
        if rows:
            df = df.loc[rows] 
        if columns:
            df = df.ix[:,columns]            
            
        if row_index: df = df.set_index(row_index) 

        if column_index:
            col_index=[]
            for col in df.columns:
                col_tuple = (' ', col)
                for term in column_index:
                    if len(col)>len(term):
                        if col[:len(term)] == term:
                            col_tuple = (term, col[len(term)+1:])
                            continue
                col_index.append(col_tuple)
            df.columns = pd.MultiIndex.from_tuples(col_index)
            
        if head:
            return df.head(head)
        else:
            return df
        
    def remove_rows(self, rows):
        
        self._df.drop(rows, inplace=True)

        return self.get_table()
        
    def remove_columns(self, columns):

        self._df.drop(columns, axis=1, inplace=True) 

        return self.get_table()

    _basic_properties={'nbasis':'get_basis_funcs',
                        'basis':'get_basis_descript',
                       'optimised':'is_optimised',
                       'conformer': 'is_conformer'}
                       
    def get_basic_property(self, prop):
        """returns a series of a basic run property

        prop = 'basis', 'nbasis', 'optimised' or 'conformer'        
        """
        if prop not in self._basic_properties.keys():
            raise ValueError('{0} not a molecule property'.format(prop))
        return self._df.Molecule.map(
            lambda m: getattr(m, self._basic_properties[prop])())

    def add_basic_properties(self, props=['basis', 'nbasis', 
                                          'optimised', 'conformer']):
        """adds columns giving info of basic run properties """
        for prop in props:
            series = self.get_basic_property(prop)
            self._df[prop.capitalize()] = series  
        
        return self.get_table()
    
    def remove_non_optimised(self):
        """removes runs that were not optimised """
        non_optimised = self._df[self.get_basic_property('optimised')==False].copy()
        self._df = self._df[self.get_basic_property('optimised')==True]
        return non_optimised
        
    def remove_non_conformers(self):
        """removes runs with negative frequencies """
        non_conformers = self._df[self.get_basic_property('conformer')==False].copy()
        self._df = self._df[self.get_basic_property('conformer')==True]
        return non_conformers

    def add_mol_property(self, name, method, *args, **kwargs):
        """compute molecule property for all rows and create data column """
                
        if type(name) is tuple or type(name) is list:
            for idx, n in enumerate(name):
                func = lambda m: getattr(m, method)(*args, **kwargs)[idx]
                self._df[n] = self._df.Molecule.map(func)
        else:
            func = lambda m: getattr(m, method)(*args, **kwargs)
            self._df[name] = self._df.Molecule.map(func)
        
        return self.get_table()

    def add_mol_property_subset(self, name, method, 
                                     rows=[], filters={}, 
                                     args=[], kwargs={}):
        """compute molecule property for a subset of rows and create/add-to data column """
        
        df = self.get_table(rows=rows, filters=filters, mol=True)

        if type(name) is tuple or type(name) is list:
            for idx, n in enumerate(name):
                func = lambda m: getattr(m, method)(*args, **kwargs)[idx]
                if n in self._df.columns:
                    self._df[n] = df.Molecule.map(func).combine_first(self._df[n])
                else:                
                    self._df[n] = df.Molecule.map(func)
    
        else:
            func = lambda m: getattr(m, method)(*args, **kwargs)
            if name in self._df.columns:
                self._df[name] = df.Molecule.map(func).combine_first(self._df[name])
            else:                
                self._df[name] = df.Molecule.map(func)
        
        return self.get_table()
            
    ## TODO will active work?
    def yield_mol_images(self, rows=[], filters={}, mtype='optimised',
                       align_to=[],
                        gbonds=True, ball_stick=True, rotations=[[0., 0., 0.]], 
                       zoom=1., width=300, height=300, axis_length=0,
                       relative=False, minval=-1, maxval=1,
                       highlight=[], active=False):
        """show molecules
        
        mtype = 'initial', 'optimised', 'nbo' or 'highlight'
        """
        df = self.get_table(columns=['Molecule'], rows=rows, 
                       filters=filters, mol=True)
        
        for mol in df.Molecule:
            if align_to: mol.set_alignment_atoms(*align_to)
            if mtype == 'initial':
                yield mol.show_initial(gbonds=gbonds, ball_stick=ball_stick, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length)
            elif mtype == 'optimised':
                yield mol.show_optimisation(gbonds=gbonds, ball_stick=ball_stick, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length)
            elif mtype == 'nbo':
                yield mol.show_nbo_charges(gbonds=gbonds, ball_stick=ball_stick, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length,
                                       relative=relative, 
                                       minval=minval, maxval=maxval)
            elif mtype == 'highlight':
                yield mol.show_highlight_atoms(highlight, gbonds=gbonds, 
                                               ball_stick=ball_stick, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length)
            else:
                raise ValueError('mtype must be initial, optimised, nbo or highlight')                
                
    def plot_radviz_comparison(self, category_column, 
                               columns=[], rows=[], filters={}, point_size=30):
        """return plot axis of radviz graph
        
        RadViz is a way of visualizing multi-variate data. 
        It is based on a simple spring tension minimization algorithm. 
        Basically you set up a bunch of points in a plane. In our case they are 
        equally spaced on a unit circle. Each point represents a single attribute. 
        You then pretend that each sample in the data set is attached to each 
        of these points by a spring, the stiffness of which is proportional to 
        the numerical value of that attribute (they are normalized to unit 
        interval). The point in the plane, where our sample settles to (where 
        the forces acting on our sample are at an equilibrium) is where a dot 
        representing our sample will be drawn. Depending on which class that 
        sample belongs it will be colored differently.
        """
        col_names = self._df.drop('Molecule', axis=1).columns.tolist()
        if category_column not in col_names:
            raise ValueError('{0} not in columns'.format(category_column))
        
        if columns and category_column not in columns:
            if all(isinstance(item, int) for item in columns):
                columns.append(col_names.index(category_column))
            else:
                columns.append(category_column) 
            
        df = self.get_table(rows, columns, filters)

        ax = radviz(df, category_column, s=point_size)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        #ax.set_frame_on(False)
        return ax
    
    def calc_kmean_groups(self, category_column, category_name, groups,
                    columns=[], rows=[], filters={}):
        """calculate the kmeans grouping of rows 
        
        The KMeans algorithm clusters data by trying to separate samples in n 
        groups of equal variance, minimizing a criterion known as the inertia 
        or within-cluster sum-of-squares. This algorithm requires the number of 
        clusters to be specified. It scales well to large number of samples and 
        has been used across a large range of application areas in many 
        different fields.
        """
        col_names = self._df.drop('Molecule', axis=1).columns.tolist()
        if category_column not in col_names:
            raise ValueError('{0} not in columns'.format(category_column))

        filters[category_column] = category_name
        df = self.get_table(rows, columns, filters)
        
        k_means = KMeans(n_clusters=groups)
        k_means.fit(df)
        cats = k_means.predict(df)
        
        return pd.DataFrame({'Name':category_name, 'Category':cats}, 
                            index=df.index)
                
if __name__ == '__main__':
    pass
                                                                        
                                  
                                  
                                  
                                  
                                  
                                  