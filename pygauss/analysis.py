# -*- coding: utf-8 -*-
import os, glob
from itertools import product
import copy
import re

import pandas as pd
from pandas.tools.plotting import radviz
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans

from IPython.core.display import clear_output

from .molecule import Molecule
from .utils import df_to_img

import paramiko
import socket
import errno

class Analysis(object):
    def __init__(self, folderpath='', 
                 headers=[], ssh_server='', ssh_username='', ssh_passwrd=''):
        """analysis class 

        folderpath : str
            the folder directory storing the file to be analysed
        headers : list
            the variable categories for each computation
        """    
        self._folderpath = False
        self._ssh_server = False
        self._ssh_username = False
        #TODO encrypt?
        self._ssh_passwrd = False
        if folderpath: self.set_folderpath(folderpath, ssh_server, 
                                           ssh_username, ssh_passwrd)
            
        heads = headers[:]+['Molecule']
        
        self._df = pd.DataFrame(columns=heads)
        self._next_index = 0
        
    def __repr__(self):
        return self.get_table().to_string()
        
    def copy(self):
        clone = copy.deepcopy(self)
        return clone        

    def get_folderpath(self):
        if not self._folderpath:
            raise IOError('folder path not set')
        return self._folderpath, self._ssh_server, self._ssh_username

    def _connect_ssh(self, ssh_server, ssh_username, ssh_passwrd):
        """ connect and verify ssh connection """
        ssh = paramiko.SSHClient() 
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        try:
            ssh.connect(ssh_server, username=ssh_username, password=ssh_passwrd)
        except socket.error, e:
            raise IOError(
            'could not connect to the ssh server: \n {0} \n {1}'.format(ssh_server, e))
        except paramiko.ssh_exception.AuthenticationException, e:
            raise IOError(
            'username or password authentication error \n {0}'.format(e))
        except Exception, e:
            raise IOError('error connectin to server: \n {0}'.format(e))

        return ssh
        
    def set_folderpath(self, folderpath, ssh_server='', 
                       ssh_username='', ssh_passwrd=''): 
                   
        if not ssh_server:
            if not os.path.exists(str(folderpath)):
                raise IOError("{0} does not exist".format(folderpath))
            self._folderpath = folderpath
        else:
            if not ssh_passwrd:
                ssh_passwrd = raw_input('Please enter ssh server password:')
                try:
                    clear_output()    
                except:
                    pass

            ssh = self._connect_ssh(ssh_server, ssh_username, ssh_passwrd)

            sftp = ssh.open_sftp()
            try:
                sftp.stat(folderpath)
            except IOError, e:
                ssh.close()
                if e.errno == errno.ENOENT:
                    raise IOError("{0} does not exist on server: {1}".format(folderpath, 
                                                                          ssh_server))
                else:
                    IOError('error trying to validate folder \n {0}'.format(e))
            
            ssh.close()
            
            self._folderpath = folderpath
            self._ssh_server = ssh_server
            self._ssh_username = ssh_username
            self._ssh_passwrd = ssh_passwrd
            
    folderpath = property(get_folderpath, set_folderpath, 
                        doc="The folderpath for gaussian runs")      
    
    def list_files_in_dir(self, pattern=None):
        if not self._folderpath:
            raise IOError('no folder set')
        if self._ssh_server:

            ssh = self._connect_ssh(self._ssh_server, 
                              self._ssh_username, self._ssh_passwrd)

            sftp = ssh.open_sftp()
            files = sftp.listdir(self._folderpath)
            
            ssh.close()
            if pattern:
                pattern = "".join(
                [ c if c.isalnum() or c=='*' else "["+c+"]" for c in pattern]
                ).replace('*', '.*')
                files = filter(lambda x: re.match(pattern,x), files)
        else:
            if not pattern:
                pattern = '*'
            filepaths = glob.glob(os.path.join(self._folderpath, pattern))
            files = [os.path.basename(f) for f in filepaths]
        
        return files

    def add_run(self, identifiers={}, 
                      init_fname=None, opt_fname=None, 
                      freq_fname=None, nbo_fname=None,
                      add_if_error=False,
                      alignto=[], sftp=None):
        """add single Gaussian run input/outputs """             

        if not sftp and self._ssh_server:
    
            ssh = self._connect_ssh(self._ssh_server, 
                              self._ssh_username, self._ssh_passwrd)
                              
            mol_sftp = ssh.open_sftp()
        else:
            mol_sftp = sftp

        molecule = Molecule(self._folderpath,
                            init_fname, opt_fname, 
                            freq_fname, nbo_fname,
                            sftp=mol_sftp,
                            alignto=alignto,
                            fail_silently=True)
        
        num_files = filter(lambda x:x, [init_fname, opt_fname, 
                                        freq_fname, nbo_fname])
        if not len(molecule.get_init_read_errors()) == num_files:
            if not molecule.get_init_read_errors() or add_if_error:
                    
                identifiers['Molecule'] = molecule
                series = pd.DataFrame(identifiers, 
                                      index=[self._next_index])
                self._df = self._df.copy().append(series)#, ignore_index=True)
                self._next_index += 1
        
        if not sftp and self._ssh_server:
            ssh.close()
        
        return molecule.get_init_read_errors() 
             
    def add_runs(self, headers=[], values=[], 
                 init_pattern=None, opt_pattern=None, 
                 freq_pattern=None, nbo_pattern=None,
                 add_if_error=False,
                 alignto=[], ipython_print=False):
        """add multiple Gaussian run inputs/outputs """             
        if self._ssh_server:
    
            ssh = self._connect_ssh(self._ssh_server, 
                              self._ssh_username, self._ssh_passwrd)
                              
            sftp = ssh.open_sftp()
        else:
            sftp=None

        read_errors=[]
        for idents in product(*values):
            identifiers = dict(zip(headers, idents))
            if ipython_print: print identifiers
            init = init_pattern.format(*idents) if init_pattern else None
            if type(opt_pattern) is str:
                opt = opt_pattern.format(*idents) if opt_pattern else None
            elif type(opt_pattern) is list or type(opt_pattern) is tuple:
                opt = [o.format(*idents) for o in opt_pattern]
            else:
                opt = None
            freq = freq_pattern.format(*idents) if freq_pattern else None
            nbo = nbo_pattern.format(*idents) if nbo_pattern else None
            
            file_read_errs = self.add_run(identifiers, init, opt, freq, nbo,
                        add_if_error=add_if_error, alignto=alignto, sftp=sftp)
            
            for fname, msg in file_read_errs:
                idents = identifiers.copy()
                idents.pop('Molecule', '_')
                idents['File'] = fname
                idents['Error_Message'] = msg
                read_errors.append(idents)

            if ipython_print: 
                try:
                    clear_output(wait=True)    
                except:
                    pass
        
        if self._ssh_server:
            ssh.close()
                        
        err_df = pd.DataFrame(read_errors)
        if read_errors:
            cols = err_df.columns.tolist()
            cols.remove('File')
            cols.append('File')
            cols.remove('Error_Message')
            cols.append('Error_Message')
            err_df = err_df[cols]
        
        return err_df
                
    def get_table(self, rows=[], columns=[],  filters={},
                  precision=4, head=False, mol=False, 
                  row_index=[], column_index=[], 
                  as_image=False, na_rep='-', font_size=None,
                  width=None, height=None, unconfined=False):
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
        as_image : bool
            output the table as an image (used pygauss.utils.df_to_img)
        na_rep : str
            how to represent empty (nan) cells (if outputting image)
        width, height, unconfined : int, int, bool
            args for IPy Image
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
            if type(row_index) is list:
                columns += row_index
            else:
                columns.append(row_index)
            columns = list(set(columns))
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
            df = df.head(head)
        
        if as_image:
            return df_to_img(df, na_rep=na_rep, font_size=font_size,
                             width=width, height=height, unconfined=unconfined)            
            
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
        """returns a series of a basic run property or nan if it is not available

        prop = 'basis', 'nbasis', 'optimised' or 'conformer'        
        """
        if prop not in self._basic_properties.keys():
            raise ValueError('{0} not a molecule property'.format(prop))
        
        def get_prop(m):
            method = getattr(m, self._basic_properties[prop])
            try: 
                out = method()
            except:
                out = pd.np.nan
            return out
            
        return self._df.Molecule.map(get_prop)

    def add_basic_properties(self, props=['basis', 'nbasis', 
                                          'optimised', 'conformer']):
        """adds columns giving info of basic run properties """
        for prop in props:
            try:
                series = self.get_basic_property(prop)
            except Exception:
                print 'error reading {0} \n setting to NaN'.format(prop)
                series = pd.np.nan
            self._df[prop.capitalize()] = series  
        
        return self.get_table()
    
    def remove_non_optimised(self):
        """removes runs that were not optimised """
        non_optimised = self._df[self.get_basic_property('optimised')!=True].copy()
        self._df = self._df[self.get_basic_property('optimised')==True]
        return non_optimised
        
    def remove_non_conformers(self):
        """removes runs with negative frequencies """
        non_conformers = self._df[self.get_basic_property('conformer')!=True].copy()
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
                vals = df.Molecule.map(func)
                if n in self._df.columns:
                    self._df[n] = vals.combine_first(self._df[n])
                else:                
                    self._df[n] = vals
    
        else:
            func = lambda m: getattr(m, method)(*args, **kwargs)
            vals = df.Molecule.map(func)
            if name in self._df.columns:
                self._df[name] = vals.combine_first(self._df[name])
            else:                
                self._df[name] = vals
        
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
        
        for indx, mol in zip(df.index, df.Molecule):
            if align_to: mol.set_alignment_atoms(*align_to)
            if mtype == 'initial':
                yield indx, mol.show_initial(gbonds=gbonds, ball_stick=ball_stick, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length)
            elif mtype == 'optimised':
                yield indx, mol.show_optimisation(gbonds=gbonds, ball_stick=ball_stick, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length)
            elif mtype == 'nbo':
                yield indx, mol.show_nbo_charges(gbonds=gbonds, ball_stick=ball_stick, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length,
                                       relative=relative, 
                                       minval=minval, maxval=maxval)
            elif mtype == 'highlight':
                yield indx, mol.show_highlight_atoms(highlight, gbonds=gbonds, 
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
        
        columns = columns[:]
        if columns and category_column not in columns:
            if all(isinstance(item, int) for item in columns):
                columns.append(col_names.index(category_column))
            else:
                columns.append(category_column) 
            
        df = self.get_table(rows, columns, filters)

        f, ax = plt.subplots()
        ax = radviz(df, category_column, ax=ax, s=point_size)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.set_frame_on(False)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        
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
                                                                        
                                  
                                  
                                  
                                  
                                  
                                  