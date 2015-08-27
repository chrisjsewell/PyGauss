# -*- coding: utf-8 -*-
from itertools import product, imap
import copy
import math
import string
import multiprocessing
import platform

import numpy as np
import pandas as pd
from pandas.tools.plotting import radviz
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from sklearn.cluster import KMeans

from IPython.core.display import clear_output

from .molecule import Molecule
from .utils import df_to_img
from .file_io import Folder

def unpack_and_make_molecule(val_dict):      
    if val_dict.has_key('args'):
        args = val_dict.pop('args')
    else:
        args = []            
    return Molecule(*args, **val_dict)

class Analysis(object):
    """a class to analyse multiple computations """
    
    def __init__(self, folderpath='', server=None, username=None, passwrd=None, 
                 folder_obj=None, headers=[]):
        """a class to analyse multiple computations

        Parameters
        ----------
        folderpath : str
            the folder directory storing the files to be analysed
        server : str
            the name of the server storing the files to be analysed
        username : str
            the username to connect to the server
        passwrd : str
            server password, if not present it will be asked for during initialisation
        headers : list
            the variable categories for each computation
            
        """  
        self._folder = None  
        if folder_obj:
            self._folder = folder_obj                              
        elif folderpath or server: 
            self.set_folder(folderpath, server, username, passwrd)
            
        heads = headers[:]+['Molecule']
        
        self._df = pd.DataFrame(columns=heads)
        self._df.index.name = 'ID'
        self._next_index = 0
        
    def __repr__(self):
        return self.get_table().to_string()
        
    def copy(self):
        clone = copy.deepcopy(self)
        return clone        
    
    def get_folder(self):
        return self._folder
        
    def set_folder(self, folderpath='', server=None, 
                       username=None, passwrd=None): 
                   
        self._folder = Folder(folderpath, server, username, passwrd)
            
    folder = property(get_folder, set_folder, 
                        doc="The folder for gaussian runs")      

    def count_runs(self):
        """ get number of runs held in analysis """
        return len(self._df.index)        
        
    def _add_molecule(self, molecule, identifiers):
        """add molecule to internal dataframe """

        identifiers['Molecule'] = molecule
        series = pd.DataFrame(identifiers, 
                              index=[self._next_index])
        self._df = self._df.copy().append(series)
        self._next_index += 1
        
        return True
    
    def add_run(self, identifiers={}, 
                      init_fname=None, opt_fname=None, 
                      freq_fname=None, nbo_fname=None,
                      alignto=[], atom_groups={},
                      add_if_error=False, folder_obj=None):
        """add single Gaussian run input/outputs """             
        if not folder_obj:
            folder_obj = self._folder
            
        molecule = Molecule(init_fname=init_fname, 
                            opt_fname=opt_fname, 
                            freq_fname=freq_fname, 
                            nbo_fname=nbo_fname,
                            folder_obj=folder_obj,
                            alignto=alignto, atom_groups=atom_groups,
                            fail_silently=True)
        
        num_files = filter(lambda x:x, [init_fname, opt_fname, 
                                        freq_fname, nbo_fname])
        read_errors = molecule.get_init_read_errors()
        if len(read_errors) != num_files and (not read_errors or add_if_error):                
                self._add_molecule(molecule, identifiers)
                
        return molecule.get_init_read_errors() 
             
    def _get_molecules(self, mol_inputs, folder_obj, identifiers, ipython_print=False):
        """ get molecules """

        if folder_obj.islocal() and not platform.system() == 'Windows':    
            pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
            mapping = pool.imap
        else:
            mapping = imap
        
        with folder_obj:
            molecules=[]
            all_read_errors = []
            for molecule in mapping(unpack_and_make_molecule, mol_inputs):
                molecules.append(molecule)

                read_errors = []
                for typ, fname, msg in molecule.get_init_read_errors():
                    idents = identifiers[len(molecules)-1].copy()
                    idents.pop('Molecule', '_')
                    idents['Type'] = typ
                    idents['File'] = fname
                    idents['Error_Message'] = msg
                    read_errors.append(idents)
                all_read_errors.append(read_errors)

                if ipython_print: 
                    print 'Reading data {0} of {1}'.format(len(molecules), 
                                                           len(mol_inputs))
                    try:
                        clear_output(wait=True)    
                    except:
                        pass

        if folder_obj.islocal() and not platform.system() == 'Windows':
            pool.close()
            pool.join()
                    
        return molecules, all_read_errors                    

    def add_runs(self, headers=[], values=[], 
                 init_pattern=None, opt_pattern=None, 
                 freq_pattern=None, nbo_pattern=None,
                 add_if_error=False,
                 alignto=[], atom_groups={},
                 ipython_print=False, folder_obj=None):
        """add multiple Gaussian run inputs/outputs """ 
        # set folder oject
        if not folder_obj:
            folder_obj = self._folder

        #get variables for each run
        mol_inputs = []
        identifiers = []
        for idents in product(*values):
            mol_input = {}
            identifiers.append(dict(zip(headers, idents)))
            mol_input['init_fname'] = init_pattern.format(*idents) if init_pattern else None
            if type(opt_pattern) is str:
                mol_input['opt_fname'] = opt_pattern.format(*idents) if opt_pattern else None
            elif type(opt_pattern) is list or type(opt_pattern) is tuple:
                mol_input['opt_fname'] = [o.format(*idents) for o in opt_pattern]
            else:
                mol_input['opt_fname'] = None
            mol_input['freq_fname'] = freq_pattern.format(*idents) if freq_pattern else None
            mol_input['nbo_fname'] = nbo_pattern.format(*idents) if nbo_pattern else None
            
            mol_input['folder_obj'] = folder_obj
            mol_input['alignto'] = alignto 
            mol_input['atom_groups'] = atom_groups
            mol_input['fail_silently'] = True
                            
            mol_inputs.append(mol_input)
            
        #create the molecules
        molecules, read_errors = self._get_molecules(mol_inputs, folder_obj, 
                                                     identifiers, ipython_print)

        #add the molecules to the internal table  
        for molecule, idents, inputs, read_error in zip(molecules, identifiers, 
                                                         mol_inputs, read_errors):
            num_files = filter(lambda x:x, [inputs['init_fname'], inputs['opt_fname'], 
                                            inputs['freq_fname'], inputs['nbo_fname']])
            if read_error != num_files and (not read_error or add_if_error):
                self._add_molecule(molecule, idents)
        
        #collate read errors into a dataframe to return  
        read_errors = filter(len, read_errors)                         
        err_df = pd.DataFrame([item for sublist in read_errors for item in sublist])
        if read_errors:
            cols = err_df.columns.tolist()
            #rearrange columns headers
            cols.remove('Type'); cols.append('Type')
            cols.remove('File'); cols.append('File')
            cols.remove('Error_Message'); cols.append('Error_Message')
            err_df = err_df[cols]
        
        return err_df
                
    def get_table(self, rows=[], columns=[],  filters={},
                  precision=4, head=False, mol=False, 
                  row_index=[], column_index=[], 
                  as_image=False, na_rep='-', font_size=None,
                  width=None, height=None, unconfined=False):
        """return pandas table of requested data in requested format

        Parameters
        -----------
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
        
        Returns
        -------
        df : pandas.DataFrame
            a table of data            
            
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
            cols = columns[:]
            if type(row_index) is list:
                cols += row_index
            else:
                cols.append(row_index)
            if mol:
                cols.append('Molecule')
            unique_cols = []
            [unique_cols.append(x) for x in cols if x not in unique_cols]
            df = df.ix[:,unique_cols]            
            
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
        """remove one or more rows of molecules

        Parameters
        ----------
        rows : int or list of ints:
            the rows to remove
        
        """
        self._df.drop(rows, inplace=True)

        return self.get_table()
        
    def remove_columns(self, columns):

        self._df.drop(columns, axis=1, inplace=True) 

        return self.get_table()

    _basic_properties={'nbasis':'get_basis_funcs',
                        'basis':'get_basis_descript',
                       'optimised':'is_optimised',
                       'opt_error': 'get_run_error',
                       'conformer': 'is_conformer'}
                       
    def get_basic_property(self, prop, *args, **kwargs):
        """returns a series of a basic run property or nan if it is not available

        Parameters
        ----------
        prop : str
            can be 'basis', 'nbasis', 'optimised', 'opt_error' or 'conformer'        
        """
        if prop not in self._basic_properties.keys():
            raise ValueError('{0} not a molecule property'.format(prop))
        
        def get_prop(m):
            method = getattr(m, self._basic_properties[prop])
            try: 
                out = method(*args, **kwargs)
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
        
    def remove_non_conformers(self, cutoff=0.):
        """removes runs with negative frequencies """
        non_conformers = self._df[self.get_basic_property('conformer', cutoff=cutoff)!=True].copy()
        self._df = self._df[self.get_basic_property('conformer', cutoff=cutoff)==True]
        return non_conformers

    def add_mol_property(self, name, method, *args, **kwargs):
        """compute molecule property for all rows and create a data column 
        
        Parameters
        ----------
        name : str
            what to name the data column
        method : str
            what molecule method to call 
        *args : various
            arguments to pass to the molecule method
        **kwargs : various
            keyword arguments to pass to the molecule method

        """
                
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
                                     args=[], kwargs={},
                                     relative_to_rows=[]):
        """compute molecule property for a subset of rows and create/add-to data column 

        Parameters
        ----------
        name : str or list of strings
            name for output column (multiple if method outputs more than one value)
        method : str
            what molecule method to call 
        rows : list
            what molecule rows to calculate the property for
        filters : dict
            filter for selecting molecules to calculate the property for
        args : list
            the arguments to pass to the molecule method
        kwargs : dict
            the keyword arguments to pass to the molecule method
        relative_to_rows: list of ints
            compute values relative to the summated value(s) of molecule at the 
            rows listed
        
        """
        df = self.get_table(rows=rows, filters=filters, mol=True)

        if relative_to_rows:
            rel_df = self.get_table(rows=relative_to_rows, mol=True)

        if type(name) is tuple or type(name) is list:
            
            for idx, n in enumerate(name):
                func = lambda m: getattr(m, method)(*args, **kwargs)[idx]
                vals = df.Molecule.map(func)
                
                if relative_to_rows:
                    rel_val = rel_df.Molecule.map(func).sum()
                    vals = vals - rel_val
                    
                if n in self._df.columns:
                    self._df[n] = vals.combine_first(self._df[n])
                else:                
                    self._df[n] = vals
            
    
        else:
            func = lambda m: getattr(m, method)(*args, **kwargs)
            vals = df.Molecule.map(func)
                
            if relative_to_rows:
                rel_val = rel_df.Molecule.map(func).sum()
                vals = vals - rel_val
                    
            if name in self._df.columns:
                self._df[name] = vals.combine_first(self._df[name])
            else:                
                self._df[name] = vals

        return self.get_table()
            
    def get_ids(self, variable_names, variable_lists):
        """return ids of a list of unique computations """
        df = self.get_table()
        df['Index'] = df.index
        df.set_index(variable_names, inplace=True)
        df.sortlevel(inplace=True)
        
        ids = []
        for variable_lst in variable_lists:
            df1 = df.copy()
            try:                
                for v in variable_lst:
                    df1 = df1.loc[v]
            except KeyError:
                raise ValueError(
                        'could not find variable set; {}'.format(variable_lst))
            i = df1.Index
            if hasattr(i, 'values'):
                raise ValueError(
                        'variable set is not unique; {}'.format(variable_lst))
            ids.append(int(i))
        return ids

    
    def get_molecule(self, row):
        """ get molecule object coresponding to particular row """
        return copy.deepcopy(self._df.Molecule.loc[row])
    
    ## TODO will active work?
    def yield_mol_images(self, rows=[], filters={}, mtype='optimised',
                         sort_columns=[], align_to=[], rotations=[[0., 0., 0.]],
                         gbonds=True, represent='ball_stick', 
                         zoom=1., width=300, height=300, axis_length=0,
                         relative=False, minval=-1, maxval=1,
                         highlight=[], active=False, 
                         sopt_min_energy=20., sopt_cutoff_energy=0.,
                         atom_groups=[], alpha=0.5, transparent=False,
                         hbondwidth=5, eunits='kJmol-1', no_hbonds=False, 
                         ipyimg=True):
        """yields molecules

        Parameters
        ----------
        mtype :
            'initial', 'optimised', 'nbo', 'highlight', 'highlight-initial', 'sopt' or 'hbond'
        info_columns : list of str
            columns to use as info in caption
        max_cols : int
            maximum columns in plot
        label_size : int
            subplot label size (pts)
        start_letter : str
            starting (capital) letter for labelling subplots
        save_fname : str
            name of file, if you wish to save the plot to file
        rows : int or list
            index for the row of each molecule to plot (all plotted if empty)
        filters : dict
            {columns:values} to filter by
        sort_columns : list of str
            columns to sort by
        align_to : [int, int, int]
            align geometries to the plane containing these atoms
        rotations : list of [float, float, float]
            for each rotation set [x,y,z] an image will be produced 
        gbonds : bool
            guess bonds between atoms (via distance)
        represent : str
            representation of molecule ('none', 'wire', 'vdw' or 'ball_stick')
        zoom : float
            zoom level of images
        width : int
            width of original images
        height : int
            height of original images (although width takes precedent)
        axis_length : float
            length of x,y,z axes in negative and positive directions
        relative : bool
            coloring of nbo atoms scaled to min/max values in atom set (for nbo mtype)
        minval : float
            coloring of nbo atoms scaled to absolute min (for nbo mtype)
        maxval : float
            coloring of nbo atoms scaled to absolute max (for nbo mtype)
        highlight : list of lists
            atom indxes to highlight (for highlight mtype)
        eunits : str
            the units of energy to return (for sopt/hbond mtype)
        sopt_min_energy : float
            minimum energy to show (for sopt/hbond mtype)
        sopt_cutoff_energy : float
            energy below which bonds will be dashed (for sopt mtype)
        alpha : float
            alpha color value of geometry (for highlight/sopt/hbond mtypes)
        transparent : bool
            whether atoms should be transparent (for highlight/sopt/hbond mtypes)
        hbondwidth : float   
            width of lines depicting interaction (for hbond mtypes)  
        atom_groups : [list or str, list or str]
            restrict interactions to between two lists (or identifiers) of atom indexes (for sopt/hbond mtypes)
        no_hbonds : bool
            whether to ignore H-Bonds in the calculation
        ipyimg : bool
            whether to return an IPython image, PIL image otherwise 
        
        Yields
        -------
        indx : int
            the row index of the molecule
        mol : IPython.display.Image or PIL.Image
            an image of the molecule in the format specified by ipyimg            
        
        """
        df = self.get_table(columns=['Molecule']+sort_columns, rows=rows, 
                       filters=filters, mol=True)
        
        if sort_columns:
            df.sort(sort_columns, inplace=True)
        
        for indx, mol in zip(df.index, df.Molecule):
            if align_to: 
                align_atoms = mol.get_atom_group(align_to)
                mol.set_alignment_atoms(*align_atoms)
            if mtype == 'initial':
                yield indx, mol.show_initial(gbonds=gbonds, represent=represent, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length, ipyimg=ipyimg)
            elif mtype == 'optimised':
                yield indx, mol.show_optimisation(gbonds=gbonds, represent=represent, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length, ipyimg=ipyimg)
            elif mtype == 'nbo':
                yield indx, mol.show_nbo_charges(gbonds=gbonds, represent=represent, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length,
                                       relative=relative, 
                                       minval=minval, maxval=maxval, ipyimg=ipyimg)
            elif mtype == 'highlight':
                yield indx, mol.show_highlight_atoms(highlight, 
                                       alpha=alpha, optimised=True,
                                       transparent=transparent,
                                       gbonds=gbonds, 
                                       represent=represent, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length, ipyimg=ipyimg)
            elif mtype == 'highlight-initial':
                yield indx, mol.show_highlight_atoms(highlight, 
                                       alpha=alpha, optimised=False,
                                       transparent=transparent,
                                       gbonds=gbonds, 
                                       represent=represent, 
                                       rotations=rotations, zoom=zoom, 
                                       width=width, height=height, 
                                       axis_length=axis_length, ipyimg=ipyimg)
            elif mtype == 'sopt':
                yield indx, mol.show_sopt_bonds(min_energy=sopt_min_energy,
                                    cutoff_energy=sopt_cutoff_energy, no_hbonds=no_hbonds, 
                                    eunits=eunits, atom_groups=atom_groups,
                                    alpha=alpha, transparent=transparent,
                                    gbonds=gbonds, represent=represent, 
                                    rotations=rotations, zoom=zoom, 
                                    width=width, height=height, 
                                    axis_length=axis_length,
                                    relative=relative, 
                                    minval=minval, maxval=maxval, ipyimg=ipyimg)
            elif mtype == 'hbond':
                yield indx, mol.show_hbond_analysis(min_energy=sopt_min_energy,
                                    cutoff_energy=sopt_cutoff_energy, eunits=eunits,
                                    atom_groups=atom_groups, bondwidth=hbondwidth,
                                    alpha=alpha, transparent=transparent,
                                    gbonds=gbonds, represent=represent, 
                                    rotations=rotations, zoom=zoom, 
                                    width=width, height=height, 
                                    axis_length=axis_length,
                                    relative=relative, 
                                    minval=minval, maxval=maxval, ipyimg=ipyimg)
            else:
                raise ValueError(
                'mtype must be initial, optimised, nbo, highlight, highligh-initial, sopt or hbond')                
    
    def _get_letter(self, number):
        """get an uppercase letter according to a number"""
        if number < 26:
            return string.ascii_uppercase[number]
        else:
            first_letter = string.ascii_uppercase[int(number/26)-1]
            second_letter = string.ascii_uppercase[number % 26]
            return first_letter + second_letter
        
    def plot_mol_images(self, mtype='optimised', max_cols=1, padding=(1, 1),
                        sort_columns=[], info_columns=[], info_incl_id=False,
                        label_size=20, letter_prefix='', start_letter='A',
                        rows=[], filters={}, align_to=[], rotations=[[0., 0., 0.]],
                        gbonds=True, represent='ball_stick',
                        zoom=1., width=500, height=500, axis_length=0,
                        relative=False, minval=-1, maxval=1,
                        highlight=[], frame_on=False, eunits='kJmol-1',
                        sopt_min_energy=20., sopt_cutoff_energy=0.,
                        atom_groups=[], alpha=0.5, transparent=False,
                        hbondwidth=5, no_hbonds=False):
        """show molecules in matplotlib table of axes

        Parameters
        ----------
        mtype :
            'initial', 'optimised', 'nbo', 'highlight', 'highlight-initial', 'sopt' or 'hbond'
        max_cols : int
            maximum columns in plot
        padding: tuple
            padding between images (horizontally, vertically)
        sort_columns : list of str
            columns to sort by
        info_columns : list of str
            columns to use as info in caption
        info_incl_id : bool
            include molecule id number in caption
        label_size : int
            subplot label size (pts)
        letter_prefix : str
            prefix for labelling subplots
        start_letter : str
            starting (capital) letter for labelling subplots
        rows : int or list
            index for the row of each molecule to plot (all plotted if empty)
        filters : dict
            {columns:values} to filter by
        align_to : [int, int, int]
            align geometries to the plane containing these atoms
        rotations : list of [float, float, float]
            for each rotation set [x,y,z] an image will be produced 
        gbonds : bool
            guess bonds between atoms (via distance)
        represent : str
            representation of molecule ('none', 'wire', 'vdw' or 'ball_stick')
        zoom : float
            zoom level of images
        width : int
            width of original images
        height : int
            height of original images (although width takes precedent)
        axis_length : float
            length of x,y,z axes in negative and positive directions
        relative : bool
            coloring of nbo atoms scaled to min/max values in atom set (for nbo mtype)
        minval : float
            coloring of nbo atoms scaled to absolute min (for nbo mtype)
        maxval : float
            coloring of nbo atoms scaled to absolute max (for nbo mtype)
        highlight : list of lists
            atom indxes to highlight (for highlight mtype)
        eunits : str
            the units of energy to return (for sopt/hbond mtype)
        sopt_min_energy : float
            minimum energy to show (for sopt/hbond mtype)
        sopt_cutoff_energy : float
            energy below which bonds will be dashed (for sopt mtype)
        alpha : float
            alpha color value of geometry (for sopt/hbond mtypes)
        transparent : bool
            whether atoms should be transparent (for sopt/hbond mtypes)
        hbondwidth : float   
            width of lines depicting interaction (for hbond mtypes)     
        atom_groups : [list or str, list or str]
            restrict interactions to between two lists (or identifiers) of atom indexes (for sopt/hbond mtypes)
        no_hbonds : bool
            whether to ignore H-Bonds in the calculation (for sopt only)
        frame_on : bool
            whether to show frame around each image 

        Returns
        -------
        fig : matplotlib.figure.Figure
            A figure containing subplots for each molecule image
        caption : str
            A caption describing each subplot, given info_columns
        
        """
        letter_offset = string.ascii_uppercase.find(start_letter)
        if letter_offset == -1:
            raise ValueError('start_letter must be an uppercase single letter')
        
        df = self.get_table(rows=rows, columns=info_columns, filters=filters)
        num_mols = len(df)
        
        imgs = self.yield_mol_images(rows=rows, filters=filters, mtype=mtype,
                        align_to=align_to, gbonds=gbonds, represent=represent, 
                        sort_columns=sort_columns, rotations=rotations, zoom=zoom, 
                        width=width, height=height, axis_length=axis_length,
                        relative=relative, minval=minval, maxval=maxval,
                        highlight=highlight, active=False, ipyimg=False,
                        eunits=eunits, sopt_min_energy=sopt_min_energy, 
                        sopt_cutoff_energy=sopt_cutoff_energy,
                        atom_groups=atom_groups, alpha=alpha, 
                        transparent=transparent,
                        hbondwidth=hbondwidth, no_hbonds=no_hbonds)
        
        #num_rows = int(math.ceil(num_mols/float(max_cols)))
        #num_cols = min([max_cols, num_mols])
        num_cols=int(max_cols)
        num_rows=int(math.ceil(num_mols/float(num_cols)))
 
        fig, axes = plt.subplots(num_rows, num_cols, squeeze=False,
                                 gridspec_kw={'width_ratios':[1]*num_cols})
                                 
        for ax in fig.get_axes():
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)
            ax.set_anchor('NW')
            ax.set_frame_on(False)

        mol_num = 0
        caption = []                       
        for indx, img in imgs:

            ax = axes[int(math.ceil((mol_num+1)/float(max_cols)))-1,
                      mol_num % max_cols]
            ax.imshow(img)#, aspect='equal')
            ax.set_frame_on(frame_on)

            if label_size:
                ax.text(0,0.8,letter_prefix+self._get_letter(mol_num+letter_offset), 
                        size=label_size, weight="bold")
            
            info = ', '.join(df[info_columns].loc[indx].fillna('-').astype(str))
            if info_incl_id:
                info = str(indx) + ', ' + info
                                
            caption.append(
                '(' + letter_prefix+self._get_letter(mol_num+letter_offset) + ') ' + info)
            
            mol_num += 1            

        #resize extra axes to be same as last img
        while mol_num < num_rows*num_cols:
            ax = axes[int(math.ceil((mol_num+1)/float(max_cols)))-1,
                      mol_num % max_cols]
            ax.imshow(img)
            ax.clear()
            mol_num += 1               

        fig.tight_layout(w_pad=padding[0], h_pad=padding[1])
                
        caption = ', '.join(caption)
        #insert newline character every 80 charaters
        #caption = re.sub("(.{80})", "\\1\n", caption, 0, re.DOTALL)
        
        return fig, caption
        
                
    def plot_mol_graphs(self, gtype='energy', share_plot=False, max_cols=1, 
                    padding=(1,1), tick_rotation=0,
                    rows=[], filters={}, sort_columns=[], info_columns=[], 
                    info_incl_id=False, letter_prefix='', start_letter='A',
                    grid=True, sharex=True, sharey=True, legend_size=10,
                    color_scheme='jet', eunits='eV',
                    per_energy=1., lbound=None, ubound=None,
                    color_homo='g', color_lumo='r', 
                    homo_lumo_lines=True,homo_lumo_values=True,band_gap_value=True):
        """get a set of data plots for each molecule
        
        Parameters
        ----------
        gtype : str
            the type of plot, 
            energy = optimisation energies, 
            freq = frequency analsis,
            dos = Densty of States,
        share_plot : bool
            whether to plot all data on the same or separate axes
        max_cols : int
            maximum columns on plots (share_plot=False only)
        padding: tuple
            padding between images (horizontally, vertically)
        tick_rotation : int
            rotation of x-axis labels
        rows : int or list
            index for the row of each molecule to plot (all plotted if empty)
        filters : dict
            {columns:values} to filter by
        sort_columns : list of str
            columns to sort by
        info_columns : list of str
            columns to use as info in caption
        info_incl_id : bool
            include molecule id number in labels
        letter_prefix : str
            prefix for labelling subplots (share_plot=False only)
        start_letter : str
            starting (capital) letter for labelling subplots (share_plot=False only)
        grid : bool
            whether to include a grid in the axes
        sharex : bool
            whether to align x-axes (share_plot=False only)
        sharey : bool
            whether to align y-axes (share_plot=False only)
        legend_size : int
            the font size (in pts) for the legend
        color_scheme : str
            the scheme to use for each molecule (share_plot=True only)
            according to http://matplotlib.org/examples/color/colormaps_reference.html 
        eunits : str
            the units of energy to use
        per_energy : float
            energy interval to group states by (DoS only)
        lbound : float
            lower bound energy (DoS only)
        ubound: float
            upper bound energy (DoS only)
        color_homo : matplotlib.colors
            color of homo in matplotlib format
        color_lumo : matplotlib.colors
            color of lumo in matplotlib.colors
        homo_lumo_lines : bool
            draw lines at HOMO and LUMO energies
        homo_lumo_values : bool
            annotate HOMO and LUMO lines with exact energy values
        band_gap_value : bool
            annotate inbetween HOMO and LUMO lines with band gap value
        
        Returns
        -------
        data : matplotlib.figure.Figure
            plotted frequency data
        caption : str
            A caption describing each subplot, given info_columns
            
        """
        df = self.get_table(columns=list(set(info_columns+sort_columns)), 
                            rows=rows, 
                            filters=filters, mol=True)
        num_plots = df.index.shape[0]
        
        if sort_columns:
            df.sort(sort_columns, inplace=True)        

        if gtype == 'energy':
            mol_func = 'plot_opt_energy'
            x_label = 'Optimisation Step'
            y_label = 'Energy ({0})'.format(eunits)
            all_plot_kwargs = {'units':eunits}
            per_plot_kwargs = {'linecolor':getattr(cm,color_scheme)(
                                            np.linspace(0.1, 0.9, num_plots))}
        elif gtype == 'freq':
            mol_func = 'plot_freq_analysis'
            x_label = 'Frequency ($cm^{-1}$)'
            y_label = 'IR Intensity ($km/mol$)' 
            all_plot_kwargs = {}
            per_plot_kwargs = {'color':getattr(cm,color_scheme)(
                                                np.linspace(0, 1, num_plots)),
                               'alpha':np.linspace(1, 0.5, num_plots),
                                'marker_size':np.linspace(25, 15, num_plots)}
        elif gtype == 'dos':
            if share_plot:
                raise ValueError('share_plots not available for Density of States')
            mol_func = 'plot_dos'
            x_label = 'Density of States (per {0} {1})'.format(per_energy, eunits)
            y_label = 'Energy ({})'.format(eunits)
            all_plot_kwargs = {'eunits':eunits, 'per_energy':per_energy, 
                               'lbound':lbound, 'ubound':ubound,
                 'color_homo':color_homo, 'color_lumo':color_lumo, 
                 'homo_lumo_lines':homo_lumo_lines, 'homo_lumo_values':homo_lumo_values,
                 'band_gap_value':band_gap_value, 'legend_size':legend_size}
        else:
            raise ValueError('gtype; {0}, not available'.format(gtype))

        ax_num = 0
        caption = []  

        if share_plot:
            fig, ax = plt.subplots()
            
            legend = []
            for indx, row in df.iterrows():
                plot_kwargs = all_plot_kwargs.copy()
                for k, v in per_plot_kwargs.iteritems():
                    plot_kwargs[k] = v[ax_num]
                getattr(row.Molecule, mol_func)(ax=ax, **plot_kwargs)

                label = ', '.join(row[info_columns].fillna('-').astype(str))
                if info_incl_id:
                    label = str(indx) + ', ' + label
                legend.append(label)
                
                ax_num += 1

            ax.grid(grid)
            for tick in ax.get_xticklabels():
                tick.set_rotation(tick_rotation)
            ax.legend(legend, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                      prop={'size':legend_size})
        else:
        
            num_rows = int(math.ceil(num_plots/float(max_cols)))
            num_cols = min([max_cols, num_plots])
            fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, squeeze=False,
                                     sharex=sharex, sharey=sharey)

            letter_offset = string.ascii_uppercase.find(start_letter)
            if letter_offset == -1:
                raise ValueError('start_letter must be an uppercase single letter')

            for indx, row in df.iterrows():
                i = int(math.ceil((ax_num+1)/float(max_cols)))-1
                j = ax_num % max_cols               
                
                getattr(row.Molecule, mol_func)(ax=axes[i,j], **all_plot_kwargs) 
                axes[i,j].grid(grid)
                for tick in axes[i,j].get_xticklabels():
                    tick.set_rotation(tick_rotation)

                info = ', '.join(row[info_columns].fillna('-').astype(str))
                if info_incl_id:
                    info = str(indx) + ', ' + info
                letter = self._get_letter(ax_num+letter_offset)
                axes[i,j].set_title(letter_prefix+letter, fontweight="bold")
                    
                caption.append('(' + letter_prefix+letter + ') ' + info)
                
                ax_num += 1
            
            #hide extraneous axes
            for extra_ax in range(ax_num, num_rows*num_cols):
                i = int(math.ceil((extra_ax+1)/float(max_cols)))-1
                j = extra_ax % max_cols  
                axes[i,j].axis('off')


            ax = fig.add_subplot(111)    # The big subplot
            ax.tick_params(top='off', bottom='off', left='off', right='off',
                           labelbottom='on', labelleft='on', pad=25)
            
            ax.set_xticklabels([])
            ax.set_yticklabels([])           
            ax.set_frame_on(False)
        
        ax.set_xlabel(x_label)  
        ax.set_ylabel(y_label)
        
        fig.tight_layout(w_pad=padding[0], h_pad=padding[1])
        
        caption = ', '.join(caption)
        
        return fig, caption
        
        
    def plot_radviz_comparison(self, category_column, 
                               columns=[], rows=[], filters={}, point_size=30,
                                **kwargs):
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
        df = df.sort(category_column)

        f, ax = plt.subplots()
        ax = radviz(df, category_column, ax=ax, s=point_size, **kwargs)
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
                                                                        
                                  
                                  
                                  
                                  
                                  
                                  
