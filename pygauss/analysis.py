# -*- coding: utf-8 -*-
from itertools import product
import copy
import math
import string
import re

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
        if not len(molecule.get_init_read_errors()) == num_files:
            if not molecule.get_init_read_errors() or add_if_error:
                    
                identifiers['Molecule'] = molecule
                series = pd.DataFrame(identifiers, 
                                      index=[self._next_index])
                self._df = self._df.copy().append(series)#, ignore_index=True)
                self._next_index += 1
                
        return molecule.get_init_read_errors() 
             
    def add_runs(self, headers=[], values=[], 
                 init_pattern=None, opt_pattern=None, 
                 freq_pattern=None, nbo_pattern=None,
                 add_if_error=False,
                 alignto=[], atom_groups={},
                 ipython_print=False):
        """add multiple Gaussian run inputs/outputs """             
        with self._folder as folder:
            
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
                            alignto=alignto, atom_groups=atom_groups,
                            add_if_error=add_if_error, folder_obj=folder)
                
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
            
    def get_molecule(self, row):
        """ get molecule object coresponding to particular row """
        return copy.deepcopy(self._df.Molecule.loc[row])
    
    ## TODO will active work?
    def yield_mol_images(self, rows=[], filters={}, mtype='optimised',
                         align_to=[], rotations=[[0., 0., 0.]],
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
            'initial', 'optimised', 'nbo', 'highlight', 'sopt' or 'hbond'
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
        df = self.get_table(columns=['Molecule'], rows=rows, 
                       filters=filters, mol=True)
        
        for indx, mol in zip(df.index, df.Molecule):
            if align_to: mol.set_alignment_atoms(*align_to)
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
                                       alpha=alpha,
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
                'mtype must be initial, optimised, nbo, highlight, sopt or hbond')                
    
    def plot_mol_images(self, mtype='optimised', info_columns=[], info_incl_id=False,
                        max_cols=1, label_size=20, start_letter='A', save_fname=None,
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
            'initial', 'optimised', 'nbo', 'highlight', 'sopt' or 'hbond'
        info_columns : list of str
            columns to use as info in caption
        info_incl_id : bool
            include molecule id number in caption
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
                        rotations=rotations, zoom=zoom, 
                        width=width, height=height, axis_length=axis_length,
                        relative=relative, minval=minval, maxval=maxval,
                        highlight=highlight, active=False, ipyimg=False,
                        eunits=eunits, sopt_min_energy=sopt_min_energy, 
                        sopt_cutoff_energy=sopt_cutoff_energy,
                        atom_groups=atom_groups, alpha=alpha, 
                        transparent=transparent,
                        hbondwidth=hbondwidth, no_hbonds=no_hbonds)
        
        num_rows = int(math.ceil(num_mols/float(max_cols)))
        num_cols = min([max_cols, num_mols])

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
                ax.text(0,0.8,string.ascii_uppercase[mol_num+letter_offset], 
                        size=label_size, weight="bold")
            
            info=[]
            if info_incl_id:
                info.append(str(indx))
            for col in info_columns:
                try:
                    isnan = pd.np.isnan(df[col].loc[indx])
                except TypeError:
                    isnan = False
                if isnan:
                    value = '-'
                else:
                    value = df[col].loc[indx]
                    
                info.append(str(value))
                                
            caption.append(
                '(' + string.ascii_uppercase[mol_num+letter_offset] + ') ' + ', '.join(info))
            
            mol_num += 1                            

        fig.tight_layout(h_pad=2.0)
        
        if save_fname:
            self._folder.save_mplfig(fig, save_fname)
        
        caption = 'Figure: ' + ', '.join(caption)
        #insert newline character every 80 charaters
        caption = re.sub("(.{80})", "\\1\n", caption, 0, re.DOTALL)
        
        return fig, caption
        
    def get_freq_analysis(self, info_columns=[], rows=[], filters={}):
        """return frequency analysis 

        Parameters
        ----------
        info_columns : list of str
            columns to use as info in caption
        rows : int or list
            index for the row of each molecule to plot (all plotted if empty)
        filters : dict
            {columns:values} to filter by
        
        Returns
        -------
        data : pd.DataFrame
            frequency data
        """
        
        df = self.get_table(columns=info_columns, rows=rows, 
                            filters=filters, mol=True)
        
        main = pd.DataFrame()
        for indx, row in df.iterrows():
            df = row.Molecule.get_freq_analysis()
            
            df['row']=indx
            for col in info_columns:
                df[col] = row[col]
            main = main.append(df)
        
        return main
        
    def plot_freq_analysis(self, info_columns=[], rows=[], filters={}, 
                           share_plot=True, include_row=False):
        """plot frequency analysis 

        Parameters
        ----------
        info_columns : list of str
            columns to use as info in caption
        rows : int or list
            index for the row of each molecule to plot (all plotted if empty)
        filters : dict
            {columns:values} to filter by
        share_plot : bool
            whether to share a single plot or have multiple ones
        include_row : bool
            include row number in legend labels
        
        Returns
        -------
        data : matplotlib.figure.Figure
            plotted frequency data
        """

        df = self.get_freq_analysis(info_columns=info_columns, rows=rows, 
                                    filters=filters)
        
        if share_plot:
            fig, ax = plt.subplots()
        
            colors = cm.rainbow(np.linspace(0, 1, df.row.unique().shape[0]))
            alphas = np.linspace(1, 0.5, df.row.unique().shape[0])
            m_sizes = np.linspace(25, 15, df.row.unique().shape[0])
            
            for data, c, a, s in zip(df.groupby('row'), colors, alphas, m_sizes):
                row, group = data
                label = ', '.join([str(group[col].iloc[0]) for col in info_columns])
                if include_row:
                    label += ' ({0})'.format(row)
                ax.bar(group['Frequency ($cm^{-1}$)'], group['IR Intensity ($km/mol$)'], 
                         align='center', width=30, linewidth=0,  alpha=a,color=c, label=label)
                ax.scatter(group['Frequency ($cm^{-1}$)'], group['IR Intensity ($km/mol$)'], 
                                 marker='o', alpha=a, color=c, s=s)

            ax.grid()
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            ax.set_ybound(-10)
            ax.set_xlabel('Frequency ($cm^{-1}$)')
            ax.set_ylabel('IR Intensity ($km/mol$)')       
        else:
            fig, axes = plt.subplots(df.row.unique().shape[0], sharex=True, sharey=True)            
            for data, ax in zip(df.groupby('row'), axes):
                row, group = data
                label = ', '.join([str(group[col].iloc[0]) for col in info_columns])
                if include_row:
                    label += ' ({0})'.format(row)
                ax.bar(group['Frequency ($cm^{-1}$)'], group['IR Intensity ($km/mol$)'], 
                         align='center', width=30, linewidth=0)
                ax.scatter(group['Frequency ($cm^{-1}$)'], group['IR Intensity ($km/mol$)'], 
                                 marker='o')
                ax.set_title(label)
                ax.grid()
                ax.set_ybound(-10)

            ax = fig.add_subplot(111)    # The big subplot
            ax.tick_params(top='off', bottom='off', left='off', right='off',
                           labelbottom='on', labelleft='on', pad=25)
            
            ax.set_xticklabels([])
            ax.set_yticklabels([])           
            ax.set_frame_on(False)
            ax.set_ylabel('IR Intensity ($km/mol$)')  
            ax.set_xlabel('Frequency ($cm^{-1}$)')
            
        fig.tight_layout()
        return fig
        
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
                                                                        
                                  
                                  
                                  
                                  
                                  
                                  