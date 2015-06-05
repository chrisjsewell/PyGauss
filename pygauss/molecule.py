# -*- coding: utf-8 -*-
"""
Created on Fri May 01 21:24:31 2015

@author: chris
"""
from io import BytesIO
import PIL
from PIL import Image, ImageChops
import copy 
import warnings

from math import degrees, atan2, sqrt, acos

import numpy as np
from scipy.signal import argrelextrema

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from mpl_toolkits.mplot3d import Axes3D

import pandas as pd

from IPython.display import Image as ipy_Image

from .chemlab_patch.io.handlers._cclib  import _create_cclib_handler

from chemlab.graphics.qtviewer import QtViewer
#have to add method to instances of chemlab.graphics.camera.Camera
#from chemlab.graphics.transformations import rotation_matrix
from .transformations import rotation_matrix
def orbit_z(self, angle):
    # Subtract pivot point
    self.position -= self.pivot        
    # Rotate
    rot = rotation_matrix(-angle, self.c)[:3,:3]
    self.position = np.dot(rot, self.position)        
    # Add again the pivot point
    self.position += self.pivot
    
    self.a = np.dot(rot, self.a)
    self.b = np.dot(rot, self.b)
    self.c = np.dot(rot, self.c)     
from chemlab.graphics.camera import Camera
Camera.orbit_z = orbit_z

from .chemlab_patch.graphics.renderers.atom import AtomRenderer
from .chemlab_patch.graphics.renderers.ballandstick import BallAndStickRenderer
from .chemlab_patch.graphics.renderers.line import LineRenderer
from .chemlab_patch.graphics.renderers.triangles import TriangleRenderer
from chemlab.graphics.renderers.wireframe import WireframeRenderer
#from chemlab.graphics.postprocessing import SSAOEffect # Screen Space Ambient Occlusion
from chemlab.utils import cartesian
from cclib.parser.utils import convertor

from chemlab.graphics.colors import get as str_to_colour
from chemlab.qc import molecular_orbital

#improvement to function
#TODO not making this a depedency until it works
from chemlab.qc.pgbf import pgbf
try:
    import numexpr as ne
    def __call__(self,x,y,z):
        "Compute the amplitude of the PGBF at point x,y,z"
        I,J,K = self.powers
        dx,dy,dz = x-self.origin[0],y-self.origin[1],z-self.origin[2]
        n = self.norm
        e = self.exponent
        return ne.evaluate(
                'n*(dx**I)*(dy**J)*(dz**K) * exp(-e*(dx*dx + dy*dy + dz*dz))')
    pgbf.__call__ = __call__
except ImportError:
    pass

#instead of chemview MolecularViewer to add defined colouring
#also ignore; 'FutureWarning: IPython widgets are experimental and may change in the future.'
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from .chemview_patch.viewer import MolecularViewer

from .utils import circumcenter
from .file_io import Folder
from .isosurface import get_isosurface

class Molecule(object):
    
    def __init__(self, folderpath='', 
                 init_fname=False, opt_fname=False, 
                 freq_fname=False, nbo_fname=False, 
                 pes_fname=False, 
                 fail_silently=False,
                 atom_groups={}, alignto=[],
                 server=None, username=None, passwrd=None,
                 folder_obj=None):
        """a class to analyse gaussian input/output of a single molecular geometry 
        
        Parameters
        ----------
        folderpath : str
            the folder path
        init_fname : str
            the intial geometry (.com) file
        opt_fname : str or list of str
            the optimisation log file
        freq_fname : str
            the frequency analysis log file
        nbo_fname : str
            the population analysis logfile
        pes_fname : str
            the potential energy scan logfile
        fail_silently : bool
            whether to raise an error if a file read fails (if True can use get_init_read_errors to see errors)
        atom_groups: {str:[int, ...]}
            groups of atoms that can be selected as a subset
        alignto: [int, int, int]
            the atom numbers to align the geometry to
        
        any of the file names can have wildcards (e.g. 'filename*.log) in them, 
        as long as this resolves to a single path in the directory 
        
        NB: nbo population analysis must be run with the GFInput flag to ensure 
        data is output to the log file 
        """
        if folder_obj:
            self._folder = folder_obj
        else:
            self._folder = Folder(folderpath, 
                                  server, username, passwrd)            
        
        self._init_data = None
        self._prev_opt_data = []
        self._opt_data = None
        self._freq_data = None
        self._nbo_data = None
        self._pes_data = []
        self._alignment_atom_indxs = ()
        self._t_matrix = None
        if alignto: 
            self.set_alignment_atoms(*alignto)
        self._atom_groups = atom_groups
                
        parts=[[init_fname, self.add_initialgeom],
               [opt_fname,  self.add_optimisation],
               [freq_fname, self.add_frequency],
               [nbo_fname, self.add_nbo_analysis],
               [pes_fname, self.add_pes_analysis]]
        self._init_read_errors = []
        for fname, method in parts:
            if fname:
                if fail_silently:
                    try:
                        method(fname)
                    except Exception, e:
                        self._init_read_errors.append([fname, str(e)])   
                else:
                    method(fname)
    
    def get_folder(self):
        """ return the Folder instance """
        return self._folder
        
    def get_init_read_errors(self):
        """ get read errors, recorded if fail_silently was set to True on initialise """
        return self._init_read_errors[:]
                     
    def __repr__(self):
        return '<PyGauss Molecule>'
            
    def __deepcopy__(self, memo):
        if not self._folder.islocal():
            warnings.warn('Cannot deepcopy a molecule created via non-local IO')
            return copy.copy(self)
        else:
            cls = self.__class__
            result = cls.__new__(cls)
            memo[id(self)] = result
            for k, v in self.__dict__.items():
                setattr(result, k, copy.deepcopy(v, memo))
            return result
        
    def _get_data(self, file_name, ftype='gaussian'):

        gaussian_handler = _create_cclib_handler(ftype)
        
        if not self._folder.active():
            with self._folder as folder:
                with folder.read_file(file_name) as fd:            
                    data = gaussian_handler(fd)
        else:
            with self._folder.read_file(file_name) as fd:            
                data = gaussian_handler(fd)
        
        return data       
    
    def add_initialgeom(self, file_name):
        
        self._init_data = self._get_data(file_name, ftype='gausscom')

    def add_optimisation(self, file_name):
                
        if type(file_name) is list or type(file_name) is tuple:
            self._opt_data = self._get_data(file_name[-1])
            self._prev_opt_data = []
            for f in file_name[:-1]:
                self._prev_opt_data.append(self._get_data(f))
        else:
            self._opt_data = self._get_data(file_name)

    def add_frequency(self, file_name):
        
        self._freq_data = self._get_data(file_name)
    
    def add_nbo_analysis(self, file_name):
        
        self._nbo_data = self._get_data(file_name)

    def add_pes_analysis(self, file_names):
        
        if type(file_names) is str:
            file_names = [file_names]
        self._pes_data = [self._get_data(fname) for fname in file_names]
    
    def _read_data(self, ftype, dtype):
        """ read data """
        if not getattr(self, ftype):
            raise ValueError(
                '{0} has not been set for this molecule'.format(ftype))
        return getattr(self, ftype).read(dtype)

    def get_basis_descript(self):
        
        return self._read_data('_opt_data', 'basis_descript')

    def get_basis_funcs(self):
        
        return self._read_data('_opt_data', 'nbasis')
        
    def get_run_error(self, rtype='opt'):
        """True if there were errors in the computation, else False """
        return getattr(self, '_{0}_data'.format(rtype)).read('run_error')
        
    def is_optimised(self):
        """ was the geometry optimised """
        return self._read_data('_opt_data', 'optdone')

    def get_optimisation_E(self, units='eV', final=True):
        """ return the SCF optimisation energy 
        
        Parameters
        ----------
        units : str
            the unit type of the energy
        final : bool
            return only the final optimised energy if True, else for all steps            
        
        Returns
        -------
        out : float or list of floats
            dependant on final
        """
        
        if not self._opt_data:
            return np.nan
        
        energies = self._read_data('_opt_data', 'scfenergies')
        
        if energies.shape[0] == 0:
            return np.nan if final else energies
        
        if not units == 'eV':
            energies = convertor(energies, 'eV', units)
        
        return energies[-1] if final else energies  
            
    def plot_optimisation_E(self, units='eV'):
        """ plot SCF optimisation energy """
        energies = self._read_data('_opt_data', 'scfenergies')
        for data in reversed(self._prev_opt_data):
            energies = np.concatenate([data.read('scfenergies'), energies])
            
        if not units == 'eV':
            energies = convertor(energies, 'eV', units)
        
        f, ax = plt.subplots()
        ax.plot(energies)
        ax.set_ylabel('Energy ({0})'.format(units))
        ax.set_xlabel('Optimisation Step')
        ax.grid(True)
        
        return ax

    def is_conformer(self, cutoff=0.):
        """False if any frequencies in the frequency analysis were negative"""
        imgaginary_freqs = self._read_data('_freq_data', 'vibfreqs') < cutoff 
        return not imgaginary_freqs.any()
        
    def get_freq_analysis(self):
        """return frequency analysis
        
        Returns
        -------
        data : pd.DataFrame
            frequency data
        """
        frequencies = self._read_data('_freq_data', 'vibfreqs')
        irs = self._read_data('_freq_data', 'vibirs')
       
        return pd.DataFrame(zip(frequencies, irs), 
                             columns=['Frequency ($cm^{-1}$)', 
                             'IR Intensity ($km/mol$)'])

    def plot_freq_analysis(self):
        """plot frequency analysis 

        Returns
        -------
        data : matplotlib.axes._subplots.AxesSubplot
            plotted frequency data
        
        """

        df = self.get_freq_analysis()
        
        fig, ax = plt.subplots()
                
        ax.bar(df['Frequency ($cm^{-1}$)'], df['IR Intensity ($km/mol$)'], 
                 align='center', width=30, linewidth=0)
        ax.scatter(df['Frequency ($cm^{-1}$)'], df['IR Intensity ($km/mol$)'] , 
                      marker='o',alpha=0.7)
        ax.grid()
        ax.set_ybound(-10)
        ax.set_xlabel('Frequency ($cm^{-1}$)')
        ax.set_ylabel('IR Intensity ($km/mol$)')    
        
        return ax

    def set_alignment_atoms(self, idx1, idx2, idx3):
        
        assert type(idx1) is int and type(idx2) is int and type(idx3) is int
        
        self._alignment_atom_indxs = (idx1, idx2, idx3)    
        
    def remove_alignment_atoms(self):
        
        self._alignment_atom_indxs = ()   

    def _midpoint_coordinates(self, coord_list):

        return np.mean(np.array(coord_list), axis=0)                                 

    def _midpoint_atoms(self, molecule, atom_ids):

        return np.mean(molecule.r_array[atom_ids], axis=0)                                 

    def _create_transform_matrix(self, c1, c2, c3):
        """
        A function to take three coordinates and creates a transformation matrix
        that aligns their plane with the standard axes 
        
        there centre point will be at (0, 0, 0)
        c1 will be aligned to the x-axis
        the normal to the plane will be aligned to the z-axis
        """
        # find midpoint of coords
        c0 = circumcenter([c1, c2, c3])        
        #c0 = self._midpoint_coordinates([c1, c2, c3])
        
        #translate c0 to the origin [0,0,0] and pick two vectors
        v1=c1-c0; v2=c2-c0; v3=c3-c0
    
        #now find the orthonormal basis set   
        # a plane is a*x+b*y+c*z+d=0 where[a,b,c] is the normal and d is 0  
        # (since the origin now intercepts the plane). Thus, we calculate;
        normal = np.cross(v2,v3)
        #a, b, c = normal
        vf3 = normal/np.linalg.norm(normal)
    
        vf1 =  v1/np.linalg.norm(v1)        
         
        vf2 = np.cross(vf3, vf1)
        vf2 = vf2/np.linalg.norm(vf2)
    
        #create the translation matrix that moves the new basis to the origin
        ident=np.vstack((np.identity(3), np.zeros(3)))
        translate_matrix = np.hstack((ident, np.array(np.append(-c0, 1))[np.newaxis].T))
        #create the rotation matrix that rotates the new basis onto the standard basis
        rotation_matrix = np.hstack((np.array([vf1, vf2, vf3, np.zeros(3)]), 
                                     np.array(np.append([0, 0, 0], 1))[np.newaxis].T))
        # translate before rotating
        transform_matrix = np.dot(rotation_matrix, translate_matrix)
        
        return transform_matrix
    
    def _apply_transfom_matrix(self, transform_matrix, coords):
        
        if transform_matrix is None:
            return coords
            
        t_coords = [np.dot(transform_matrix, 
                    np.array(np.append(coord, 1))[np.newaxis].T)[:-1].flatten()
                    for coord in coords]
        return np.array(t_coords)
    
    def _create_molecule(self, optimised=True, opt_step=False, scan_step=False, 
                         gbonds=True, data=None, alignment_atoms=None):
        if not optimised:
            molecule = self._read_data('_init_data', 'molecule')            
        else:
            indata = data if data else self._opt_data
            if not type(opt_step) is bool:
                molecule = indata.read('molecule', step=opt_step) 
            elif not type(scan_step) is bool:
                molecule = indata.read('molecule', scan=scan_step)
            else:
                molecule = indata.read('molecule') 
            
        if gbonds: molecule.guess_bonds()
        
        t_matrix = None
        if alignment_atoms:
            a, b, c = alignment_atoms
            a0, b0, c0 = molecule.r_array[[a-1,b-1,c-1]]
            t_matrix = self._create_transform_matrix(a0,b0,c0)
        elif self._alignment_atom_indxs:
            a, b, c = self._alignment_atom_indxs
            a0, b0, c0 = molecule.r_array[[a-1,b-1,c-1]]
            t_matrix = self._create_transform_matrix(a0,b0,c0)
            
        molecule.r_array = self._apply_transfom_matrix(t_matrix, molecule.r_array)
        
        self._t_matrix = t_matrix

        return molecule
    
    #instead of from chemlab.notebook import display_molecule to add ball_stick
    def _view_molecule(self, molecule, represent='vdw', colorlist=[]):
        """active representationion of molecule using chemview 
        
        """
        allowed_rep = ['none', 'wire', 'vdw', 'ball_stick']
        if represent not in allowed_rep:
            raise ValueError(
            'unknown molecule representation: {0}, must be in {1}'.format(
            represent, allowed_rep))
        
        topology = {
            'atom_types': molecule.type_array,
            'bonds': molecule.bonds
        }
    
        mv = MolecularViewer(molecule.r_array, topology)
        
        if molecule.n_bonds != 0:
            if represent=='ball_stick':
                mv.ball_and_sticks(colorlist=colorlist)
            elif represent=='wire':
                mv.points(size=0.15, colorlist=colorlist)
                mv.lines(colorlist=colorlist)
            else:
                raise NotImplementedError('none and vdw not implemented in active view')
        else:
            mv.points()
    
        return mv

    def _trim_image(self, im):
        """
        a simple solution to trim whitespace on the image
        
        1. It gets the border colour from the top left pixel, using getpixel, 
        so you don't need to pass the colour.
        2. Subtracts a scalar from the differenced image, 
        this is a quick way of saturating all values under 100, 100, 100 to zero. 
        So is a neat way to remove any 'wobble' resulting from compression.
        """
        bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
        diff = ImageChops.difference(im, bg)
        diff = ImageChops.add(diff, diff, 2.0, -100)
        bbox = diff.getbbox()
        if bbox:
            return im.crop(bbox)

    def _image_molecule(self, molecule, represent='vdw', 
                        colorlist=[], transparent=False,
                         rotation=[0., 0., 0.], width=300, height=300, zoom=1.,
                        lines=[], linestyle='impostors',
                        surfaces=[]):
        """create image of molecule 

        Parameters
        ----------
        molecule : chemlab.core.molecule.Molecule
            the molecule to image
        represent : str
            representation of molecule ('none', 'wire', 'vdw' or 'ball_stick')
        colorlist : list
            color override for each of the atoms (if empty colored by atom type)
        transparent=True : 
            whether atoms should be transparent (based on alpha value) 
        zoom : float
            zoom level of images
        width : int
            width of original images
        height : int
            height of original images (although width takes precedent)
        lines : list
            lines to add to the image in the form; 
            [start_coord, end_coord, start_color, end_color, width, dashed]
        surfaces : list
            surfaces to add to the image in the format;
            [vertices, normals, colors, transparent, wireframe]

        Returns
        -------
        image : PIL.Image
            an image of the system
        
        """
        allowed_rep = ['none', 'wire', 'vdw', 'ball_stick']
        if represent not in allowed_rep:
            raise ValueError(
            'unknown molecule representation: {0}, must be in {1}'.format(
            represent, allowed_rep))
            
        v = QtViewer()
        w = v.widget        
        w.initializeGL()

        if represent=='ball_stick':
            r = v.add_renderer(BallAndStickRenderer,
                                molecule.r_array,
                                molecule.type_array,
                                molecule.bonds,
                                rgba_array=colorlist,
                                linestyle=linestyle,
                                transparent=transparent)
        elif represent=='vdw':
            r = v.add_renderer(AtomRenderer, 
                                molecule.r_array, 
                                molecule.type_array, 
                                rgba_array=colorlist,
                                transparent=transparent)
        elif represent=='wire':
            r = v.add_renderer(WireframeRenderer,
                                molecule.r_array,
                                molecule.type_array,
                                molecule.bonds)
        elif represent=='none':
            r = None
        
        for line in lines:
            #line = [start_coord, end_coord, start_color, end_color, width, dashed]
            #for some reason it didn't like unpacking them to named variables
            v.add_renderer(LineRenderer, [line[0], line[1]], 
                           [[str_to_colour(line[2]), str_to_colour(line[3])]], 
                             width=line[4], dashed=line[5])
        
        for surface in surfaces:
            vertices, normals, colors, transparent, wireframe = surface
            v.add_renderer(TriangleRenderer, vertices, normals, colors, 
                           transparent=transparent,
                           wireframe=wireframe)

        #v.add_post_processing(SSAOEffect)
        w.camera.autozoom(molecule.r_array*1./zoom)
        w.camera.orbit_x(rotation[0]*np.pi/180.)
        w.camera.orbit_y(rotation[1]*np.pi/180.)
        w.camera.orbit_z(rotation[2]*np.pi/180.)
        
        image = w.toimage(width, height)

        # Cleanup
        v.clear()
        del v
        del w
        del r
        
        return self._trim_image(image)

    def _concat_images_horizontal(self, images, gap=10):
        """ concatentate one or more PIL images horizontally 

        Parameters
        ----------
        images : list of pil.Image
            the images to concatenate
        gap : int
            the pixel gap between images
        """
        if len(images) == 1: return images[0]
        
        total_width = sum([img.size[0] for img in images]) + len(images)*gap
        max_height = max([img.size[1] for img in images])
        
        final_img = PIL.Image.new("RGBA", (total_width, max_height), color='white')
        
        horizontal_position = 0
        for img in images:
            final_img.paste(img, (horizontal_position, 0))
            horizontal_position += img.size[0] + gap
        
        return final_img

    def _concat_images_vertical(self, images, gap=10):
        """ concatentate one or more PIL images vertically 

        Parameters
        ----------
        images : list of pil.Image
            the images to concatenate
        gap : int
            the pixel gap between images
        """
        if len(images) == 1: return images[0]
        
        total_width = sum([img.size[0] for img in images]) 
        max_height = max([img.size[1] for img in images]) + len(images)*gap
        
        final_img = PIL.Image.new("RGBA", (total_width, max_height), color='white')
        
        vertical_position = 0
        for img in images:
            final_img.paste(img, (0, vertical_position))
            vertical_position += img.size[1] + gap
        
        return final_img

    def _color_to_transparent(self, image, color=(255, 255, 255)):
        """ sets alpha to 0 for specific colour in PIL image 
        
        Parameters
        ----------
        image : pil.Image
            the images to process
        color : (int, int, int)
            the RGB (0 to 255) color to set alpha to 0
        """
        datas = image.getdata()
    
        newData = []
        for item in datas:
            if item[0] == color[0] and item[1] == color[1] and item[2] == color[2]:
                newData.append((color[0], color[1], color[2], 0))
            else:
                newData.append(item)
    
        image.putdata(newData)
        
        return image

    def _show_molecule(self, molecule, active=False, 
                       represent='vdw', rotations=[[0., 0., 0.]],
                       colorlist=[], transparent=False, 
                       axis_length=0, lines=[], linestyle='impostors',
                       surfaces=[],
                       zoom=1., width=300, height=300, ipyimg=True):
        """show the molecule
        
        Parameters
        ----------
        molecule : chemlab.core.molecule.Molecule
            the molecule to image
        active : bool
            whether the molecule representation should be interactive 
            (ipython notebook only)
        represent : str
            representation of molecule ('none', 'wire', 'vdw' or 'ball_stick')
        colorlist : list
            color override for each of the atoms (if empty colored by atom type)
        transparent=True : 
            whether atoms should be transparent (based on alpha value) 
        axis_length :float
            if non-zero lines will be drawn along each axis to +/- axis_length
        lines : list
            lines to add to the image in the form; 
            [start_coord, end_coord, start_color, end_color, width, dashed]
        surfaces : list
            surfaces to add to the image in the format;
            [vertices, normals, colors, transparent, wireframe]
        zoom : float
            zoom level of images
        width : int
            width of original images
        height : int
            height of original images (although width takes precedent)
        ipyimg : bool
            whether to return an IPython image, PIL image otherwise 

        Returns
        -------
        mol : IPython.display.Image or PIL.Image
            an image of the molecule in the format specified by ipyimg 
            or an active representation
        
        """                
        if active:
            return self._view_molecule(molecule, represent=represent, 
                                          colorlist=colorlist)
        else:
            drawlines=lines[:]
            if axis_length:
                if type(axis_length) is list or type(axis_length) is tuple:
                    neg_length, pos_length = axis_length
                else:
                    neg_length = pos_length = axis_length
                drawlines.append([(-1*neg_length,0,0), (pos_length,0,0), 
                              'red', 'dark_red', 3, True])
                drawlines.append([(0,-1*neg_length,0), (0,pos_length,0), 
                              'light_green', 'dark_green', 3, True])
                drawlines.append([(0,0,-1*neg_length), (0,0,pos_length), 
                              'light_blue', 'dark_blue', 3, True])

            images = []
            for rotation in rotations:
                images.append(self._image_molecule(molecule, 
                                    represent=represent, colorlist=colorlist, 
                                    rotation=rotation, zoom=zoom,
                                    width=width, height=width,
                                    lines=drawlines, linestyle=linestyle,
                                    transparent=transparent, surfaces=surfaces))  
            image = self._concat_images_horizontal(images)
            del images
            
            if ipyimg:
                b = BytesIO()
                image.save(b, format='png')    
                return ipy_Image(data=b.getvalue())
            else:
                return image
                                
    def show_initial(self, gbonds=True, active=False, represent='vdw', 
                     rotations=[[0., 0., 0.]], zoom=1., width=300, height=300,
                     axis_length=0, lines=[], ipyimg=True):
        """show initial geometry (before optimisation) of molecule """
        molecule = self._create_molecule(optimised=False, gbonds=gbonds)
        
        return self._show_molecule(molecule, active=active, 
                                   represent=represent, 
                                   rotations=rotations, zoom=zoom, 
                                   lines=lines, axis_length=axis_length, ipyimg=ipyimg)      
       
    def show_optimisation(self, opt_step=False, gbonds=True, active=False,
                          represent='vdw', rotations=[[0., 0., 0.]], zoom=1.,
                          width=300, height=300, axis_length=0, lines=[], 
                          ipyimg=True):
        """show optimised geometry of molecule """       
        molecule = self._create_molecule(optimised=True, opt_step=opt_step, 
                                         gbonds=gbonds)

        return self._show_molecule(molecule, active=active, 
                                  represent=represent, 
                                  rotations=rotations, zoom=zoom,
                                  lines=lines, axis_length=axis_length,
                                  width=width, height=height, ipyimg=ipyimg)             

    def _rgb_to_hex(self, rgb):
        """convert RGB color to hex format"""        
        return int('0x%02x%02x%02x' % rgb[:3], 16)
       
    def _get_highlight_colors(self, natoms, atomlists, active=False, alpha=0.7):

        norm = mpl.colors.Normalize(vmin=1, vmax=len(atomlists))
        cmap = cm.jet_r
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        
        colorlist = [(211, 211, 211, int(255*alpha)) for n in range(natoms)]
        
        for n in range(natoms):
            for group, atomlist in enumerate(atomlists):
                if n+1 in atomlist:
                    colorlist[n] = m.to_rgba(group+1, bytes=True)
                    break
          
        if active:           
            colorlist = [self._rgb_to_hex(col) for col in colorlist]
        
        return colorlist

    def show_highlight_atoms(self, atomlists, transparent=False, alpha=0.7,
                             gbonds=True, active=False, optimised=True,
                        represent='vdw', rotations=[[0., 0., 0.]], zoom=1.,
                        width=300, height=300, axis_length=0, lines=[], ipyimg=True):
               
        if optimised:
            natoms = self._read_data('_opt_data', 'natom')        
        else:
            natoms = self._read_data('_init_data', 'natom')
        
        atomlists=[self._atom_groups[grp] if type(grp) is str else grp for grp in atomlists]

        colorlist = self._get_highlight_colors(natoms, atomlists, active,
                                               alpha=alpha)
        
        molecule = self._create_molecule(optimised=optimised, gbonds=gbonds)

        if transparent:
            linestyle='lines'
        else:
            linestyle='impostors'
            
        return self._show_molecule(molecule, active=active, 
                                   transparent=transparent,
                                  represent=represent, 
                                  rotations=rotations, zoom=zoom,
                                  colorlist=colorlist, linestyle=linestyle,
                                  lines=lines, axis_length=axis_length,
                                  width=width, height=height, ipyimg=ipyimg) 
                                  
    def _write_init_file(self, molecule, file_name, descript='', 
                         overwrite=False, decimals=8,
                         charge=0, multiplicity=1,
                         folder_obj=None):
        """ write a template gaussian input file to folder
        
                
        """
        if not type(charge) is int or not type(multiplicity) is int:
            raise ValueError('charge and multiplicity of molecule must be defined')
        
        if not folder_obj:
            folder_obj = self._folder           
            
        with folder_obj as folder:
            with folder.write_file(file_name+'_init.com', overwrite) as f:
                f.write('%chk={0}_init.chk \n'.format(file_name))
                f.write('# opt b3lyp/3-21g \n')
                f.write('\n')
                f.write('{0} \n'.format(descript))
                f.write('\n')
                f.write('{0} {1} \n'.format(charge, multiplicity))
                for t, c in zip(molecule.type_array, molecule.r_array*10.): # nanometers to angstrom
                    x, y, z = c.round(decimals)
                    f.write(' {0}\t{1}\t{2}\t{3} \n'.format(t, x, y, z))
                f.write('\n')
        
        return True

    def _array_transformation(self, array, rotations, transpose=[0,0,0]):
        """ 3D rotation around x-axis, then y-axis, then z-axis,
            then transposition """
        if rotations == [0,0,0]: 
            new = array
        else:            
            x, y, z = rotations
            rot_x = rotation_matrix(x*np.pi/180., [1, 0, 0])[:3,:3]
            rot_y = rotation_matrix(y*np.pi/180., [0, 1, 0])[:3,:3]
            rot_z = rotation_matrix(z*np.pi/180., [0, 0, 1])[:3,:3]
            
            rot = np.dot(rot_z, np.dot(rot_y, rot_x))
            
            new = np.array([np.dot(rot, coord) for coord in array])   
        
        new[:,0] += transpose[0]
        new[:,1] += transpose[1]
        new[:,2] += transpose[2]
        
        return new
        
    def combine_molecules(self, other_mol, self_atoms=False, other_atoms=False,
                          self_rotation=[0,0,0], other_rotation=[0,0,0],
                          self_transpose=[0,0,0], other_transpose=[0,0,0],
                          self_opt=True, other_opt=True,
                          charge=None, multiplicity=None,
                          out_name=False, descript='', overwrite=False,
                          active=False,
                          represent='ball_stick', rotations=[[0., 0., 0.]], zoom=1.,
                          width=300, height=300, axis_length=0, ipyimg=True,
                          folder_obj=None):
        """ transpose in nanometers """                      
        mol = self._create_molecule(optimised=self_opt)
        if self_atoms:
            if type(self_atoms) is str:
                self_atoms = self._atom_groups[self_atoms]
            self_indxs = np.array(self_atoms) - 1
            mol.r_array = mol.r_array[self_indxs]
            mol.type_array = mol.type_array[self_indxs]
        mol.r_array = self._array_transformation(mol.r_array, 
                                                 self_rotation, self_transpose)
        
        mol_atoms = [i+1 for i in range(len(mol.type_array))]
        
        other = other_mol._create_molecule(optimised=other_opt)
        if other_atoms:
            if type(other_atoms) is str:
                other_atoms = other_mol._atom_groups[other_atoms]
            other_indxs = np.array(other_atoms) - 1
            other.r_array = other.r_array[other_indxs]
            other.type_array = other.type_array[other_indxs]
        other.r_array = self._array_transformation(other.r_array, 
                                                   other_rotation, other_transpose)

        other_atoms = [i+1+len(mol.type_array) for i in range(len(other.type_array))]        
        
        mol.r_array = np.concatenate([mol.r_array, other.r_array])
        mol.type_array = np.concatenate([mol.type_array, other.type_array])
        mol.guess_bonds()
        
        if out_name:
            self._write_init_file(mol, out_name, descript, overwrite,
                                  charge=charge, multiplicity=multiplicity,
                                  folder_obj=folder_obj)
            
        colorlist = self._get_highlight_colors(len(mol.type_array), 
                                               [mol_atoms, other_atoms], active)
            
        return self._show_molecule(mol, active=active, 
                                  represent=represent, 
                                  rotations=rotations, zoom=zoom,
                                  colorlist=colorlist,
                                  axis_length=axis_length,
                                  width=width, height=height, ipyimg=ipyimg,
                                  ) 
        
    def _get_charge_colors(self, relative=False, minval=-1, maxval=1, alpha=None):
        
        charges = self._read_data('_nbo_data', 'atomcharges')['natural']
        if relative: minval, maxval = (min(charges), max(charges))
        norm = mpl.colors.Normalize(vmin=minval, vmax=maxval)
        cmap = cm.bwr
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors=m.to_rgba(charges, alpha=alpha, bytes=True)
        
        return colors

    def show_nbo_charges(self, gbonds=True, active=False, 
                         relative=False, minval=-1, maxval=1,
                         represent='vdw', rotations=[[0., 0., 0.]], zoom=1.,
                         width=300, height=300, axis_length=0, lines=[], ipyimg=True):
        
        colorlist = self._get_charge_colors(relative, minval, maxval)

        molecule = self._create_molecule(optimised=True, gbonds=gbonds)

        return self._show_molecule(molecule, active=active, 
                                  represent=represent, 
                                  rotations=rotations, zoom=zoom,
                                  colorlist=colorlist,
                                  lines=lines, axis_length=axis_length,
                                  width=width, height=height, ipyimg=ipyimg) 
    
    def get_orbital_count(self):
        """return number of orbitals """
        moenergies = self._read_data('_nbo_data', "moenergies")[0]
        return int(moenergies.shape[0])

        
    def _find_nearest_above(self, my_array, target):
        diff = my_array - target
        mask = np.ma.less_equal(diff, 0)
        # We need to mask the negative differences and zero
        # since we are looking for values above
        if np.all(mask):
            return None # returns None if target is greater than any value
        masked_diff = np.ma.masked_array(diff, mask)
        return masked_diff.argmin()
    
    def _find_nearest_below(self, my_array, target):
        diff = my_array - target
        mask = np.ma.greater_equal(diff, 0)
        # We need to mask the positive differences and zero
        # since we are looking for values above
        if np.all(mask):
            return None # returns None if target is lower than any value
        masked_diff = np.ma.masked_array(diff, mask)
        return masked_diff.argmax()

    def get_orbital_homo_lumo(self):
        """return orbital numbers of homo and lumo """
        moenergies = self._read_data('_nbo_data', "moenergies")[0]
        homo = self._find_nearest_below(moenergies, 0.) + 1
        lumo = self._find_nearest_above(moenergies, 0.) + 1
        
        return homo, lumo
                
    def get_orbital_energies(self, orbitals, eunits='eV'):
        """the orbital energies for listed orbitals
        
        Parameters
        ----------
        orbitals : int or iterable of ints
            the orbital(s) to return energies for (starting at 1) 
        eunits : str
            the units of energy
        Returns
        -------
        moenergies : np.array
            energy for each orbital
        """
        orbitals = np.array(orbitals, ndmin=1, dtype=int)
        assert np.all(orbitals>0) and np.all(orbitals<=self.get_orbital_count()), (
            'orbitals must be in range 1 to number of orbitals')
        moenergies = self._read_data('_nbo_data', "moenergies")[0]
        if not eunits=='eV':
            moenergies = convertor(moenergies, 'eV', eunits)
        return moenergies[orbitals-1]
        
    def yield_orbital_images(self, orbitals, iso_value=0.02, extents=(2,2,2),
                     transparent=True, alpha=0.5, wireframe=True,
                     bond_color=(255, 0, 0), antibond_color=(0, 255, 0),
                     resolution=100, gbonds=True, represent='ball_stick', 
                     rotations=[[0., 0., 0.]], zoom=1.,
                     width=300, height=300, axis_length=0, lines=[], ipyimg=True):
        """yield orbital images

        Parameters
        ----------
        orbitals : int or list of ints
            the orbitals to show (in range 1 to number of orbitals)
        iso_value : float
            The value for which the function should be constant.
        extents : (float, float, float)
            +/- x,y,z to extend the molecule geometrty when constructing the surface
        transparent=True : 
            whether iso-surface should be transparent (based on alpha value) 
        alpha : 
            alpha value of iso-surface
        wireframe : 
            whether iso-surface should be wireframe (or solid)
        bond_color : 
            color of bonding orbital surface in RGB format
        antibond_color : 
            color of anti-bonding orbital surface in RGB format
        resolution : int
            The number of grid point to use for the surface. An high value will 
            give better quality but lower performance.
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
        lines : [start_coord, end_coord, start_color, end_color, width, dashed]
            lines to add to image
        ipyimg : bool
            whether to return an IPython image, PIL image otherwise 

        Returns
        -------
        mol : IPython.display.Image or PIL.Image
            an image of the molecule in the format specified by ipyimg 
            
        """
        warnings.warn('Orbitals are currently an experimental feature')
        
        orbitals = np.array(orbitals, ndmin=1, dtype=int)
        assert np.all(orbitals>0) and np.all(orbitals<=self.get_orbital_count()), (
            'orbitals must be in range 1 to number of orbitals')
        
        r, g, b = bond_color
        bond_rgba = (r, g, b, int(255*alpha))
        r, g, b = antibond_color
        antibond_rgba = (r, g, b, int(255*alpha))

        #To fix issue with rotations
        #TODO could probably do this better (no self._t_matrix)
        alignto = self._alignment_atom_indxs[:]
        self._alignment_atom_indxs = None                
        r_array = self._create_molecule(optimised=True).r_array
        self._alignment_atom_indxs = alignto

        molecule = self._create_molecule(optimised=True, gbonds=gbonds)
        mocoeffs = self._read_data('_nbo_data', "mocoeffs")
        gbasis = self._read_data('_nbo_data', "gbasis")

        for orbital in orbitals:
            coefficients = mocoeffs[0][orbital-1]
            f = molecular_orbital(r_array.astype('float32'), 
                                  coefficients.astype('float32'), 
                                  gbasis)
                
            surfaces = []
            b_iso = get_isosurface(r_array, f, iso_value, bond_rgba,
                                   resolution=resolution)
            if b_iso:
                verts, normals, colors = b_iso
                verts = self._apply_transfom_matrix(self._t_matrix, verts)
                normals = self._apply_transfom_matrix(self._t_matrix, normals)
                surfaces.append([verts, normals, colors, transparent, wireframe])                                        
            a_iso = get_isosurface(molecule.r_array, f, -iso_value, antibond_rgba,
                                   resolution=resolution)
            if a_iso:
                averts, anormals, acolors = a_iso
                averts = self._apply_transfom_matrix(self._t_matrix, averts)
                anormals = self._apply_transfom_matrix(self._t_matrix, anormals)
                surfaces.append([averts, anormals, acolors, transparent,wireframe])                                        
            
            yield self._show_molecule(molecule, 
                                      represent=represent, 
                                      rotations=rotations, zoom=zoom,
                                      surfaces=surfaces, transparent=False,
                                      lines=lines, axis_length=axis_length,
                                      width=width, height=height, ipyimg=ipyimg) 

    def show_active_orbital(self, orbital, iso_value=0.03, alpha=0.5,
                            bond_color=(255, 0, 0), antibond_color=(0, 255, 0),
                            gbonds=True):
        """get interactive representation of orbital
        
        Parameters
        ----------
        orbital : int
            the orbital to show (in range 1 to number of orbitals)
        iso_value : float
            The value for which the function should be constant.
        alpha : 
            alpha value of iso-surface
        bond_color : 
            color of bonding orbital surface in RGB format
        antibond_color : 
            color of anti-bonding orbital surface in RGB format
        gbonds : bool
            guess bonds between atoms (via distance)
        """
        orbital = np.array(orbital, ndmin=1, dtype=int)
        assert np.all(orbital>0) and np.all(orbital<=self.get_orbital_count()), (
            'orbital must be in range 1 to number of orbitals')
        orbital = orbital[0]
        r, g, b = bond_color
        bond_rgba = (r, g, b, int(255*alpha))
        r, g, b = antibond_color
        antibond_rgba = (r, g, b, int(255*alpha))

        molecule = self._create_molecule(optimised=True, gbonds=gbonds)
        mocoeffs = self._read_data('_nbo_data', "mocoeffs")
        gbasis = self._read_data('_nbo_data', "gbasis")

        coefficients = mocoeffs[0][orbital-1]
        f = molecular_orbital(molecule.r_array.astype('float32'), 
                              coefficients.astype('float32'), 
                              gbasis)

        mv = MolecularViewer(molecule.r_array, { 'atom_types': molecule.type_array,
                                                 'bonds': molecule.bonds })
        mv.wireframe()
        mv.add_isosurface(f, isolevel=iso_value, color=self._rgb_to_hex(bond_rgba))
        mv.add_isosurface(f, isolevel=-iso_value, color=self._rgb_to_hex(antibond_rgba))
        return mv

    def _converter(self, val, unit1, unit2):
        
        multiple = {('nm', 'nm') : 1.,
                    ('nm', 'Angstrom') : 0.1}

        return val * multiple[(unit1, unit2)]      

    def calc_min_dist(self, idx_list1, idx_list2, optimisation=True, units='nm',
                      ignore_missing=True):
        """ indexes start at 1 """

        if optimisation:
            molecule = self._read_data('_opt_data', 'molecule')  
        else:
            molecule = self._read_data('_init_data', 'molecule')
            
        if type(idx_list1) is str:
            idx_list1 = self._atom_groups[idx_list1]
        if type(idx_list2) is str:
            idx_list2 = self._atom_groups[idx_list2]
        
        # remove atoms not in molecule
        if ignore_missing:
            idx_list1 = [idx for idx in idx_list1[:] if idx <= molecule.n_atoms]
            idx_list2 = [idx for idx in idx_list2[:] if idx <= molecule.n_atoms]
        
            if not idx_list1 or not idx_list2:
                return np.nan

        indx_combis = cartesian([idx_list1, idx_list2])
        c1 = molecule.r_array[indx_combis[:, 0]-1]
        c2 = molecule.r_array[indx_combis[:, 1]-1]

        dist =  np.min(np.linalg.norm(c1-c2, axis=1))       
        
        return self._converter(dist, 'nm', units)
                                  
    def calc_bond_angle(self, indxs, optimisation=True, mol=None):
        """ Returns the angle in degrees between three points    """

        if mol:
            molecule = mol
        elif optimisation:
            molecule = self._read_data('_opt_data', 'molecule')  
        else:
            molecule = self._read_data('_init_data', 'molecule')

        v1 = molecule.r_array[indxs[0]-1] - molecule.r_array[indxs[1]-1]
        v2 = molecule.r_array[indxs[2]-1] - molecule.r_array[indxs[1]-1]
        cosang = np.dot(v1, v2)
        sinang = np.linalg.norm(np.cross(v1, v2))
        
        return np.degrees(np.arctan2(sinang, cosang))

    def calc_dihedral_angle(self, indxs, optimisation=True, mol=None):
        """ Returns the angle in degrees between four points  """

        if mol:
            molecule = mol
        elif optimisation:
            molecule = self._read_data('_opt_data', 'molecule')  
        else:
            molecule = self._read_data('_init_data', 'molecule')

        p = np.array([molecule.r_array[indxs[0]-1], molecule.r_array[indxs[1]-1], 
                      molecule.r_array[indxs[2]-1], molecule.r_array[indxs[3]-1]])
        b = p[:-1] - p[1:]
        b[0] *= -1
        v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
        # Normalize vectors
        v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
        b1 = b[1] / np.linalg.norm(b[1])
        x = np.dot(v[0], v[1])
        m = np.cross(v[0], b1)
        y = np.dot(m, v[1])
        angle = np.degrees(np.arctan2( y, x ))
        
        return angle #np.mod(angle, 360)
                                  
    def calc_polar_coords_from_plane(self, p1, p2, p3, c, optimisation=True, 
                                   units='nm'):
        """ returns the distance r and angles theta, phi of atom c 
        to the circumcenter of the plane formed by [p1, p2, p3]
        
        the plane formed will have;
            x-axis along p1, 
            y-axis anticlock-wise towards p2,
            z-axis normal to the plane
        
        theta (azimuth) is the in-plane angle from the x-axis towards the y-axis
        phi (inclination) is the out-of-plane angle from the x-axis towards 
        the z-axis
        """
        
        alignto = self._alignment_atom_indxs[:]    
        self._alignment_atom_indxs = (p1, p2, p3)
        
        if optimisation:
            molecule = self._create_molecule(optimised=True)
        else:
            molecule = self._create_molecule(optimised=False)
        
        if len(molecule.r_array)<c:
            self._alignment_atom_indxs = alignto
            return np.nan, np.nan, np.nan
            
        x, y, z = molecule.r_array[c-1]
        
        r = self._converter(sqrt(x*x+y*y+z*z), 'nm', units)
        theta = degrees(atan2(y, x))
        phi = degrees(atan2(z, x))
                
        self._alignment_atom_indxs = alignto
        
        return r, theta, phi  

    def calc_2plane_angle(self, p1, p2, optimisation=True):
        """return angle of planes """
        a1, a2, a3 = p1
        b1, b2, b3 = p2        
        
        if optimisation:
            molecule = self._read_data('_opt_data', 'molecule')  
        else:
            molecule = self._read_data('_init_data', 'molecule')

        v1a = molecule.r_array[a2-1] - molecule.r_array[a1-1]
        v2a = molecule.r_array[a3-1] - molecule.r_array[a1-1]

        v1b = molecule.r_array[b2-1] - molecule.r_array[b1-1]
        v2b = molecule.r_array[b3-1] - molecule.r_array[b1-1]
        
        vnormala = np.cross(v1a,v2a)
        vnormalb = np.cross(v1b,v2b)
        
        cos_theta = np.dot(vnormala, vnormalb)/(
                    np.linalg.norm(vnormala)*np.linalg.norm(vnormalb))
        #account for rounding errors
        if cos_theta > 1.: cos_theta = 1.
        if cos_theta < -1.: cos_theta = -1.                   
        
        return degrees(acos(cos_theta))

    def calc_opt_trajectory(self, atom, plane=[]):
        """ calculate the trajectory of an atom as it is optimised,
        relative to a plane of three atoms """
        alignto = self._alignment_atom_indxs[:]
        self._alignment_atom_indxs = plane

        #get coord from init
        mol = self._create_molecule(optimised=False)
        init = mol.r_array[atom-1]
        #get coords from opt
        opts=[]
        for data in self._prev_opt_data + [self._opt_data]:
            run = []
            for n in range(len(data.read('atomcoords'))):
                mol = self._create_molecule(data=data, opt_step=n)
                run.append(mol.r_array[atom-1])
            opts.append(np.array(run))
        
        self._alignment_atom_indxs = alignto
        
        return init, opts

    _SUFFIXES = {1: 'st', 2: 'nd', 3: 'rd'}
    def _ordinal(self, num):
        # I'm checking for 10-20 because those are the digits that
        # don't follow the normal counting scheme. 
        if 10 <= num % 100 <= 20:
            suffix = 'th'
        else:
            # the second parameter is a default.
            suffix = self._SUFFIXES.get(num % 10, 'th')
        return str(num) + suffix

    def plot_opt_trajectory(self, atom, plane=[], ax_lims=None, ax_labels=False):
        """plot the trajectory of an atom as it is optimised,
        relative to a plane of three atoms """        
        init, opts = self.calc_opt_trajectory(atom, plane)
  
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(init[0], init[1], init[2], c='r', 
                   s=30, label='Initial Position')
        ax.scatter(opts[-1][-1,0], opts[-1][-1,1], opts[-1][-1,2], c=['g'], 
                   s=30, label='Optimised Position')
        for i, opt in enumerate(opts):
            ax.plot(opt[:,0], opt[:,1], opt[:,2], 
                    label='{0} optimisation'.format(self._ordinal(i+1)))

        mol = self._create_molecule().r_array  
        a,b,c=plane
        ax.scatter(*mol[a-1], c='k', marker='^', s=30, label='Atom {0}'.format(a))
        ax.scatter(*mol[b-1], c='k', marker='o', s=30, label='Atom {0}'.format(b))
        ax.scatter(*mol[c-1], c='k', marker='s', s=30, label='Atom {0}'.format(c))

        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  
        
        if ax_lims:
            x, y, z = ax_lims
            ax.set_xlim3d(-x, x)
            ax.set_ylim3d(-y, y)
            ax.set_zlim3d(-z, z)
        
        if ax_labels:        
            ax.set_xlabel('x (nm)')
            ax.set_ylabel('y (nm)')
            ax.set_zlabel('z (nm)')
    
        return ax
        
    def calc_nbo_charge(self, atoms=[]):
        """ returns total charge of the atoms """
        charges = self._read_data('_nbo_data', 'atomcharges')['natural']
        if atoms==[]: 
            return np.sum(charges)
            
        if type(atoms) is str:
            atoms = self._atom_groups[atoms]
            
        atoms = np.array(atoms) -1 # 1->0 base        
        try:
            subcharges = charges[atoms]
        except IndexError:
            return np.nan
        
        return np.sum(subcharges)
        
    def calc_nbo_charge_center(self, p1, p2, p3, positive=True, units='nm', 
                               atoms=[]):
        """ returns the distance r amd angles theta, phi of the positive/negative
        charge center to the circumcenter of the plane formed by [p1, p2, p3]
        
        the plane formed will have;
            x-axis along p1, 
            y-axis anticlock-wise towards p2,
            z-axis normal to the plane
        
        theta (azimuth) is the in-plane angle from the x-axis towards the y-axis
        phi (inclination) is the out-of-plane angle from the x-axis towards 
        the z-axis
        """

        molecule = self._create_molecule(alignment_atoms=(p1, p2, p3))
        charges = self._read_data('_nbo_data', 'atomcharges')['natural']
        coords = molecule.r_array   
        
        if type(atoms) is str:
            atoms = self._atom_groups[atoms]

        if atoms:
            atoms = np.array(atoms) -1 # 1->0 base
            charges = charges[atoms]
            coords = coords[atoms]

        if positive:
            weighted_coords = charges[charges>0] * coords[charges>0].T            
        else:
            weighted_coords = -1*charges[charges<0] * coords[charges<0].T
        charge_center = np.mean(weighted_coords.T, axis=0)
        x, y, z = charge_center
        
        r = self._converter(sqrt(x*x+y*y+z*z), 'nm', units)
        theta = degrees(atan2(y, x))
        phi = degrees(atan2(z, x))
        
        return r, theta, phi                     

    def get_sopt_analysis(self, eunits='kJmol-1', atom_groups=[], 
                          charge_info=False):
        """interactions between "filled" (donor) Lewis-type 
        Natural Bonding Orbitals (NBOs) and "empty" (acceptor) non-Lewis NBOs,
        using Second Order Perturbation Theory (SOPT)
        
        Parameters
        ----------
        eunits : str
            the units of energy to return
        atom_groups : [list or str, list or str]
            restrict interactions to between two lists (or identifiers) of atom indexes
        charge_info : bool
            include charge info for atoms 
            (under headings 'A_Charges' and 'D_Charges')
        Returns
        -------
        analysis : pandas.DataFrame
            a table of interactions
        """
        
        #sopt = copy.deepcopy(self._read_data('_nbo_data', 'sopt'))
        sopt = self._read_data('_nbo_data', 'sopt')
        
        
        df = pd.DataFrame(sopt, 
                          columns=['Dtype', 'Donors', 'Atype', 'Acceptors', 'E2'])
                        
        if atom_groups:
            
            group1, group2 = atom_groups

            if type(group1) is str:
                group1 = self._atom_groups[group1]
            if type(group2) is str:
                group2 = self._atom_groups[group2]
        
            match_rows=[]
            for indx, rw in df.iterrows():
                if set(group1).issuperset(rw.Acceptors) and set(group2).issuperset(rw.Donors):
                    match_rows.append(rw)
                elif set(group2).issuperset(rw.Acceptors) and set(group1).issuperset(rw.Donors):
                    match_rows.append(rw)
            
            df = pd.DataFrame(match_rows)
                    
        if not eunits=='kcal': 
            df.E2 = convertor(df.E2, 'kcal', eunits)

        typ = self._read_data('_nbo_data', 'molecule').type_array        
        df['D_Symbols'] = df.Donors.apply(lambda x: [typ[i-1] for i in x])
        df['A_Symbols'] = df.Acceptors.apply(lambda x: [typ[i-1] for i in x])

        if charge_info:        
            chrg= self._read_data('_nbo_data', 'atomcharges')['natural']
            df['D_Charges'] = df.Donors.apply(lambda x: [chrg[i-1] for i in x])
            df['A_Charges'] = df.Acceptors.apply(lambda x: [chrg[i-1] for i in x])

            return df[['Dtype', 'Donors', 'D_Symbols', 'D_Charges', 
                       'Atype', 'Acceptors', 'A_Symbols', 'A_Charges', 
                       'E2']] 
        else:
            return df[['Dtype', 'Donors', 'D_Symbols', 
                       'Atype', 'Acceptors', 'A_Symbols', 
                       'E2']] 
            
        
    def get_hbond_analysis(self, min_energy=0., atom_groups=[], eunits='kJmol-1'):
        """EXPERIMENTAL! hydrogen bond analysis (DH---A), 
        using Second Order Bond Perturbation Theiry

        Parameters
        ----------
        min_energy : float
            the minimum interaction energy to report
        eunits : str
            the units of energy to return
        atom_groups : [list or str, list or str]
            restrict interactions to between two lists (or identifiers) of atom indexes
        Returns
        -------
        analysis : pandas.DataFrame
            a table of interactions

        uses a strict definition of a hydrogen bond as:
        interactions between "filled" (donor) Lewis-type Lone Pair (LP) NBOs 
        and "empty" (acceptor) non-Lewis Bonding (BD) NBOs
        """

        df = self.get_sopt_analysis(atom_groups=atom_groups, eunits=eunits)
        df = df[df.E2 >= min_energy]
        df = df[df.A_Symbols.apply(lambda x: 'H' in x) & 
                df.Dtype.str.contains('LP') & 
                df.Atype.str.contains('BD*')]
        
        return df

    def calc_sopt_energy(self, atom_groups=[], eunits='kJmol-1', no_hbonds=False):
        """calculate total energy of interactions between "filled" (donor) Lewis-type 
        Natural Bonding Orbitals (NBOs) and "empty" (acceptor) non-Lewis NBOs,
        using Second Order Perturbation Theory 

        Parameters
        ----------
        eunits : str
            the units of energy to return
        atom_groups : [list or str, list or str]
            restrict interactions to between two lists (or identifiers) of atom indexes
        no_hbonds : bool
            whether to ignore H-Bonds in the calculation
        Returns
        -------
        analysis : pandas.DataFrame
            a table of interactions
        """
        df = self.get_sopt_analysis(atom_groups=atom_groups, eunits=eunits)
        
        if no_hbonds:
            dfh = df[df.A_Symbols.apply(lambda x: 'H' in x) & 
                    df.Dtype.str.contains('LP') & 
                    df.Atype.str.contains('BD*')]
            df = df.loc[set(df.index).difference(dfh.index)]          
        
        return df.E2.sum()

    def show_sopt_bonds(self, min_energy=20., cutoff_energy=0., atom_groups=[],
                        bondwidth=5, eunits='kJmol-1', no_hbonds=False,
                        gbonds=True, active=False,
                        represent='ball_stick', rotations=[[0., 0., 0.]], zoom=1.,
                        width=300, height=300, axis_length=0, lines=[],
                        relative=False, minval=-1, maxval=1,
                        alpha=0.5, transparent=True,
                        ipyimg=True):
        """visualisation of interactions between "filled" (donor) Lewis-type 
        Natural Bonding Orbitals (NBOs) and "empty" (acceptor) non-Lewis NBOs,
        using Second Order Perturbation Theory
                
        """
        df = self.get_sopt_analysis(atom_groups=atom_groups, eunits=eunits)
        df = df[df.E2 >= min_energy]
        
        if no_hbonds:
            dfh = self.get_hbond_analysis(min_energy=min_energy, eunits=eunits,
                                          atom_groups=atom_groups)
            df = df.loc[set(df.index).difference(dfh.index)]  
        
        molecule = self._create_molecule(gbonds=gbonds)
        
        drawlines = lines[:]
        for i, rw in df.iterrows():
            d_coord = np.mean([molecule.r_array[d-1] for d in rw.Donors], axis=0)
            a_coord = np.mean([molecule.r_array[a-1] for a in rw.Acceptors], axis=0)
            
            dashed = rw.E2 < cutoff_energy
            drawlines.append([d_coord, a_coord, 'blue', 'red', 
                              max([1, bondwidth-1]), dashed])
        
        colorlist = self._get_charge_colors(relative, minval, maxval, alpha=alpha)
        
        return self._show_molecule(molecule, active=active, 
                                  represent=represent, 
                                  rotations=rotations, zoom=zoom,
                                  colorlist=colorlist,
                                  lines=drawlines, axis_length=axis_length,
                                  width=width, height=height, linestyle='lines', 
                                  transparent=transparent,
                                  ipyimg=ipyimg) 
    
    def calc_hbond_energy(self, atom_groups=[], eunits='kJmol-1'):
        
        df = self.get_hbond_analysis(atom_groups=atom_groups, eunits=eunits)
        
        return df.E2.sum()
    
    def show_hbond_analysis(self, min_energy=0., atom_groups=[], 
                        cutoff_energy=0., eunits='kJmol-1', bondwidth=5, 
                        gbonds=True, active=False,
                        represent='ball_stick', rotations=[[0., 0., 0.]], zoom=1.,
                        width=300, height=300, axis_length=0, lines=[],
                        relative=False, minval=-1, maxval=1,
                        alpha=0.5, transparent=True, ipyimg=True):
        """EXPERIMENTAL! hydrogen bond analysis DH---A
        
        For a hydrogen bond to occur there must be both a hydrogen donor and an 
        acceptor present. The donor in a hydrogen bond is the atom to which the 
        hydrogen atom participating in the hydrogen bond is covalently bonded, 
        and is usually a strongly electronegative atom such as N, O, or F. The 
        hydrogen acceptor is the neighboring electronegative ion or molecule, 
        and must posses a lone electron pair in order to form a hydrogen bond.
        
        Since the hydrogen donor is strongly electronegative, it pulls the 
        covalently bonded electron pair closer to its nucleus, and away from 
        the hydrogen atom. The hydrogen atom is then left with a partial 
        positive charge, creating a dipole-dipole attraction between the 
        hydrogen atom bonded to the donor, and the lone electron pair on the acceptor.
        """
        df = self.get_hbond_analysis(min_energy=min_energy, eunits=eunits, 
                                     atom_groups=atom_groups)
        
        molecule = self._create_molecule(gbonds=gbonds)
        
        drawlines = lines[:]
        for i, rw in df.iterrows():
            d_coord = np.mean([molecule.r_array[d-1] for d in rw.Donors], axis=0)
            
            h_indx = rw.A_Symbols.index('H')
            a_coord = molecule.r_array[rw.Acceptors[h_indx]-1]
            
            dashed = rw.E2 < cutoff_energy
            drawlines.append([d_coord, a_coord, 'blue', 'red', 
                              max([1, bondwidth-1]), dashed])
        
        colorlist = self._get_charge_colors(relative, minval, maxval, alpha=alpha)
        
        return self._show_molecule(molecule, active=active, 
                                  represent=represent, 
                                  rotations=rotations, zoom=zoom,
                                  colorlist=colorlist,
                                  lines=drawlines, axis_length=axis_length,
                                  width=width, height=height, linestyle='lines', 
                                  transparent=transparent,
                                  ipyimg=ipyimg) 

    def _img_to_plot(self, x, y, image, ax=None, zoom=1):
        """add image to matplotlib axes at (x,y) """
        if ax is None:
            ax = plt.gca()
        im = OffsetImage(image, zoom=zoom)
        artists = []
        ab = AnnotationBbox(im, (x, y), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
        #ax.update_datalim(np.column_stack([x, y]))
        ax.autoscale(tight=False)
        return artists
    
    # TODO get fixed atoms from scan file
    def plot_pes_scans(self, fixed_atoms, eunits='kJmol-1', 
                       img_pos='', rotation=[0., 0., 0.], zoom=1, order=1):
        """plot Potential Energy Scan

        img_pos : <'','local_mins','local_maxs','global_min','global_max'>
            position image(s) of molecule conformation(s) on plot
        rotation : [float, float, float]
            rotation of molecule image(s)
        """
        scan_datas = self._pes_data
        
        if len(fixed_atoms) == 4:
            xlabel = 'Dihedral Angle'
            func = self.calc_dihedral_angle
        elif len(fixed_atoms)==3:
            xlabel = 'Valence Angle'
            func = self.calc_bond_angle
        else:
            raise Exception('not 3 or 4 fixed atoms')
            
        angles = []
        energies = []
        mols = []
        for scan in scan_datas:
            for i in range(scan.read('nscans')):
                mol = scan.read('molecule', scan=i)
                mols.append(self._create_molecule(data=scan, scan_step=i))
                angles.append(func(fixed_atoms, mol=mol))
            energies.extend(scan.read('scanenergies'))
        
        # remove duplicate angles and sort by angle 
        # so that the local max are found correctly
        df = pd.DataFrame({'energy':convertor(np.array(energies), 'eV', eunits), 
                           'angle':angles, 'mol':mols})
        df['rounded'] = df.angle.round(2) #rounding errors?
        df.drop_duplicates('rounded', inplace=True)
        df.sort('angle', inplace=True)
                
        angles = np.array(df.angle.tolist())
        energies = np.array(df.energy.tolist())
        
        fig, ax = plt.subplots()
        ax.plot(angles, energies)
        ax.scatter(angles, energies)
        ax.set_ylabel('Energy ({0})'.format(eunits))
        ax.set_xlabel(xlabel)
            
        feature_dict = {
        '':[],
        'local_maxs' : df.index[argrelextrema(energies, np.greater, mode='wrap', order=order)[0]],
        'local_mins' : df.index[argrelextrema(energies, np.less, mode='wrap', order=order)[0]],
        'global_min' : [df.energy.idxmin()],
        'global_max' : [df.energy.idxmax()]}
         
        for indx in feature_dict[img_pos]:
            img = self._image_molecule(df.mol.loc[indx], rotation=rotation, represent='ball_stick')                        
            img = self._color_to_transparent(img)
            self._img_to_plot(df.angle.loc[indx], df.energy.loc[indx], img, zoom=zoom, ax=ax)
            
        return ax
