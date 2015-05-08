# -*- coding: utf-8 -*-
"""
Created on Fri May 01 21:24:31 2015

@author: chris
"""
import os
from io import BytesIO
import PIL
from PIL import Image, ImageChops
from types import MethodType

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
from chemlab.graphics.transformations import rotation_matrix
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

from .chemlab_patch.graphics.renderers.atom import AtomRenderer
from .chemlab_patch.graphics.renderers.ballandstick import BallAndStickRenderer
from .chemlab_patch.graphics.renderers.line import LineRenderer
from chemlab.graphics.postprocessing import SSAOEffect # Screen Space Ambient Occlusion
from chemlab.utils import cartesian
from cclib.parser.utils import convertor

from chemlab.graphics.colors import get as str_to_colour

#instead of chemview MolecularViewer to add defined colouring
#also ignore; 'FutureWarning: IPython widgets are experimental and may change in the future.'
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from .chemview_patch.viewer import MolecularViewer

from .circumcenter import circumcenter

class Molecule:
    def __init__(self, folder, init_fname=False, opt_fname=False, 
                 freq_fname=False, nbo_fname=False, 
                 pes_fname=False,
                 alignto=[]):
        """a class to contain gaussian input/output 
        for a single molecular geometry 
        
        NB: nbo population analysis must be run with the GFInput flag to ensure 
        data is output to the log file 
        """
        self._folder = folder
        self._init_data = None
        self._prev_opt_data = []
        self._opt_data = None
        self._freq_data = None
        self._nbo_data = None
        self._pes_data = []
        self._alignment_atom_indxs = ()
        if alignto: self.set_alignment_atoms(*alignto)
        
        if init_fname: self.add_initialgeom(init_fname)        
        if opt_fname: self.add_optimisation(opt_fname)
        if freq_fname: self.add_frequency(freq_fname)
        if nbo_fname: self.add_nbo_analysis(nbo_fname)
        if pes_fname: self.add_pes_analysis(pes_fname)
    
    def __repr__(self):
        return '<PyGauss Molecule>'
    
    def _get_data(self, log_file, ftype='gaussian'):
        
        gaussian_handler = _create_cclib_handler(ftype)
        fd = open(os.path.join(self._folder, log_file), 'rb')
        data = gaussian_handler(fd)
        fd.close()
        return data

    def add_initialgeom(self, file_name):
        
        self._init_data = self._get_data(file_name, 'gausscom')

    def add_optimisation(self, file_name):
        
        if type(file_name) is list:
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

    def get_basis_descript(self):
        
        return self._opt_data.read('basis_descript')

    def get_basis_funcs(self):
        
        return self._opt_data.read('nbasis')
        
    def is_optimised(self):
        
        return self._opt_data.read('optdone')

    def get_optimisation_E(self, units='eV'):
        
        energies = self._opt_data.read('scfenergies')
        
        if not units == 'eV':
            energies = convertor(energies, 'eV', units)
        return energies[-1]            
            
    def plot_optimisation_E(self, units='eV'):
        
        energies = self._opt_data.read('scfenergies')
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

    def is_conformer(self):
        
        imgaginary_freqs = self._freq_data.read('vibfreqs') < 0. 
        return not imgaginary_freqs.any()
        
    def plot_IRfreqs(self):
        
        frequencies = self._freq_data.read('vibfreqs')
        irs = self._freq_data.read('vibirs')

        f, ax = plt.subplots()        
        ax.scatter(frequencies, irs)   
        ax.set_xlabel('Frequency ($cm^{-1}$)')
        ax.set_ylabel('IR Intensities ($km/mol$)')
        ax.grid(True)
        
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
        for coord in coords:
            yield np.dot(transform_matrix, 
                    np.array(np.append(coord, 1))[np.newaxis].T)[:-1].flatten()
    
    def _realign_geometry(self, r_array, align_indxs):
        """inputs coordinate array, index 1, index 2, index 3 """
        
        a1, a2, a3 = align_indxs
        t_matrix = self._create_transform_matrix(r_array[a1], r_array[a2], 
                                                 r_array[a3])
        new_array=np.array(
            [r for r in self._apply_transfom_matrix(t_matrix, r_array)])
        return new_array

    def _create_molecule(self, optimised=True, opt_step=False, scan_step=False, 
                         gbonds=True, data=None):
        if not optimised:
            molecule = self._init_data.read('molecule')            
        else:
            indata = data if data else self._opt_data
            if not type(opt_step) is bool:
                molecule = indata.read('molecule', step=opt_step) 
            elif not type(scan_step) is bool:
                molecule = indata.read('molecule', scan=scan_step)
            else:
                molecule = indata.read('molecule') 
            
        if gbonds: molecule.guess_bonds()
            
        if self._alignment_atom_indxs:
            a, b, c = self._alignment_atom_indxs
            molecule.r_array = self._realign_geometry(molecule.r_array, 
                                                      [a-1, b-1, c-1])

        return molecule
    
    #instead of from chemlab.notebook import display_molecule to add ball_stick
    def _view_molecule(self, molecule, ball_stick=False, colorlist=[]):
        
        topology = {
            'atom_types': molecule.type_array,
            'bonds': molecule.bonds
        }
    
        mv = MolecularViewer(molecule.r_array, topology)
        
        if molecule.n_bonds != 0:
            if ball_stick:
                mv.ball_and_sticks(colorlist=colorlist)
            else:
                mv.points(size=0.15, colorlist=colorlist)
                mv.lines(colorlist=colorlist)
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

    def _image_molecule(self, molecule, ball_stick=False, colorlist=[],
                         rotation=[0., 0., 0.], width=300, height=300, zoom=1.,
                        lines=[], linestyle='impostors'):

        v = QtViewer()
        w = v.widget
        w.camera.orbit_z = MethodType(orbit_z, w.camera)
        
        w.initializeGL()

        if ball_stick:
            r = v.add_renderer(BallAndStickRenderer,
                                molecule.r_array,
                                molecule.type_array,
                                molecule.bonds,
                                rgba_array=colorlist,
                                linestyle=linestyle)
        else:
            r = v.add_renderer(AtomRenderer, 
                                molecule.r_array, 
                                molecule.type_array, 
                                rgba_array=colorlist)
        
        for line in lines:
            #line = [start_coord, end_coord, start_color, end_color, width, dashed]
            #for some reason it didn't like unpacking them to named variables
            v.add_renderer(LineRenderer, [line[0], line[1]], 
                           [[str_to_colour(line[2]), str_to_colour(line[3])]], 
                             width=line[4], dashed=line[5])

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
        
        if len(images) == 1: return images[0]
        
        total_width = sum([img.size[0] for img in images]) + len(images)*gap
        max_height = max([img.size[1] for img in images])
        
        final_img = PIL.Image.new("RGBA", (total_width, max_height), color='white')
        
        horizontal_position = 0
        for img in images:
            final_img.paste(img, (horizontal_position, 0))
            horizontal_position += img.size[0] + gap
        
        return final_img

    def _color_to_transparent(self, image, colour=(255, 255, 255)):
        """ makes colour (default: white) in the image transparent """
        datas = image.getdata()
    
        newData = []
        for item in datas:
            if item[0] == colour[0] and item[1] == colour[1] and item[2] == colour[2]:
                newData.append((colour[0], colour[1], colour[2], 0))
            else:
                newData.append(item)
    
        image.putdata(newData)
        
        return image

    def _show_molecule(self, molecule, active=False, 
                       ball_stick=False, zoom=1., width=300, height=300,
                       rotations=[[0., 0., 0.]],
                       colorlist=[], lines=[], axis_length=0,
                       linestyle='impostors'):
                
        if active:
            return self._view_molecule(molecule, ball_stick=ball_stick, 
                                          colorlist=colorlist)
        else:
            drawlines=lines[:]
            if axis_length:
                drawlines.append([(-1*axis_length,0,0), (axis_length,0,0), 
                              'red', 'dark_red', 3, True])
                drawlines.append([(0,-1*axis_length,0), (0,axis_length,0), 
                              'light_green', 'dark_green', 3, True])
                drawlines.append([(0,0,-1*axis_length), (0,0,axis_length), 
                              'light_blue', 'dark_blue', 3, True])

            images = []
            for rotation in rotations:
                images.append(self._image_molecule(molecule, 
                                    ball_stick=ball_stick, colorlist=colorlist, 
                                    rotation=rotation, zoom=zoom,
                                    width=width, height=width,
                                    lines=drawlines, linestyle=linestyle))  
            image = self._concat_images_horizontal(images)
            del images
                
            b = BytesIO()
            image.save(b, format='png')    
            return ipy_Image(data=b.getvalue())
            
    def show_initial(self, gbonds=True, active=False, ball_stick=False, 
                     rotations=[[0., 0., 0.]], zoom=1., width=300, height=300,
                     axis_length=0, lines=[]):
        
        molecule = self._create_molecule(optimised=False, gbonds=gbonds)
        
        return self._show_molecule(molecule, active=active, 
                                   ball_stick=ball_stick, 
                                   rotations=rotations, zoom=zoom, 
                                   lines=lines, axis_length=axis_length)      
       
    def show_optimisation(self, opt_step=False, gbonds=True, active=False,
                          ball_stick=False, rotations=[[0., 0., 0.]], zoom=1.,
                          width=300, height=300, axis_length=0, lines=[]):
        
        molecule = self._create_molecule(optimised=True, opt_step=opt_step, 
                                         gbonds=gbonds)

        return self._show_molecule(molecule, active=active, 
                                  ball_stick=ball_stick, 
                                  rotations=rotations, zoom=zoom,
                                  lines=lines, axis_length=axis_length,
                                  width=width, height=height)             

    def _rgb_to_hex(self, rgb):
        
        return int('0x%02x%02x%02x' % rgb[:3], 16)
       
    def show_highlight_atoms(self, atomlists, gbonds=True, active=False, 
                             optimised=True,
                        ball_stick=False, rotations=[[0., 0., 0.]], zoom=1.,
                        width=300, height=300, axis_length=0, lines=[]):
               
        if optimised:
            natoms = self._opt_data.read('natom')        
        else:
            natoms = self._init_data.read('natom')

        norm = mpl.colors.Normalize(vmin=1, vmax=len(atomlists))
        cmap = cm.jet_r
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        
        colorlist = [(211, 211, 211, 150) for n in range(natoms)]
        
        for n in range(natoms):
            for group, atomlist in enumerate(atomlists):
                if n+1 in atomlist:
                    colorlist[n] = m.to_rgba(group+1, bytes=True)
                    break
          
        if active:           
            colorlist = [self._rgb_to_hex(col) for col in colorlist]
            
        molecule = self._create_molecule(optimised=optimised, gbonds=gbonds)

        return self._show_molecule(molecule, active=active, 
                                  ball_stick=ball_stick, 
                                  rotations=rotations, zoom=zoom,
                                  colorlist=colorlist,
                                  lines=lines, axis_length=axis_length,
                                  width=width, height=height) 
                                  
    def _get_charge_colors(self, relative=False, minval=-1, maxval=1, alpha=None):
        
        charges = self._nbo_data.read('atomcharges')['natural']
        if relative: minval, maxval = (min(charges), max(charges))
        norm = mpl.colors.Normalize(vmin=minval, vmax=maxval)
        cmap = cm.bwr
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors=m.to_rgba(charges, alpha=alpha, bytes=True)
        
        return colors

    def show_nbo_charges(self, gbonds=True, active=False, 
                         relative=False, minval=-1, maxval=1,
                         ball_stick=False, rotations=[[0., 0., 0.]], zoom=1.,
                         width=300, height=300, axis_length=0, lines=[]):
        
        colorlist = self._get_charge_colors(relative, minval, maxval)

        molecule = self._create_molecule(optimised=True, gbonds=gbonds)

        return self._show_molecule(molecule, active=active, 
                                  ball_stick=ball_stick, 
                                  rotations=rotations, zoom=zoom,
                                  colorlist=colorlist,
                                  lines=lines, axis_length=axis_length,
                                  width=width, height=height) 

    def _converter(self, val, unit1, unit2):
        
        multiple = {('nm', 'nm') : 1.,
                    ('nm', 'Angstrom') : 0.1}

        return val * multiple[(unit1, unit2)]      

    def calc_min_dist(self, idx_list1, idx_list2, optimisation=True, units='nm',
                      ignore_missing=True):
        """ indexes start at 1 """

        if optimisation:
            molecule = self._opt_data.read('molecule')  
        else:
            molecule = self._init_data.read('molecule')
        
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
            molecule = self._opt_data.read('molecule')  
        else:
            molecule = self._init_data.read('molecule')

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
            molecule = self._opt_data.read('molecule')  
        else:
            molecule = self._init_data.read('molecule')

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
            molecule = self._opt_data.read('molecule')  
        else:
            molecule = self._init_data.read('molecule')

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
        charges = self._nbo_data.read('atomcharges')['natural']
        if atoms==[]: return np.sum(charges)
            
        atoms = np.array(atoms) -1 # 1->0 base        
        try:
            subcharges = charges[atoms]
        except IndexError:
            return np.nan
        
        return np.sum(subcharges)
        
    # TODO alignment atoms resetting on error / calculate without creating molecule
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

        alignto = self._alignment_atom_indxs[:]
        self._alignment_atom_indxs = (p1, p2, p3)

        molecule = self._create_molecule()
        charges = self._nbo_data.read('atomcharges')['natural']
        coords = molecule.r_array   
        
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

        self._alignment_atom_indxs = alignto
        
        return r, theta, phi                     

    def calc_SOPT_bonds(self, min_energy=20.):
        """Second Order Perturbation Theory analysis """
        sopt = self._nbo_data.read('sopt')
        sopt_e = [i for i in sopt if i[4]>=min_energy]

        return sopt_e 

    def show_SOPT_bonds(self, min_energy=20., gbonds=True, active=False,
                        ball_stick=True, rotations=[[0., 0., 0.]], zoom=1.,
                        width=300, height=300, axis_length=0, lines=[],
                        relative=False, minval=-1, maxval=1,
                        alpha=0.5):
        """Second Order Perturbation Theory analysis """
        sopt_e = self.calc_SOPT_bonds(min_energy=min_energy)
        molecule = self._create_molecule(gbonds=gbonds)
        
        drawlines = lines[:]
        for dtype, donors, atype, acceptors, energy in sopt_e:
            d_coord = np.mean([molecule.r_array[d-1] for d in donors], axis=0)
            a_coord = np.mean([molecule.r_array[a-1] for a in acceptors], axis=0)
            drawlines.append([d_coord, a_coord, 'blue', 'red', 5, False])
        
        colorlist = self._get_charge_colors(relative, minval, maxval, alpha=alpha)
        
        return self._show_molecule(molecule, active=active, 
                                  ball_stick=ball_stick, 
                                  rotations=rotations, zoom=zoom,
                                  colorlist=colorlist,
                                  lines=drawlines, axis_length=axis_length,
                                  width=width, height=height, linestyle='lines') 

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
        for scan in scan_datas:
            for i in range(scan.read('nscans')):
                mol = scan.read('molecule', scan=i)
                angles.append(func(fixed_atoms, mol=mol))
            energies.extend(scan.read('scanenergies'))
        
        # remove duplicate angles and sort by angle 
        # so that the local max are found correctly
        df = pd.DataFrame({'energy':convertor(np.array(energies), 'eV', eunits), 
                           'angle':angles})
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
        'local_maxs' : argrelextrema(energies, np.greater, mode='wrap', order=order)[0],
        'local_mins' : argrelextrema(energies, np.less, mode='wrap', order=order)[0],
        'global_min' : [np.argmin(energies)],
        'global_max' : [np.argmax(energies)]}
         
        for indx in feature_dict[img_pos]:
            pscans = 0
            for scan in scan_datas:
                if indx < scan.read('nscans') + pscans:
                    mol = self._create_molecule(data=scan, scan_step=indx-pscans)
                    img = self._image_molecule(mol, rotation=rotation, ball_stick=True)
                    break
                else:
                    pscans += scan.read('nscans') 
                
            img = self._color_to_transparent(img)
            self._img_to_plot(angles[indx], energies[indx], img, zoom=zoom, ax=ax)
            
        return ax
