# -*- coding: utf-8 -*-
import os
from io import BytesIO
import PIL

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from IPython.display import Image as ipy_Image

from chemlab.io.handlers._cclib import _create_cclib_handler
from chemlab.graphics.qtviewer import QtViewer
from chemlab.graphics.renderers import AtomRenderer, BallAndStickRenderer
from chemlab.graphics.postprocessing import SSAOEffect # Screen space ambient Occlusion
from chemlab.utils import cartesian
from cclib.parser.utils import convertor

#instead of chemview MolecularViewer to add defined colouring
from .chemview_patch.viewer import MolecularViewer

class Molecule:
    def __init__(self, folder, opt_fname=False, freq_fname=False, alignto=None):
        
        self._folder = folder
        self._opt_data = None
        self._freq_data = None
        self._alignment_atom_indxs = alignto
        
        if opt_fname: self.add_optimisation(opt_fname)
        if freq_fname: self.add_frequency(freq_fname)
    
    def _get_data(self, log_file):
        
        gaussian_handler = _create_cclib_handler('gaussian')
        fd = open(os.path.join(self._folder, log_file), 'rb')
        data = gaussian_handler(fd)
        fd.close()
        return data

    def add_optimisation(self, file_name):
        
        self._opt_data = self._get_data(file_name)

    def add_frequency(self, file_name):
        
        self._freq_data = self._get_data(file_name)
    
    def is_optimised(self):
        
        return self._opt_data.read('optdone')

    def plot_optimisation_E(self, units='eV'):
        
        energies = self._opt_data.read('scfenergies')
        if not units == 'eV':
            energies = convertor(energies, 'eV', units)
        
        plt.plot(energies)   
        plt.ylabel('Energy ({0})'.format(units))
        plt.xlabel('Optimisation Step')
        plt.grid(True)

    def is_conformer(self):
        
        imgaginary_freqs = self._freq_data.read('vibfreqs') < 0. 
        return not imgaginary_freqs.any()
        
    def plot_IRfreqs(self):
        
        frequencies = self._freq_data.read('vibfreqs')
        irs = self._freq_data.read('vibirs')
        
        plt.scatter(frequencies, irs)   
        plt.xlabel('Frequency ($cm^{-1}$)')
        plt.ylabel('IR Intensities ($km/mol$)')
        plt.grid(True)

    def add_alignment_atoms(self, idx1, idx2, idx3):
        
        self._alignment_atom_indxs = (idx1, idx2, idx3)    
        
    def remove_alignment_atoms(self):
        
        self._alignment_atom_indxs = None   

    def _midpoint_atoms(self, molecule, atom_ids):

        return np.mean(molecule.r_array[atom_ids], axis=0)                                 

    def _create_transform_matrix(self, c1, c2, c3):
        """
        A function to take three coordinates and creates a transformation matrix
        that aligns their plane with the standard axes 
        """
        #translate c1 to the origin [0,0,0]
        v2=c2-c1; v3=c3-c1
    
        #now find the orthonormal basis set   
        # a plane is a*x+b*y+c*z+d=0 where[a,b,c] is the normal and d is 0  
        # (since the origin now intercepts the plane). Thus, we calculate;
        normal = np.cross(v2,v3)
        a, b, c = normal
        vf3 = normal/np.linalg.norm(normal)
    
        if c==0.:
            #just in case
            vf1=np.array([0,0.0, 1.])
        else:
            #find initial vector for (x=1, y=0)
            vf1=np.array([1.0,0.0, (-a*1.0-b*0.0)/c])
            vf1 = vf1/np.linalg.norm(vf1)
             
        vf2 = np.cross(vf3, vf1)
        vf2 = vf2/np.linalg.norm(vf2)
    
        #create the translation matrix that moves the new basis to the origin
        ident=np.vstack((np.identity(3), np.zeros(3)))
        translate_matrix = np.hstack((ident, np.array(np.append(-c1, 1))[np.newaxis].T))
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
    
    def _realign_geometry(self, molecule, atom_indxs):
        """inputs molecule, atom index 1, atom index 2, atom index 3 """
        
        a1, a2, a3 = atom_indxs
        t_matrix = self._create_transform_matrix(molecule.r_array[a1-1], 
                                           molecule.r_array[a2-1], 
                                           molecule.r_array[a3-1])
        molecule.r_array=np.array(
            [r for r in self._apply_transfom_matrix(t_matrix, molecule.r_array)])
        return molecule

    def _create_molecule(self, optimised=True, gbonds=True):

        if optimised:
            molecule = self._opt_data.read('molecule')  
        else:
            molecule = self._opt_data.read('molecule', step=0) 
            
        if gbonds: molecule.guess_bonds()
            
        if self._alignment_atom_indxs:
            molecule = self._realign_geometry(molecule, self._alignment_atom_indxs)

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

    def _image_molecule(self, molecule, ball_stick=False, colorlist=[],
                         rotation=[0., 0., 0.], width=300, height=300, zoom=1.):

        v = QtViewer()
        w = v.widget
        w.initializeGL()

        if ball_stick:
            v.add_renderer(BallAndStickRenderer,
                                molecule.r_array,
                                molecule.type_array,
                                molecule.bonds,
                                rgba_array=colorlist)
        else:
            v.add_renderer(AtomRenderer, 
                                molecule.r_array, 
                                molecule.type_array, 
                                rgba_array=colorlist)
        
        v.add_post_processing(SSAOEffect)
        w.camera.orbit_x(rotation[0]*np.pi/180.)
        w.camera.orbit_y(rotation[1]*np.pi/180.)
        w.camera.orbit_z(rotation[2]*np.pi/180.)
        w.camera.autozoom(molecule.r_array*1./zoom)
        
        image = w.toimage(width, height)

        # Cleanup
        del v
        del w
        
        return image

    def _concat_images_horizontal(self, images):
        
        if len(images) == 1: return images[0]
        
        total_width = sum([img.size[0] for img in images])
        max_height = max([img.size[1] for img in images])
        
        final_img = PIL.Image.new("RGB", (total_width, max_height))
        
        horizontal_position = 0
        for img in images:
            final_img.paste(img, (horizontal_position, 0))
            horizontal_position += img.size[0]
        
        return final_img

    def _color_to_transparent(image, colour = (255, 255, 255)):
        """ makes colour in the image transparent """
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
                       ball_stick=False, zoom=1., rotations=[[0., 0., 0.]],
                       colorlist=[]):
                
        if active:
            return self._view_molecule(molecule, ball_stick=ball_stick, 
                                          colorlist=colorlist)
        else:
            images = []
            for rotation in rotations:
                images.append(self._image_molecule(molecule, 
                                    ball_stick=ball_stick, colorlist=colorlist, 
                                    rotation=rotation, zoom=zoom))  
            image = self._concat_images_horizontal(images)
                
            b = BytesIO()
            image.save(b, format='png')    
            return ipy_Image(data=b.getvalue())
            
    def show_initial(self, gbonds=True, active=False, ball_stick=True, 
                     rotations=[[0., 0., 0.]], zoom=1.):
        
        molecule = self._create_molecule(optimised=False, gbonds=gbonds)

        return self._show_molecule(molecule, active=active, 
                                   ball_stick=ball_stick, 
                                   rotations=rotations, zoom=zoom)      
       
    def show_optimisation(self, gbonds=True, active=False, 
                          ball_stick=True, rotations=[[0., 0., 0.]], zoom=1.):
        
        molecule = self._create_molecule(optimised=True, gbonds=gbonds)

        return self._show_molecule(molecule, active=active, 
                                  ball_stick=ball_stick, 
                                  rotations=rotations, zoom=zoom)             

    def _rgb_to_hex(self, rgb):
        
        return int('0x%02x%02x%02x' % rgb[:3], 16)

    def show_highlight_atoms(self, atomlists, gbonds=True, active=False, 
                        ball_stick=True, rotations=[[0., 0., 0.]], zoom=1.):
               
        natoms = self._opt_data.read('natom')        

        norm = mpl.colors.Normalize(vmin=1, vmax=len(atomlists))
        cmap = cm.jet
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        
        colorlist = [(211, 211, 211, 150) for n in range(natoms)]
        
        for n in range(natoms):
            for group, atomlist in enumerate(atomlists):
                if n+1 in atomlist:
                    colorlist[n] = m.to_rgba(group+1, bytes=True)
                    break
          
        if active:           
            colorlist = [self._rgb_to_hex(col) for col in colorlist]
            
        molecule = self._create_molecule(optimised=True, gbonds=gbonds)

        return self._show_molecule(molecule, active=active, 
                                  ball_stick=ball_stick, 
                                  rotations=rotations, zoom=zoom,
                                  colorlist=colorlist) 
                                  
    def _converter(self, val, unit1, unit2):
        
        multiple = {('nm', 'nm') : 1.,
                    ('nm', 'Angstrom') : 0.1}

        return val * multiple[(unit1, unit2)]      

    def calc_min_dist(self, idx_list1, idx_list2, optimisation=True, units='nm'):
        """ indexes start at 1 """

        if optimisation:
            molecule = self._opt_data.read('molecule')  
        else:
            molecule = self._opt_data.read('molecule', step=0)

        indx_combis = cartesian([idx_list1, idx_list2])
        c1 = molecule.r_array[indx_combis[:, 0]-1]
        c2 = molecule.r_array[indx_combis[:, 1]-1]

        dist =  np.min(np.linalg.norm(c1-c2, axis=1))       
        
        return (self._converter(dist, 'nm', units), units)
                                  
    def calc_bond_angle(self, indxs, optimisation=True):
        """ Returns the angle in degrees between three points    """

        if optimisation:
            molecule = self._opt_data.read('molecule')  
        else:
            molecule = self._opt_data.read('molecule', step=0)

        v1 = molecule.r_array[indxs[0]-1] - molecule.r_array[indxs[1]-1]
        v2 = molecule.r_array[indxs[2]-1] - molecule.r_array[indxs[1]-1]
        cosang = np.dot(v1, v2)
        sinang = np.linalg.norm(np.cross(v1, v2))
        
        return np.degrees(np.arctan2(sinang, cosang))

    def calc_dihedral_angle(self, indxs, optimisation=True):
        """ Returns the angle in degrees between four points  """

        if optimisation:
            molecule = self._opt_data.read('molecule')  
        else:
            molecule = self._opt_data.read('molecule', step=0)

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
                                  
    def calc_angle_to_normal(self):

        return                             
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  
                                  