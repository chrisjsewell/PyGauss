# -*- coding: utf-8 -*-
"""
Created on Mon May 25 15:23:49 2015

@author: chris
taken liberally from chemview!
"""
import numpy as np

from chemview.marchingcubes import marching_cubes

##TODO not setting normals correctly (some are negative)
def get_isosurface(coordinates, function, isolevel=0.3, color=(255, 0, 0, 255),
                   resolution=32):
    '''Add an isosurface to the current scene.

    Parameters
    ----------
    coordinates : numpy.array
        the coordinates of the system
    function : func
        A function that takes x, y, z coordinates as input and is broadcastable 
        using numpy. Typically simple functions that involve standard 
        arithmetic operations and functions such as ``x**2 + y**2 + z**2`` or 
        ``np.exp(x**2 + y**2 + z**2)`` will work. If not sure, you can first
        pass the function through ``numpy.vectorize``.
        Example: ``mv.add_isosurface(np.vectorize(f))``
    isolevel : float
        The value for which the function should be constant.
    color : (int, int, int, int)
        The color given in RGBA format
    resolution : int
        The number of grid point to use for the surface. An high value will 
        give better quality but lower performance.

    '''
    if coordinates.shape[0] == 0:
        return 

    # We want to make a container that contains the whole molecule
    # and surface
    coordinates = coordinates.astype('float32')
    area_min = coordinates.min(axis=0) - 0.2
    area_max = coordinates.max(axis=0) + 0.2

    x = np.linspace(area_min[0], area_max[0], resolution)
    y = np.linspace(area_min[1], area_max[1], resolution)
    z = np.linspace(area_min[2], area_max[2], resolution)

    xv, yv, zv = np.meshgrid(x, y, z)
    spacing = np.array((area_max - area_min)/resolution)

    if isolevel >= 0:
        triangles = marching_cubes(function(xv, yv, zv), isolevel)
    else: # Wrong triangle unwinding order -- god only knows why
        triangles = marching_cubes(-function(xv, yv, zv), -isolevel)

    if len(triangles) == 0:
        ## NO surface
        return 

    verts = []
    for i, t in enumerate(triangles):
       verts.extend(t)

    verts = area_min + spacing/2 + np.array(verts)*spacing
    
    # "Normals, this is quite easy since they are the vertices"??
    normals = []
    for vertex in verts:
        normals.append(vertex/np.linalg.norm(vertex))
    normals = np.array(normals) # Numpyize  

    colors = np.array([color for _ in range(verts.shape[0])]) 

    return verts.astype('float32'), normals.astype('float32'), colors

if __name__ == '__main__':

    from chemlab.io import remotefile
    url = "https://raw.githubusercontent.com/cclib/cclib/master/data/GAMESS/basicGAMESS-US2012/water_mp2.out"
    df = remotefile(url, "gamess")
    
    molecule = df.read('molecule')
    molecule.guess_bonds()
    
    mocoeffs = df.read('mocoeffs')
    gbasis = df.read('gbasis')

    from chemlab.qc import molecular_orbital
    orbital = 2
    iso_level = -0.3
    color = (255, 0, 0, 100)
    coefficients = mocoeffs[0][orbital]
    f = molecular_orbital(molecule.r_array, coefficients, gbasis)
    verts, normals, colors = get_isosurface(molecule.r_array, f, 
                                            iso_level, color)


