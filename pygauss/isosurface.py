# -*- coding: utf-8 -*-
"""
Created on Mon May 25 15:23:49 2015

@author: chris
based on add_isosurface function from chemview
"""
import numpy as np

#TODO not making this a depedency until it works
try:
    from skimage.measure import marching_cubes, correct_mesh_orientation
except ImportError:
    pass

def get_isosurface(coordinates, function, isolevel=0.03, color=(255, 0, 0, 255),
                   extents=(5,5,5), resolution=100):
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
    extents : (float, float, float)
        +/- x,y,z to extend the molecule geometrty when constructing the surface
    resolution : int
        The number of grid point to use for the surface. An high value will 
        give better quality but lower performance.

    '''
    if coordinates.shape[0] == 0:
        return False

    # We want to make a container that contains the whole molecule
    # and surface
    coordinates = coordinates.astype('float32')
    area_min = coordinates.min(axis=0) - np.array(extents)
    area_max = coordinates.max(axis=0) + np.array(extents)

    x = np.linspace(area_min[0], area_max[0], resolution)
    y = np.linspace(area_min[1], area_max[1], resolution)
    z = np.linspace(area_min[2], area_max[2], resolution)

    xv, yv, zv = np.meshgrid(x, y, z)
    spacing = np.array((area_max - area_min)/resolution)

    if isolevel >= 0:
        sign = 1
    else:
        sign = -1
        
    volume = sign*function(xv, yv, zv)
    if volume.flatten().min() <= sign*isolevel and volume.flatten().max() >= sign*isolevel: 
        verts, faces = marching_cubes(volume, sign*isolevel) 
        #don't know if i need this                       
        faces = correct_mesh_orientation(volume, verts, faces, 
                                         gradient_direction='descent')
        verts = area_min + spacing/2 + verts*spacing

        # TODO for some reason need to invert, but no one knows what the problem is
        verts[:,[0,1]] = verts[:,[1,0]]                                
        # TODO and its messing up the normals so using double facing 
        #(kind of works if transparent)
        opposite_faces = faces.copy()
        opposite_faces[:,[0,1]] = opposite_faces[:,[1,0]]     
        faces = np.concatenate([faces,opposite_faces])                           
        
    else:
        return

    normals = calc_normals(verts, faces) 

    verts = verts[faces.flatten()]
    normals = normals[faces.flatten()]
    
    colors = np.array([color for _ in range(verts.shape[0])]) 
    
    return verts.astype('float32'), normals.astype('float32'), colors

def normalize_v3(arr):
    ''' Normalize a numpy array of 3 component vectors shape=(n,3) '''
    lens = np.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens                
    return arr

def calc_normals(verts, faces):
    """from; 
    https://sites.google.com/site/dlampetest/python/calculating-normals-of-a-triangle-mesh-using-numpy
    """
    #Create a zeroed array with the same type and shape as our vertices i.e., per vertex normal
    norm = np.zeros( verts.shape, dtype=verts.dtype )
    #Create an indexed view into the vertex array using the array of three indices for triangles
    tris = verts[faces]
    #Calculate the normal for all the triangles, by taking the cross product of the vectors v1-v0, and v2-v0 in each triangle             
    n = np.cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )
    # n is now an array of normals per triangle. The length of each normal is dependent the vertices, 
    # we need to normalize these, so that our next step weights each normal equally.
    normalize_v3(n)
    # now we have a normalized array of normals, one per triangle, i.e., per triangle normals.
    # But instead of one per triangle (i.e., flat shading), we add to each vertex in that triangle, 
    # the triangles' normal. Multiple triangles would then contribute to every vertex, so we need to normalize again afterwards.
    # The cool part, we can actually add the normals through an indexed view of our (zeroed) per vertex normal array
    norm[ faces[:,0] ] += n
    norm[ faces[:,1] ] += n
    norm[ faces[:,2] ] += n
    normalize_v3(norm)
    
    return norm

def my_calc_normals(verts, faces):
    """ doesn't work """
    normals=[]
    #for each vertex
    for i in range(verts.shape[0]):
        #find all triangles with vertex in
        mask = faces==i
        tris = faces[mask.any(axis=1)]
        # find the normal of each triangle
        vs = tris[tris!=i].reshape((tris.shape[0],2))
        v0 = verts[i]
        v1s = verts[vs[:,0]]
        v2s = verts[vs[:,1]]
        vi = v1s - v0
        vj = v2s - v0
        norms = np.cross(vi, vj)        
        scalar = 1/np.linalg.norm(norms, axis=1)
        norms = norms * scalar.reshape((scalar.shape[0],1))
        #sum all the normals and normalise        
        norm = np.sum(norms, axis=0)
        normals.append(norm / np.linalg.norm(norm))
    return np.array(normals) # Numpyize  
    
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


