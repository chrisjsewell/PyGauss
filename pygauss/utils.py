# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 01:08:30 2015

@author: chris
"""
from numpy import bmat, hstack, dot, ones, zeros, sum, asarray
from numpy.linalg import solve

from pandas.core.index import MultiIndex

def is_wellcentered(pts, tol=1e-8):
    """
    Determines whether the M points in N dimensions define a 
    well-centered simplex.
    """
    bary_coords = circumcenter_barycoords(pts)    
    return min(bary_coords) > tol

def circumcenter_barycoords(pts):
    """
    Computes the barycentric coordinates of the circumcenter M, N-dimensional
    points (1 <= M <= N + 1 and N >= 1). The points are given by the rows of 
    an (M)x(N) dimensional matrix pts.
    
    Uses an extension of the method described here:
    http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
    """    

    pts = asarray(pts)

    rows,cols = pts.shape

    assert(rows <= cols + 1)    

    A = bmat( [[ 2*dot(pts,pts.T), ones((rows,1)) ],
               [  ones((1,rows)) ,  zeros((1,1))  ]] )

    b = hstack((sum(pts * pts, axis=1),ones((1))))
    x = solve(A,b)
    bary_coords = x[:-1]  

    return bary_coords
    
def circumcenter(pts):
    """
    Computes the circumcenter and circumradius of M, N-dimensional
    points (1 <= M <= N + 1 and N >= 1). The points are given by the rows of 
    an (M)x(N) dimensional maatrix pts.  
    
    Returns a tuple (center, radius) where center is a
    column vector of length N and radius is a scalar.
        
        In the case of four points in 3D, pts is a 4x3 matrix arranged as:
            
        pts = [ x0 y0 z0 ]
              [ x1 y1 z1 ]
              [ x2 y2 z2 ]          
              [ x3 y3 z3 ]
        
        with return value ([ cx cy cz ], R)    
    
    Uses an extension of the method described here:
    http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
    """
    pts = asarray(pts)      
    bary_coords = circumcenter_barycoords(pts)
    center = dot(bary_coords,pts)
    #radius = norm(pts[0,:] - center)
    return center

import os    
import subprocess
from subprocess import PIPE
try:
    from subprocess import DEVNULL # Python 3
except ImportError:
    DEVNULL = open(os.devnull, 'r+b', 0)

from IPython.display import Image

def df_to_img(df, na_rep='-', im_exe='convert', other_temp=None):
    """ converts a pandas Dataframe to an IPython image 
    
    na_rep : str
        how to represent empty (nan) cells
    im_name : str
        the name of the imagemagick executable
    other_temp : str
        a latex template to use for the table other than the default 
        
    The function uses pandas to convert the dataframe to a latex table, 
    applies a template, converts to a PDF, converts to an image,
    and finally return the image
    
    to use this function you will need 
    the pdflatex executable from tex distribution,
    the convert executable from imagemagick, which also requires ghostscript; 
    http://www.ghostscript.com/download/gsdnld.html
    http://www.imagemagick.org/script/binary-releases.php    
    
    NB: on Windows some issues were found with convert being an already 
    existing application. To overcome this change its filename and use the 
    im_name variable.
    """
    # pandas 0.16 has a bug when using heirarchical row indexes
    use_indx = True
    if type(df.index) == MultiIndex:
        df = df.reset_index()
        use_indx = False
    latex_str = df.to_latex(index=use_indx, escape=False, na_rep=na_rep)
    
    filename = 'df_to_pdf_out.tex'
    pdffile = 'df_to_pdf_out.pdf'
    logname = 'df_to_pdf_out.log'
    auxname = 'df_to_pdf_out.aux'
    imgname = 'df_to_pdf_out.png'
    
    template = r'''\documentclass{{article}}
                \usepackage{{graphicx}}
                \usepackage{{booktabs}}
                \pagenumbering{{gobble}}
                \begin{{document}}
                \begin{{table}}[ht]
                \centering
                \resizebox{{\textwidth}}{{!}}
                {}
                \end{{table}}
                \end{{document}}
                '''
    if other_temp: template = other_temp
    
    with open(filename, 'wb') as f:
        f.write(template.format('{'+latex_str+'}'))
    
    try:    
        proc = subprocess.Popen(['pdflatex', filename], 
                        stdin=DEVNULL, stdout=PIPE, stderr=PIPE)
        out, err = proc.communicate()
    except:
        os.unlink(filename)
        if os.path.exists(logname): os.unlink(logname)
        if os.path.exists(auxname): os.unlink(auxname)
        raise RuntimeError('pdflatex not installed')
    
    os.unlink(filename)
    if os.path.exists(logname): os.unlink(logname)
    if os.path.exists(auxname): os.unlink(auxname)
        
    if err: 
        raise RuntimeError('error in pdflatex run:\n {0}'.format(err))
    if not os.path.exists(pdffile): 
        raise RuntimeError('pdflatex did not produce a pdf file')
        
    try:
        proc = subprocess.Popen([im_exe, '-density', '300', '-trim', 
                               pdffile, 
                               '-quality', '100', '-sharpen', '0x1.0', 
                               imgname],
                     stdin=DEVNULL, stdout=PIPE, stderr=PIPE)
        out, err = proc.communicate()
    except:
        os.unlink(pdffile)
        raise RuntimeError('imagemagick (convert) not installed')

    os.unlink(pdffile)
    if err: 
        raise RuntimeError('error in imagemagick run:\n {0}'.format(err))
    if not os.path.exists(imgname): 
        raise RuntimeError('imagemagick did not produce a png file')
    
    ipy_img = Image(filename=imgname)
    os.unlink(imgname)
    
    return ipy_img


if __name__ == '__main__':    
    print circumcenter([[1, 0, 0], [0, 1, 0], [0, 0, 1]])