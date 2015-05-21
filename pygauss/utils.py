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
# On OSx Conda creates its own environment with a reduced $PATH variable 
if os.path.exists(os.path.sep + os.path.join('usr', 'local', 'bin')):
    os.environ["PATH"] += os.pathsep + os.path.sep + os.path.join('usr', 'local', 'bin')
if os.path.exists(os.path.sep + os.path.join('usr', 'texbin')):
    os.environ["PATH"] += os.pathsep + os.path.sep + os.path.join('usr', 'texbin')

import subprocess
from subprocess import PIPE
try:
    from subprocess import DEVNULL # Python 3
except ImportError:
    DEVNULL = open(os.devnull, 'r+b', 0)

from IPython.display import Image as IPy_Image
#from PIL import Image
#from io import BytesIO
import random
import warnings

from scipy.interpolate import interp1d
_IMGMAGIK = 'convert'
def set_imagik_exe(name):
    assert type(name) is str
    global _IMGMAGIK
    _IMGMAGIK=name

_PTS_TO_PIX = {
10.5: 14,
  10: 13,
  11: 15,
  12: 16,
  13.5: 18,
  13: 17,
  14.5: 20,
  14: 19,
  15: 21,
  16: 22,
  17: 23,
  18: 24,
  20: 26,
  22: 29,
  24: 32,
  26: 35,
  27: 36,
  28: 37,
  29: 38,
  30: 40,
  32: 42,
  34: 45,
  36: 48,
  6: 8,
  7.5: 10,
  7: 9,
  8: 11,
  9: 12}
  
def df_to_img(df, na_rep='-', other_temp=None, font_size=None,
               width=None, height=None, unconfined=False):
    """ converts a pandas Dataframe to an IPython image 
    
    na_rep : str
        how to represent empty (nan) cells
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
    if font_size:
        #estimate height of table to give certain font size
        row_height = interp1d(_PTS_TO_PIX.keys(),
                              _PTS_TO_PIX.values())(font_size) 
        if hasattr(df.columns, 'levels'):
            header_rows = len(df.columns.levels)
        else:
            header_rows = 1
        tbl_height = int((header_rows + df.shape[0])*row_height*1.15)

    # pandas 0.16 has a bug when using heirarchical row indexes
    use_indx = True
    if type(df.index) == MultiIndex:
        df = df.reset_index()
        use_indx = False
    latex_str = df.to_latex(index=use_indx, escape=False, na_rep=na_rep)
    
    rand = random.randint(1, 100000)
    filename = 'df_to_pdf_out{0}.tex'.format(rand)
    pdffile = 'df_to_pdf_out{0}.pdf'.format(rand)
    logname = 'df_to_pdf_out{0}.log'.format(rand)
    auxname = 'df_to_pdf_out{0}.aux'.format(rand)
    imgname = 'df_to_pdf_out{0}.png'.format(rand)
    
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
        
    inargs = [_IMGMAGIK, '-trim', '-density', '600',  pdffile, 
              '-quality', '100', '-sharpen', '0x1.0', imgname]
    if font_size:
        inargs.insert(2, '-resize')
        inargs.insert(3, 'x{0}'.format(tbl_height))
    try:
        proc = subprocess.Popen(inargs,
                     stdin=DEVNULL, stdout=PIPE, stderr=PIPE)
        out, err = proc.communicate()
    except:
        os.unlink(pdffile)
        raise RuntimeError('imagemagick (convert) not installed')

    os.unlink(pdffile)
    if err: 
        if not os.path.exists(imgname):
            raise RuntimeError('error in imagemagick run:\n {0}'.format(err))
        else:
            warnings.warn('non-fatal error in imagemagick run:\n {0}'.format(err))
    if not os.path.exists(imgname): 
        raise RuntimeError('imagemagick did not produce a png file')
    
    
#    img = Image.open(imgname)
#    hpercent = (tbl_height / float(img.size[1]))
#    tbl_width = int((float(img.size[0]) * float(hpercent)))
#
#    img = img.resize((tbl_width, tbl_height), Image.ANTIALIAS)    
#
#    b = BytesIO()
#    img.save(b, format='png')    
#    
#    ipy_img = IPy_Image(data=b.getvalue(), width=width, height=height, 
#                    unconfined=unconfined)
    ipy_img = IPy_Image(filename=imgname, width=width, height=height, 
                    unconfined=unconfined)
    os.unlink(imgname)
    
    return ipy_img

import matplotlib.pyplot as plt
def imgplot_kmean_groups(analysis, category, cat_name, groups,
                columns, filters={}, output=[], max_cols=2, **kwargs):
    df = analysis.calc_kmean_groups(
        category, cat_name, groups, columns=columns, 
        filters=filters)
    for cat, gf in df.groupby('Category'):
        fig, caption = analysis.plot_mol_images(
            max_cols=max_cols, info_columns=output,
            rows=gf.index.tolist(), **kwargs)
        fig.suptitle('Category {0}:'.format(cat+1),
                    fontsize=20, x=0, color='blue')
        fig.subplots_adjust(top=0.8)        
        plt.show()  
        print caption
        
if __name__ == '__main__':    
    print circumcenter([[1, 0, 0], [0, 1, 0], [0, 0, 1]])