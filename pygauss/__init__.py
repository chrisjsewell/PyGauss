# -*- coding: utf-8 -*-

"""Python Gaussian chemical computation output analysis 

This is a layer on top of the cclib/chemlab/chemview packages for analysing 
Gaussian input/output.

"""
from ._version import __version__

# mokeypatch cclib/chemlab classes improved for this implementation 
# TODO don't think this is working, use mock?

import os, sys, platform

##TODO a weird thing sometimes happens if docx installed, 
#whereby 'import enum' installs docx/enum instead of enum (req for numba)!?
if platform.system() == 'Windows':
    import imp
    init = os.path.abspath(__file__)
    pkg = os.path.dirname(init)
    sites = os.path.dirname(pkg)
    imp.load_package('enum', os.path.join(sites, 'enum'))        
    
from .cclib_patch.parser import data

sys.modules['cclib.parser.data'] = data 

import cclib # version 1.3
import chemlab # version 0.4

from .molecule import Molecule
from .analysis import Analysis
from .file_io import Folder
from .docs import MSDocument
from .utils import df_to_img, set_imagik_exe

import cPickle as pkl
import datetime
def save_object(obj, filename, descript=''):
    filepath = filename + '.pkl'
    with open(filepath, 'wb') as f:
        pkl.dump({'date':datetime.datetime.now(),
                  'pg_version':__version__,
                  'description':descript}, 
                 f)
        pkl.dump(obj, f)
    return os.path.abspath(filepath)
def load_object(filename):
    filepath = filename + '.pkl'
    with open(filepath, 'rb') as f:
        meta = pkl.load(f)
        if not meta['pg_version'] == __version__:
            print 'Warning file ({0}) and package ({1}) versions differs'.format(
                            meta['pg_version'], __version__)
        obj = pkl.load(f)
    return obj
def load_object_meta(filename):
    filepath = filename + '.pkl'
    with open(filepath, 'rb') as f:
        meta = pkl.load(f)
    return meta

import inspect
from . import test_data
def get_test_folder():
    """return a folder obj of test data """
    return Folder(os.path.dirname(os.path.abspath(inspect.getfile(test_data))))


def run_nose(doctests=False, verbose=False):
    import pygauss, nose
    nose_argv = sys.argv
    nose_argv += ['--detailed-errors', '--exe']
    if verbose:
        nose_argv.append('-v')
    if doctests:
        nose_argv.append('--with-doctest')        
    initial_dir = os.getcwd()
    my_package_file = os.path.abspath(pygauss.__file__)
    my_package_dir = os.path.dirname(my_package_file)
    os.chdir(my_package_dir)
    try:
        nose.run(argv=nose_argv)
    finally:
        os.chdir(initial_dir)
