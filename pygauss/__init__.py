# -*- coding: utf-8 -*-

"""Python Gaussian chemical computation output analysis 

This is a layer on top of the cclib/chemlab/chemview packages for analysing 
Gaussian input/output.

"""
from ._version import __version__

# mokeypatch cclib/chemlab classes improved for this implementation 
# TODO don't think this is working, use mock?

import sys

from .cclib_patch.parser import data

sys.modules['cclib.parser.data'] = data 

import cclib # version 1.3
import chemlab # version 0.4

from .molecule import Molecule
from .analysis import Analysis
from .file_io import Folder

import os, inspect
from . import test_data
def get_test_folder():
    """return a folder obj of test data """
    return Folder(os.path.dirname(os.path.abspath(inspect.getfile(test_data))))

from .utils import df_to_img, set_imagik_exe

def run_nose(verbose=False):
    import pygauss, nose
    nose_argv = sys.argv
    nose_argv += ['--detailed-errors', '--exe']
    if verbose:
        nose_argv.append('-v')
    initial_dir = os.getcwd()
    my_package_file = os.path.abspath(pygauss.__file__)
    my_package_dir = os.path.dirname(my_package_file)
    os.chdir(my_package_dir)
    try:
        nose.run(argv=nose_argv)
    finally:
        os.chdir(initial_dir)
