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

import os, inspect
from . import test_data
def get_test_folder():
    return os.path.dirname(os.path.abspath(inspect.getfile(test_data)))

from .utils import df_to_img, ipy_img_tofile, set_imagik_exe

from . import test_molecule, test_analysis
