# -*- coding: utf-8 -*-

"""PYthon GAUSSian DFT output analysis 

This is a layer on top of the cclib/chemlab/chemview packages for analysing 
Gaussian output and outputting it to a word document.
S
"""
#import cclib # version 1.3
#import chemlab # version 0.4

# mokeypatch cclib/chemlab classes improved for this implementation 
# TODO don't think this is working, use mock?

import sys

#from .chemlab_patch.io.handlers import _cclib
#sys.modules['chemlab.io.handlers._cclib'] = _cclib

from .cclib_patch.parser import data#, gaussianparser

sys.modules['cclib.parser.data'] = data 
#sys.modules['cclib.parser.gaussianparser'] = gaussianparser

#from .chemlab_patch.graphics import camera

#sys.modules['chemlab.graphics.camera']= camera

#from .chemlab_patch.graphics.renderers import atom, ballandstick, bond, line

#sys.modules['chemlab.graphics.renderers.atom'] = atom
#sys.modules['chemlab.graphics.renderers.ballandstick'] = ballandstick
#sys.modules['chemlab.graphics.renderers.bond'] = bond
#sys.modules['chemlab.graphics.renderers.line'] = line

from . import analysis
from . import molecule

import os, inspect
from . import test_data
def get_test_folder():
    return os.path.dirname(inspect.getfile(test_data))
