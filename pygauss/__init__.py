# -*- coding: utf-8 -*-

""" 

This is a layer on top of the cclib/chemlab/chemview packages for analysing 
Gaussian output and outputting it to a word document.

"""


import chemlab # version 0.4
import cclib # version 1.3

# mokeypatch classes improved for this implementation 

from .cclib_patch.parser import data, gaussianparser

cclib.parser.data = data 
cclib.parser.gaussianparser = gaussianparser

from .chemlab_patch.io.handlers import _cclib

chemlab.io.handlers._cclib = _cclib

from .chemlab_patch.graphics import camera

chemlab.graphics.camera = camera

from .chemlab_patch.graphics.renderers import atom, ballandstick, bond

chemlab.graphics.renderers.atom = atom
chemlab.graphics.renderers.ballandstick = ballandstick
chemlab.graphics.renderers.bond = bond

from . import analysis