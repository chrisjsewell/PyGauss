'''Hidden module to abstract the cclib interface and convert it to chemlab

'''

# CJS changed relative paths to cclib ones

from chemlab.io.handlers.base import FormatNotSupported, IOHandler
from chemlab.core import Molecule
from chemlab.db import ChemlabDB
import numpy as np

import logging

from cclib.parser import (ADF,
                          GAMESS, 
                          GAMESSUK, Gaussian, Jaguar, Molpro, 
                          #NWChem, 
                          ORCA, 
                          # Psi, 
                          #QChem
                          )
from ....cclib_patch.parser import gausscomparser
from ....cclib_patch.parser import gaussianparser

_types = { 'gamess' : GAMESS,
           'gamessuk': GAMESSUK, 
           'gaussian':gaussianparser.Gaussian, 
           'gausscom': gausscomparser.Gausscom,
           'jaguar': Jaguar, 
           'molpro': Molpro, 
           #'nwchem': NWChem, 
           'orca': ORCA
         }

cdb = ChemlabDB()
symbols = cdb.get('data', 'symbols')

def _create_cclib_handler(type):
    if type not in _types:
        raise FormatNotSupported(type)

    class _Handler(IOHandler):

        def __init__(self, fd, filetype=type):
            super(_Handler, self).__init__(fd)
            self.filetype = filetype
            self.data = _types[filetype](fd, loglevel=logging.ERROR).parse()

        def read(self, feature, *args, **kwargs):
            if feature == 'molecular orbitals':
                n = args[0]
                basis_functions = self.data.gbasis
            if feature == 'molecule':
                # Angstrom to nanometers
                # CJS added option to get molecule from steps of optimisation
                # or optimisations for PES scans
                # TODO error checking
                if kwargs.has_key('step'):
                    return Molecule.from_arrays(r_array=self.data.atomcoords[kwargs['step']]/10, 
                                        type_array=np.array([symbols[a] for a in self.data.atomnos]))
                elif kwargs.has_key('scan'):
                    return Molecule.from_arrays(r_array=self.data.scancoords[kwargs['scan']]/10, 
                                        type_array=np.array([symbols[a] for a in self.data.atomnos]))
                else:
                    return Molecule.from_arrays(r_array=self.data.atomcoords[-1]/10, 
                                        type_array=np.array([symbols[a] for a in self.data.atomnos]))
            else:
                return getattr(self.data, feature)

        def available_features(self):
            return set(self.data._attrlist) & set(dir(self.data))

    # Setting consistent class Name
    _Handler.__name__ = _types[type].__name__ + 'Handler'

    return _Handler

_cclib_handlers = [(_create_cclib_handler(format), format) for format in _types]
