# This file is part of cclib (http://cclib.github.io), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006-2014, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.

"""Parser for Gaussian output files"""


from __future__ import print_function
import re

import numpy

# CJS changed . to cclib.parser
from cclib.parser import logfileparser
from cclib.parser import utils

from chemlab.db import ChemlabDB
cdb = ChemlabDB()
symbols = cdb.get('data', 'symbols')


class Gausscom(logfileparser.Logfile):
    """A Gaussian 98/03 com file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(Gausscom, self).__init__(logname="Gausscom", *args, **kwargs)
        
    def __str__(self):
        """Return a string representation of the object."""
        return "Gaussian com file %s" % (self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'Gausscom("%s")' % (self.filename)
    

    def before_parsing(self):

        self.found_options = False
        self.found_comment = False
        self.found_charge = False
        #self.found_geometry = False
        self.found_e_correl = False
        return

    def after_parsing(self):

        return        
            
    def extract(self, inputfile, line):             
        """Extract information from the file object inputfile."""

        if line[0] == '#' and not self.found_options:
            
            self.found_options = True
        
        elif line.strip() == '' and self.found_options:
            
            if not self.found_comment:
                self.found_comment = True
            elif not self.found_charge:
                self.found_charge = True
                
                line = next(inputfile)
                self.set_attribute('charge', int(line.split()[0]))
                self.set_attribute('mult', int(line.split()[1]))
                
                if not hasattr(self, "atomcoords"):
                    self.atomcoords = []
                self.inputatoms = []
                atomnames=[]
                atomcoords = []

                line = next(inputfile)
                while line.strip() != '':
                    broken = line.split()
                    atomnames.append(symbols.index(broken[0]))
                    atomcoords.append(list(map(float, broken[1:4])))
                    line = next(inputfile)
    
                self.atomcoords.append(atomcoords)
    
                self.inputatoms = atomnames[:]
                self.set_attribute('atomnos', self.inputatoms)
                self.set_attribute('natom', len(atomnames))
                
                    
            else:
                self.found_e_correl = True

if __name__ == "__main__":
    import doctest, gausscomparser, sys

    if len(sys.argv) == 1:
        doctest.testmod(gausscomparser, verbose=False)

    if len(sys.argv) >= 2:
        parser = gausscomparser.Gausscom(sys.argv[1])
        data = parser.parse()

    if len(sys.argv) > 2:
        for i in range(len(sys.argv[2:])):
            if hasattr(data, sys.argv[2 + i]):
                print(getattr(data, sys.argv[2 + i]))

