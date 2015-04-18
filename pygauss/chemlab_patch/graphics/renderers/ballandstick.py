# CJS changed relative paths to chemlab ones

from chemlab.graphics.renderers.base import AbstractRenderer
from chemlab.graphics.renderers.atom import AtomRenderer
from chemlab.graphics.renderers.bond import BondRenderer

from chemlab.db import ChemlabDB
cdb = ChemlabDB()

# CJS added option to color atoms by predefined colors
class BallAndStickRenderer(AbstractRenderer):
    '''Render a ball and stick representation of a series of
    coordinates and bonds.
    
    .. image:: /_static/ballandstick_renderer.png
    
    **Parameters**
    
    widget:
        The parent QChemlabWidget
    r_array: np.ndarray((NATOMS, 3), dtype=float)
        The coordinate array
    type_array: np.ndarray((NATOMS, 3), dtype=object)
        An array containing all the atomic symbols like `Ar`, `H`, `O`.
        If the atomic type is unknown, use the `Xx` symbol.
    bonds: np.ndarray((NBONDS, 2), dtype=int)
        An array of integer pairs that represent the bonds.


    '''
    def __init__(self, widget, r_array, type_array, bonds, shading='phong', rgba_array=[]):
        super(BallAndStickRenderer, self).__init__(widget)
        vdw_dict = cdb.get("data", 'vdwdict')        
        
        scale = 0.3
        for k in vdw_dict:
            vdw_dict[k] = vdw_dict[k] * scale
        
        self.has_bonds = len(bonds) > 0
        
        self.ar = AtomRenderer(widget, r_array, type_array,
                               radii_map = vdw_dict, shading=shading, rgba_array=rgba_array)
        
        if self.has_bonds:
            self.br = BondRenderer(widget, bonds, r_array, type_array, style='impostors',
                                   shading=shading, rgba_array=rgba_array)

        
    def draw(self):
        self.ar.draw()
        
        if self.has_bonds:
            self.br.draw()

    def update_positions(self, r_array):
        '''Update the coordinate array r_array'''
        self.ar.update_positions(r_array)
        
        if self.has_bonds:
            self.br.update_positions(r_array)