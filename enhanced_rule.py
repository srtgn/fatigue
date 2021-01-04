import traits.api as tr
from bmcs_beam.beam_config.boundary_conditions import BoundaryConditions
from bmcs_cross_section.cs_design import CrossSectionDesign
from bmcs_utils.api import InteractiveModel, \
    Int, Item, View, Float, Range, Button, ButtonEditor, mpl_align_yaxis, FloatRangeEditor
from sympy.physics.continuum_mechanics.beam import Beam
import sympy as sp
import numpy as np
from numbers import Number

class EnhancedRule(InteractiveModel):

    name = 'Enhanced Rule'

    n_x = Int(100)
    L = Int(100)
    F = Int(100)
    F = Int(100)

    ipw_view = View(
        Item('L', latex='L \mathrm{[mm]}'),
        Item('F', latex='F \mathrm{[N]}'),
        Item('n_x', latex='n_x'),
           )

    def get_array(self):
        sum = self.L + self.n_x + self.F
        array = range (sum)
        return array

    def subplots(self, fig):
        ax = fig.subplots(1, 1)

        return ax
#
    def update_plot(self, ax):
        ax.plot(self.get_array(), color='blue', label='$w$ [mm]')
        ax.set_ylabel(r'$w [\mathrm{mm}]$')
        ax.legend(loc='lower right')
