import traits.api as tr
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

    A = Float(-10.66)
    B = Float(6.1)
    C = Float(2)

    eta_uparrow = Float(0.74)
    beta_uparrow = Float(74.7)
    eta_downarrow = Float(0.590)
    beta_downarrow = Float(60.5)

    s_max_i = Float(0.90)
    s_min_i = Float(0.20)
    s_max_ip1 = Float(0.85)
    s_min_ip1 = Float(0.80)

    ipw_view = View(
        Item('A', latex='A'),
        Item('B', latex='B'),
        Item('C', latex='C'),
        Item('eta_uparrow', latex='\eta_{up}'),
        Item('beta_uparrow', latex='\beta_{up}'),
        Item('eta_downarrow', latex='\eta_{down}'),
        Item('beta_downarrow', latex='\beta_{down}'),
        Item('s_max_i', latex='s_{max}_i'),
        Item('s_min_i', latex='s_{min}_i'),
        Item('s_max_ip1', latex='s_{max}_ip1'),
        Item('s_min_ip1', latex='s_{min}_ip1')
    )

    def get_eta(self):
        eta = sum(self.eta_i) + (self.delta_eta_i)
        return eta

    def get_eta_i(self):
        eta_i = self.n_i / self.n_i_f
        return eta_i

    def get_s_bar_i(self):
        s_bar_i = (self.s_m_i + self.s_m_ip1) / 2
        return s_bar_i

    def get_s_m_i(self):
        s_m_i = (self.s_max_i + self.s_min_i) / 2
        return s_m_i

    def get_s_m_ip1(self):
        s_m_ip1 = (self.s_max_ip1 + self.s_min_ip1) / 2
        return s_m_ip1

    def get_delta_s_max_i(self):
        delta_s_max_i = (self.s_max_ip1 - self.s_min_i)
        return delta_s_max_i

    def get_delta_s_min_i(self):
        delta_s_min_i = (self.s_min_ip1 - self.s_min_i)
        return delta_s_min_i

    def get_delta_eta_i(self):
        if (self.tilde_eta_i > 0) & (self.tilde_eta_i <= self.eta_x):
            delta_eta_i = (self.delta_eta_max * (1 - ((self.eta_x - self.tilde_eta_i) / self.eta_x)))
        elif (self.tilde_eta_i > self.eta_x) & (self.ilde_eta_i < 1):
            delta_eta_i = (self.delta_eta_max * ((self.tilde_eta_i - 1) / (self.eta_x - 1)))
        return delta_eta_i

    def get_delta_eta_max_i(self):
        delta_eta_max_i = (self.f_1 * (self.delta_s_max_i + self.f_2)) * np.sign(self.delta_s_max_i)
        return delta_eta_max_i

    def get_f_1(self):
        f_1 = self.A * (self.delta_s_max_i) ** 2 + self.B * self.delta_s_max_i * np.sign(self.delta_s_max_i)
        return f_1

    def get_f_2(self):
        f_2 = self.C * (0.475 - self.s_bar_i)
        return f_2

    def get_eta_x(self):
        if self.delta_eta_max > 0:
            eta_x = (self.eta_uparrow + self.delta_eta_max / np.tan(self.beta_uparrow))
        else:
            eta_x = (self.eta_downarrow + self.delta_eta_max / np.tan(self.beta_downarrow))
        return eta_x

    def get_array(self):
        sum = self.L + self.n_x + self.F
        array = range(sum)
        return array

    def subplots(self, fig):
        ax = fig.subplots(1, 1)
        return ax

    #
    def update_plot(self, ax):
        ax.plot(self.get_array(), color='blue', label='$w$ [mm]')
        ax.set_ylabel(r'$w [\mathrm{mm}]$')
        ax.legend(loc='lower right')
