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

    A = Float(-10.66)
    B = Float(6.1)
    C = Float(2)

    eta_uparrow = Float(0.74)
    beta_uparrow = Float(74.7)
    eta_downarrow = Float(0.590)
    beta_downarrow = Float(60.5)

    s = {"s_max_i" : 0.90, "s_min_i" : 0.20, "s_max_ip1" : 0.85, "s_min_ip1" : 0.20}
    s_max_i = list(s.values())[0]
    s_min_i = list(s.values())[1]
    s_max_ip1 = list(s.values())[2]
    s_min_ip1 = list(s.values())[3]

    n_list = [500, 5000, 50000]

    ipw_view = View(
        Item('A', latex='A'),
        Item('B', latex='B'),
        Item('C', latex='C'),
        Item('eta_uparrow', latex='\eta_{up}'),
        Item('beta_uparrow', latex='\beta_{up}'),
        Item('eta_downarrow', latex='\eta_{down}'),
        Item('beta_downarrow', latex='\beta_{down}'),
        Item('s_max_i', latex='\s_{max}_i'),
        Item('s_min_i', latex='\s_{min}_i'),
        Item('s_max_ip1', latex='\s_{max}_ip1'),
        Item('s_min_ip1', latex='\s_{min}_ip1')
    )


    ''' The variable of first part of the eq27'''

    def get_eta_i(self):
        eta_i = []
        for i, n in enumerate(self.n_list):
            n_i = n
            n_i_f = sum(self.n_list)
            eta_i.append(n_i / n_i_f)
        return eta_i

    
    ''' The variables of the second part of the eq27'''
    
    #s_bar_i
    def get_s_m_i(self):
        s_m_i = (self.s_max_i + self.s_min_i) / 2
        return s_m_i

    def get_s_m_ip1(self):
        s_m_ip1 = (self.s_max_ip1 + self.s_min_ip1) / 2
        return s_m_ip1
    
    def get_s_bar_i(self):
        s_bar_i = (self.get_s_m_i() + self.get_s_m_ip1()) / 2
        return s_bar_i


    #delta_s_max_i & delta_s_min_i
    def get_delta_s_max_i(self):
        delta_s_max_i = (self.s_max_ip1 - self.s_min_i)
        return delta_s_max_i

    def get_delta_s_min_i(self):
        delta_s_min_i = (self.s_min_ip1 - self.s_min_i)
        return delta_s_min_i
    
    def get_tilde_eta_i(self):
        tilde_eta_i = self.get_eta_i()[0]
        return tilde_eta_i

    def get_delta_eta_max_i(self):
        delta_s_max_i = self.get_delta_s_max_i()
        f_1 = self.A * (delta_s_max_i) ** 2 + self.B * delta_s_max_i * np.sign(delta_s_max_i)
        f_2 = self.C * (0.475 - self.get_s_bar_i())
        delta_eta_max_i = (f_1 * (delta_s_max_i + f_2)) * np.sign(delta_s_max_i)
        return delta_eta_max_i

    def get_eta_x(self):
        if self.get_delta_eta_max_i() > 0:
            eta_x = (self.eta_uparrow + self.get_delta_eta_max_i() / np.tan(self.beta_uparrow))
        else:
            eta_x = (self.eta_downarrow + self.get_delta_eta_max_i() / np.tan(self.beta_downarrow))
        return eta_x

    def get_delta_eta_i(self):
        eta_x = self.get_eta_x()
        delta_eta_max = self.get_delta_eta_max_i()
        if (self.get_tilde_eta_i() > 0) & (self.get_tilde_eta_i() <= eta_x):
            delta_eta_i = (delta_eta_max * (1 - ((eta_x - self.get_tilde_eta_i()) / eta_x)))
        elif (self.get_tilde_eta_i() > eta_x) & (self.get_tilde_eta_i() < 1):
            delta_eta_i = (delta_eta_max * ((self.get_tilde_eta_i() - 1) / (eta_x - 1)))
        return delta_eta_i


    ''' eq 27:
        The expression for the consumed fatigue life under a loading scenario consisting 
        of n varying loading ranges is proposed in the form '''
    
    def get_eta(self):
        eta = sum(self.get_eta_i()) + (self.get_delta_eta_i())
        return eta
    
    def get_array(self):
        sum = self.L + self.n_x + self.F
        array = range(sum)
        return array

    def subplots(self, fig):
        ax = fig.subplots(1, 1)
        return ax

    def update_plot(self, ax):
        ax.bar(10, self.get_eta(), 1, color='blue', label='$w$ [mm]')
        ax.set_ylabel(r'$w [\mathrm{mm}]$')
        ax.legend(loc='lower right')
