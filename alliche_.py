'''
Implementation of the fatigue model for plain concrete [A.Alliche, 2004] under uniaxial compressive loading
(stress driven algorithm)
'''

import matplotlib.pyplot as plt
import numpy as np
import traits.api as tr
from bmcs_utils.api import InteractiveModel, \
    Int, Item, View, Float, Range, Button, ButtonEditor, mpl_align_yaxis, FloatRangeEditor
import sympy as sp

class AllicheModel(InteractiveModel):

    name = 'Alliche Model'

    # Loading parameters
    m = Int(50)       # number of increments
    n_cycle = Int(24500)
    n1 = Int(100)
    n2 = Int(50)
    n3 = Int(100)
    n4 = Int(50)
    n5 = Int(100)
    n6 = Int(50)
    sigma_u = Float(-120)
    max_H_level = Float(0.95)
    max_L_level = Float(0.90)
    order = Int(1)     # 1 for 'H-L' and 2 for 'L-H'
    b = Int(10)        # number_of_repeted_blocks

    #Alliche model parameters  (C120)
    lamda = Float(12500)
    mu = Float(18750)
    alpha = Float(2237.5)
    beta = Float(-2216.5)
    g = Float(-10.0)
    C0 = Float(0.00)
    C1 = Float(0.0019)
    K = Float(0.00485)
    n = Float(10)

    ipw_view = View(
        Item('m', latex='m'),
        Item('n_cycle', latex='n_cycle'),
        Item('n1', latex='n1'),
        Item('n2', latex='n2'),
        Item('n3', latex='n3'),
        Item('n4', latex='n2'),
        Item('n5', latex='n5'),
        Item('n6', latex='n6'),
        Item('sigma_u', latex='sigma_u'),
        Item('max_H_level', latex='max_H_level'),
        Item('max_L_level', latex='max_L_level'),
        Item('order', latex='order'),
        Item('b', latex='b'),

        Item('lamda', latex='lamda'),
        Item('mu', latex='mu'),
        Item('alpha', latex='alpha'),
        Item('beta', latex='beta'),
        Item('g', latex='g'),
        Item('C0', latex='C0'),
        Item('C1', latex='C1'),
        Item('K', latex='K'),
        Item('n', latex='n'),

    )

    def get_load(self):
        if self.order == 1:
        #   for H-L
            stress_level_1_max = self.max_H_level * self.sigma_u
            stress_level_2_max = self.max_L_level * self.sigma_u
            stress_level_3_max = self.max_L_level * self.sigma_u
            stress_level_4_max = self.max_L_level * self.sigma_u
            stress_level_5_max = self.max_L_level * self.sigma_u
            stress_level_6_max = self.max_L_level * self.sigma_u

        elif self.order == 2:
        #     for L-H
            stress_level_1_max = self.max_L_level * self.sigma_u
            stress_level_2_max = self.max_H_level * self.sigma_u
            stress_level_3_max = self.max_L_level * self.sigma_u
            stress_level_4_max = self.max_H_level * self.sigma_u
            stress_level_5_max = self.max_L_level * self.sigma_u
            stress_level_6_max = self.max_H_level * self.sigma_u

        stress_level_1_min = 0.2 * self.sigma_u
        stress_level_2_min = 0.2 * self.sigma_u
        stress_level_3_min = 0.2 * self.sigma_u
        stress_level_4_min = 0.2 * self.sigma_u
        stress_level_5_min = 0.2 * self.sigma_u
        stress_level_6_min = 0.2 * self.sigma_u
        # unloading_ratio = 0.0

        d_0 = np.zeros(1)

        d_1 = np.linspace(0, stress_level_1_max, self.n1 * 2)
        d_1.reshape(-1, 2)[:, 0] = stress_level_1_max
        # d_1[0] = 0.0
        d_1.reshape(-1, 2)[:, 1] = stress_level_1_min
        d_history_1 = d_1.flatten()
        sig_1_arr = np.hstack([np.linspace(d_history_1[i], d_history_1[i + 1], self.m, dtype=np.float_)
                               for i in range(len(d_1) - 1)])

        d_2 = np.linspace(0, stress_level_2_max, self.n2 * 2)
        d_2.reshape(-1, 2)[:, 0] = stress_level_2_max
        d_2.reshape(-1, 2)[:, 1] = stress_level_2_min
        d_history_2 = d_2.flatten()
        sig_2_arr = np.hstack([np.linspace(d_history_2[i], d_history_2[i + 1], self.m, dtype=np.float_)
                               for i in range(len(d_2) - 1)])

        d_3 = np.linspace(0, stress_level_3_max, self.n3 * 2)
        d_3.reshape(-1, 2)[:, 0] = stress_level_3_max
        d_3.reshape(-1, 2)[:, 1] = stress_level_3_min
        d_history_3 = d_3.flatten()
        sig_3_arr = np.hstack([np.linspace(d_history_3[i], d_history_3[i + 1], self.m, dtype=np.float_)
                               for i in range(len(d_3) - 1)])

        d_4 = np.linspace(0, stress_level_4_max, self.n4 * 2)
        d_4.reshape(-1, 2)[:, 0] = stress_level_4_max
        d_4.reshape(-1, 2)[:, 1] = stress_level_4_min
        d_history_4 = d_4.flatten()
        sig_4_arr = np.hstack([np.linspace(d_history_4[i], d_history_4[i + 1], self.m, dtype=np.float_)
                               for i in range(len(d_4) - 1)])

        d_5 = np.linspace(0, stress_level_5_max, self.n5 * 2)
        d_5.reshape(-1, 2)[:, 0] = stress_level_5_max
        d_5.reshape(-1, 2)[:, 1] = stress_level_5_min
        d_history_5 = d_5.flatten()
        sig_5_arr = np.hstack([np.linspace(d_history_5[i], d_history_5[i + 1], self.m, dtype=np.float_)
                               for i in range(len(d_5) - 1)])

        d_6 = np.linspace(0, stress_level_6_max, self.n6 * 2)
        d_6.reshape(-1, 2)[:, 0] = stress_level_6_max
        d_6.reshape(-1, 2)[:, 1] = stress_level_6_min
        d_history_6 = d_6.flatten()
        sig_6_arr = np.hstack([np.linspace(d_history_6[i], d_history_6[i + 1], self.m, dtype=np.float_)
                               for i in range(len(d_6) - 1)])

        sig_0_1_arr = np.linspace(
            d_0[-1], d_history_1[0], self.m, dtype=np.float_)
        sig_1_2_arr = np.linspace(
            d_history_1[-1], d_history_2[0], self.m, dtype=np.float_)
        sig_2_3_arr = np.linspace(
            d_history_2[-1], d_history_3[0], self.m, dtype=np.float_)
        sig_3_4_arr = np.linspace(
            d_history_3[-1], d_history_4[0], self.m, dtype=np.float_)
        sig_4_5_arr = np.linspace(
            d_history_4[-1], d_history_5[0], self.m, dtype=np.float_)
        sig_5_6_arr = np.linspace(
            d_history_5[-1], d_history_6[0], self.m, dtype=np.float_)

        sig_1_1_arr = np.linspace(
            sig_1_arr[-1], sig_1_arr[0], self.m, dtype=np.float_)
        sig_6_1_arr = np.linspace(
            sig_6_arr[-1], sig_1_arr[0], self.m, dtype=np.float_)
        sig_6_6_arr = np.linspace(
            sig_6_arr[-1], sig_6_arr[0], self.m, dtype=np.float_)

        sigma_1_arr = np.hstack(
            (d_0, sig_0_1_arr, sig_1_arr, sig_1_2_arr, sig_2_arr, sig_2_3_arr, sig_3_arr, sig_3_4_arr, sig_4_arr,
             sig_4_5_arr, sig_5_arr, sig_5_6_arr, sig_6_arr))

        sigma_arr = sigma_1_arr

        t_arr = np.linspace(0, 1, len(sigma_arr))

        return sigma_arr, t_arr

    def get_stress_strain(self, sigma_1_arr, lamda, mu, alpha, beta, g, C0, C1, K, n):

        # -----------------------------------------------------------------------
        # material parameters
        # -----------------------------------------------------------------------
        # lame constants [MPa]
        sigma_1_arr = self.get_load()[0]
        lamda = lamda
        mu = mu
        # fatigue model material parameter
        alpha = alpha
        beta = beta
        g = g
        C0 = C0
        C1 = C1
        K = K
        n = n

        # -----------------------------------------------------------------------
        # arrays to store the values
        # -----------------------------------------------------------------------
        # normal strain
        eps_1_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
        # lateral strain
        eps_2_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
        # damage factor
        w_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
        f_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
        D_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
        phi_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
        # -----------------------------------------------------------------------
        # state variables
        # -----------------------------------------------------------------------
        # sigma_1_arr[0] = 0
        eps_1_i = 0.0
        eps_2_i = 0.0
        w_i = 0.0
        D_i = 0.0

        for i in range(1, len(sigma_1_arr)):

            sigma_1_i = sigma_1_arr[i]

            eps_2_i = -1.0 * ((lamda + alpha * w_i) * sigma_1_i + g * w_i * (lamda + 2.0 * mu)) / \
                      ((lamda + 2.0 * mu) * (2.0 * (lamda + mu) + 4.0 *
                                             w_i * (alpha + beta)) - 2.0 * (lamda + alpha * w_i) ** 2)

            eps_1_i = sigma_1_i / \
                      (lamda + 2.0 * mu) - 2.0 * eps_2_i * \
                      (lamda + alpha * w_i) / (lamda + 2.0 * mu)

            f_i = abs(g) * eps_2_i - (C0 + 2 * C1 * w_i)

            kappa_i = (lamda + 2.0 * mu) * (2.0 * (lamda + mu) +
                                            4.0 * w_i * (alpha + beta) -
                                            alpha * (g / (2.0 * C1)) *
                                            (2.0 * eps_2_i + eps_1_i) -
                                            (g ** 2.0 / (2.0 * C1))) - 2.0 * (lamda + alpha * w_i) ** 2

            d_sigma_1 = sigma_1_arr[i] - sigma_1_arr[i - 1]
            m = -1.0 * ((lamda + alpha * w_i) / kappa_i) * d_sigma_1

            # loading stage (evolve of the fatigue damage based on (Marigo.85)
            # model)
            if m > 0:
                d_w = m * abs(g) / (2.0 * C1) * (f_i / K) ** n
            else:  # unloading stage (no fatigue damage)
                d_w = 0

            w_i = w_i + d_w

            # Energy release rate
            Y_norm = np.sqrt((-g * eps_1_i - alpha * (eps_1_i + 2. * eps_2_i) * eps_1_i
                              - 2. * beta * (eps_1_i ** 2.0)) ** 2.0 + 2.0 * (
                                         -g * eps_2_i - alpha * (eps_1_i + 2. * eps_2_i) * eps_1_i
                                         - 2 * beta * (eps_1_i ** 2.0)) ** 2.0)
            d_D = Y_norm  # * d_w
            D_i += d_D
            # print 'Y=', Y_norm
            # print('D=', D_i)

            # Helmholtz free energy
            phi_i = 0.5 * lamda * (eps_1_i + 2.0 * eps_2_i) ** 2.0 + mu * (
                        (eps_1_i) ** 2.0 + 2.0 * eps_2_i ** 2.0) + 2.0 * g * w_i * eps_2_i + alpha * \
                    (2.0 * w_i * eps_1_i * eps_2_i + 4.0 * w_i *
                     eps_2_i ** 2.0) + 4.0 * beta * w_i * eps_2_i ** 2.0

            if w_i > 5.0:
                print(' ----------> No Convergence any more')
                print(i)
                break

            if abs(eps_1_i) > 0.005:
                print(' ----------> No Convergence any more')
                print(i)
                break

            eps_1_arr[i] = eps_1_i
            eps_2_arr[i] = eps_2_i
            w_arr[i] = w_i
            f_arr[i] = f_i
            D_arr[i] = D_i
            phi_arr[i] = phi_i

        return sigma_1_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, D_arr, i, phi_arr


    # C120

    def get_output(self):

        sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, D_arr, inc, phi_arr = \
            self.get_stress_strain(self.get_load()[0], self.lamda, self.mu, \
                                   self.alpha, self.beta, self.g, self.C0, \
                                   self.C1, self.K, self.n)

        n = self.b * (self.n1 + self.n2 + self.n3 + self.n4 + self.n5 + self.n6)
        eps_1_max = np.zeros(n)
        eps_1_min = np.zeros(n)
        cycle = np.zeros(n)
        for i in range(0, n, 1):
            idx_1 = self.m + 2 * i * self.m - 1
            idx_2 = 2 * i * self.m
            if idx_1 <= len(eps_1_arr[0:inc]):
                idx_1 = idx_1
            else:
                idx_1 = self.m + 2 * (i - 1.0) * self.m - 1
                break

            if idx_2 <= len(eps_1_arr[0:inc]):
                idx_2 = idx_2
            else:
                idx_2 = 1 * (i - 1.0) * self.m
                break

            eps_1_max[i] = eps_1_arr[int(idx_1)]
            eps_1_min[i] = eps_1_arr[int(idx_2)]
            cycle[i] = i + 1
        # print('Number of cycles at Failure', i)
        return eps_1_max, eps_1_min, inc

    def subplots(self, fig):
        ax = fig.subplots(1, 1)
        return ax

    def update_plot(self, ax):
        inc = self.get_output()[2]
        ax.plot(self.get_load()[1][0:inc], abs(self.get_load()[0][0:inc]), 'k', color='blue', linewidth=0.1, alpha=1.0)
        # ax.title('loading history')
    #
    # # -----------------------------------------------------------------------
    # # plot 1
    # # -----------------------------------------------------------------------
    # plt.subplot(231)
    # plt.plot(get_load()[1][0:inc], abs(sigma_arr[0:inc]), 'k', color='blue', linewidth=0.1, alpha=1.0)
    # plt.title('loading history')
    #
    # plt.xlabel('Time')
    # plt.ylabel('$\sigma_{1}$')
    # # plt.legend(loc=4)
    #
    # # -----------------------------------------------------------------------
    # # plot 2
    # # -----------------------------------------------------------------------
    # plt.subplot(232)
    # plt.plot(abs(eps_1_arr[0:inc]), abs(sigma_arr[0:inc]), 'k', color='gray', linewidth=0.1, alpha=1.0)
    # plt.title('$ \epsilon_{11} - \sigma_{11}$')
    # plt.xlabel('$\epsilon_{11}$')
    # plt.ylabel('$\sigma_{11}$[MPa]')
    # # plt.legend(loc=4)
    #
    # # -----------------------------------------------------------------------
    # # plot 3
    # # -----------------------------------------------------------------------
    # plt.subplot(233)
    # plt.plot(abs(sigma_arr[0:inc]), w_arr[0:inc], 'k', color='black', linewidth=0.1, alpha=1.0)
    # plt.title('$ \epsilon_{11} - \sigma_{11}$')
    # plt.xlabel('$\epsilon_{11}$')
    # plt.ylabel('$\sigma_{11}$[MPa]')
    # # plt.legend(loc=4)
    #
    # PlotCalibration_n = b * (n1 + n2 + n3 + n4 + n5 + n6)
    # # -----------------------------------------------------------------------
    # # plot 4
    # # -----------------------------------------------------------------------
    # plt.subplot(234)
    #
    # eps_1_max = np.zeros(PlotCalibration_n)
    # eps_1_min = np.zeros(PlotCalibration_n)
    # cycle = np.zeros(PlotCalibration_n)
    # for i in range(0, PlotCalibration_n, 1):
    #     idx_1 = m + 2 * i * m - 1
    #     idx_2 = 2 * i * m
    #     if idx_1 <= len(eps_1_arr[0:inc]):
    #         idx_1 = idx_1
    #     else:
    #         idx_1 = m + 2 * (i - 1.0) * m - 1
    #         break
    #
    #     if idx_2 <= len(eps_1_arr[0:inc]):
    #         idx_2 = idx_2
    #     else:
    #         idx_2 = 1 * (i - 1.0) * m
    #         break
    #
    #     eps_1_max[i] = eps_1_arr[int(idx_1)]
    #     eps_1_min[i] = eps_1_arr[int(idx_2)]
    #     cycle[i] = i + 1
    #     f = open("eps_1_max.txt", "w+")
    #     #         D_eps_1_max = (str(eps_1_max[0:i])) - str((eps_1_max[0:i]))
    #     f.write(str(abs(eps_1_max[0:i])))
    # plt.plot(cycle[0:i], abs(eps_1_max[0:i]), 'k', color='pink', linewidth=1, alpha=1)
    #
    # # plt.plot(cycle, eps_2 / 1.5, color='gray', linewidth=1.0, alpha=0.5)
    # # plt.fill_between(cycle, eps_2, eps_2 / 1.5, facecolor='gray', alpha=0.5)
    # plt.xlabel('number of cycles')
    # plt.ylabel('max $\epsilon_{11}$')
    # # plt.ylim(0, 1)
    # # plt.legend()
    #
    # # -----------------------------------------------------------------------
    # # plot 5
    # # -----------------------------------------------------------------------
    # plt.subplot(235)
    # PlotCalibration_n = b * (n1 + n2 + n3 + n4 + n5 + n6)
    # w = np.zeros(PlotCalibration_n)
    # cycle = np.zeros(PlotCalibration_n)
    # for i in range(0, PlotCalibration_n, 1):
    #     idx = m + 2 * i * m - 1
    #     if idx <= len(w_arr[0:inc]):
    #         idx = idx
    #     else:
    #         idx = m + 2 * (i - 1.0) * m - 1
    #         break
    #
    #     w[i] = w_arr[idx]
    #     cycle[i] = i + 1
    #
    # plt.plot(cycle[0:i], w[0:i], 'k', color='green', linewidth=1, alpha=1)
    # plt.xlabel('number of cycles')
    # plt.ylabel('Damage')
    # # plt.ylim(0, 1)
    # # plt.legend()
    #
    # # -----------------------------------------------------------------------
    # # plot 6
    # # -----------------------------------------------------------------------
    # plt.subplot(236)
    # D = np.zeros(PlotCalibration_n)
    # cycle = np.zeros(PlotCalibration_n)
    # for i in range(0, PlotCalibration_n, 1):
    #     idx = m + 2 * i * m - 1
    #     if idx <= len(D_arr[0:inc]):
    #         idx = idx
    #     else:
    #         idx = m + 2 * (i - 1.0) * m - 1
    #         break
    #
    #     D[i] = D_arr[idx]
    #     cycle[i] = i + 1
    #
    # plt.plot(cycle[0:i], D[0:i], 'k', color='black', linewidth=1, alpha=1)
    # plt.xlabel('number of cycles')
    # plt.ylabel('$Dissipation$')
    #
    # fig = plt.gcf()
    # #     fig.set_size_inches(8, 8)
    #
    # #     manager = plt.get_current_fig_manager()
    # #     manager.window.showMaximized()
    # plt.grid(True)
    # plt.tight_layout()
    # #     plt.tight_layout(pad=1.08, h_pad=1 , w_pad=1 , rect=(1, 1, 1, 1))
    #
    # plt.show()