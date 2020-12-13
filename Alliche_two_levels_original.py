'''
Created on 14.06.2017

@author: abaktheer
'''

'''
Implementation of the fatigue model for plain concrete [A.Alliche, 2004] under uniaxial compressive loading
(stress driven algorithm)
'''


'''
To do#

1. model class
2. loading scenario class (reduce repetition)
3. improve printing
'''

import matplotlib.pyplot as plt
import numpy as np
import pathlib

# # printing the results
# d_c_65 = np.loadtxt(
#     r'D:\IMB\DrBaktheer\reports\06102020\damage_65_Constant.txt')
# d_c_95 = np.loadtxt(
#     r'D:\IMB\DrBaktheer\reports\06102020\damage_95_Constant.txt')
# d_6595 = np.loadtxt(
#     r'D:\IMB\DrBaktheer\reports\06102020\damage_6595.txt')
# d_9565 = np.loadtxt(
#     r'D:\IMB\DrBaktheer\reports\06102020\damage_9565.txt')
#
# s2_d_c_65 = np.loadtxt(
#     r'D:\IMB\DrBaktheer\reports\06102020\s2_damage_65_Constant.txt')
# s2_d_c_95 = np.loadtxt(
#     r'D:\IMB\DrBaktheer\reports\06102020\s2_damage_95_Constant.txt')
# s2_d_6595 = np.loadtxt(
#     r'D:\IMB\DrBaktheer\reports\06102020\s2_damage_6595.txt')
# s2_d_9565 = np.loadtxt(
#     r'D:\IMB\DrBaktheer\reports\06102020\s2_damage_9565.txt')
#
# plt.plot(np.arange(0,len(d_c_65),1)/len(d_c_65),d_c_65, label='d_c_65')
# plt.plot(np.arange(0,len(d_c_95),1)/len(d_c_95),d_c_95, label='d_c_95')
# # plt.plot(np.arange(0,len(d_6595),1)/len(d_6595),d_6595, label='d_6595')
#
# plt.xlabel('number of cycles')
# plt.ylabel('$Dissipation$')
# plt.title('Study 1')
# plt.legend()
# plt.show()
#
# plt.plot(np.arange(0,len(s2_d_c_65),1)/len(s2_d_c_65),s2_d_c_65, label='s2_d_c_65')
# plt.plot(np.arange(0,len(s2_d_c_95),1)/len(s2_d_c_95),s2_d_c_95, label='s2_d_c_95')
# plt.xlabel('number of cycles')
# plt.ylabel('$Dissipation$')
# plt.title('Study 2')
# plt.legend()
# plt.show()
# #
# # value = input("Stop?! :")


def get_stress_strain(sigma_1_arr, lamda, mu, alpha, beta, g, C0, C1, K, n):

    #-----------------------------------------------------------------------
    # arrays to store the values
    #-----------------------------------------------------------------------
    # normal strain
    eps_1_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    # lateral strain
    eps_2_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    # damage factor
    w_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    f_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    D_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)
    phi_arr = np.zeros_like(sigma_1_arr, dtype=np.float_)

    #-----------------------------------------------------------------------
    # material parameters
    #-----------------------------------------------------------------------
    # lame constants [MPa]
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

    #-----------------------------------------------------------------------
    # state variables
    #-----------------------------------------------------------------------
    #sigma_1_arr[0] = 0
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
            (lamda + 2.0 * mu) - 2.0 * eps_2_i *  \
            (lamda + alpha * w_i) / (lamda + 2.0 * mu)

        f_i = abs(g) * eps_2_i - (C0 + 2 * C1 * w_i)

        kappa_i = (lamda + 2.0 * mu) * (2.0 * (lamda + mu) +
                                        4.0 * w_i * (alpha + beta) -
                                        alpha * (g / (2.0 * C1)) *
                                        (2.0 * eps_2_i + eps_1_i) -
                                        (g**2.0 / (2.0 * C1))) - 2.0 * (lamda + alpha * w_i)**2

        d_sigma_1 = sigma_1_arr[i] - sigma_1_arr[i - 1]
        m = -1.0 * ((lamda + alpha * w_i) / kappa_i) * d_sigma_1

        # loading stage (evolve of the fatigue damage based on (Marigo.85)
        # model)
        if m > 0:
            d_w = m * abs(g) / (2.0 * C1) * (abs(f_i / K)**n)
            
#             d_w = m * abs(g) / (2.0 * C1) * (abs(f_i / K)**n)*np.sign(f_i / K)

#             d_w = m * abs(g) / (2.0 * C1) * (f_i / K)**n
        else:  # unloading stage (no fatigue damage)
            d_w = 0

        w_i = w_i + d_w

        # Energy release rate
        Y_norm = np.sqrt((-g * eps_1_i - alpha * (eps_1_i + 2. * eps_2_i) * eps_1_i
                          - 2. * beta * (eps_1_i**2.0))**2.0 + 2.0 * (-g * eps_2_i - alpha * (eps_1_i + 2. * eps_2_i) * eps_1_i
                                                                      - 2 * beta * (eps_1_i**2.0))**2.0)
        d_D = Y_norm  * d_w
        D_i += d_D
        # print 'Y=', Y_norm
        #print('D=', D_i)

        # Helmholtz free energy
        phi_i = 0.5 * lamda * (eps_1_i + 2.0 * eps_2_i)**2.0 + mu * ((eps_1_i)**2.0 + 2.0 * eps_2_i**2.0) + 2.0 * g * w_i * eps_2_i + alpha * \
            (2.0 * w_i * eps_1_i * eps_2_i + 4.0 * w_i *
             eps_2_i**2.0) + 4.0 * beta * w_i * eps_2_i**2.0

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


# H-L with eta-H (20, 25, 30, 35, 40, 45 ,50 )
# L-H with eta-L (20, 25, 30, 35, 40, 45 ,50 )

Mult =[20, 25, 30, 35, 40, 45 ,50]

N_H = [element * int(957 / 100) for element in Mult]
N_L = [element * int(9406 / 100) for element in Mult]


for num in N_H:
    print(num)

    if __name__ == '__main__':
    
        m = 50  # number of increments in each cycle
    
        n1 = num
        n2 = 10

        sigma_u = - 64.5

        max_s_level = 0.85
        min_s_level = 0.75
        
        stress_level_1_max = max_s_level * sigma_u
        stress_level_2_max = min_s_level * sigma_u
    
        stress_level_1_min = 0.05 * sigma_u
        stress_level_2_min = 0.05 * sigma_u
    
    #==========================================================
        d_0 = np.zeros(1)
    
        d_1 = np.linspace(0, stress_level_1_max, n1 * 2)
        d_1.reshape(-1, 2)[:, 0] = stress_level_1_max
        d_1.reshape(-1, 2)[:, 1] = stress_level_1_min
        d_history_1 = d_1.flatten()
    
        d_2 = np.linspace(0, stress_level_2_max, n2 * 2)
        d_2.reshape(-1, 2)[:, 0] = stress_level_2_max
        d_2.reshape(-1, 2)[:, 1] = stress_level_2_min
        d_history_2 = d_2.flatten()
    
        d_arr = np.hstack((d_0, d_history_1, d_history_2))
    
        sigma_arr = np.hstack([np.linspace(d_arr[i], d_arr[i + 1], m,  endpoint=False)
                               for i in range(len(d_arr) - 1)])
    
        t_arr = np.linspace(0, 1, len(sigma_arr))
    
        #================================================================
    
        # C120
        sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, D_arr, inc, phi_arr = get_stress_strain(sigma_arr,
        lamda = 8888.889,
        mu = 13333.33,
        alpha = 3000,
        beta = -2350,
        g = -9.8,
        C0 = 0.00,
        C1 = 0.002,
        K = 0.00504,
        n = 10
        )
        
    #         sigma_arr, lamda=8888.889, mu=13333.33, alpha=2150, beta=-2300, g=-9.78,
    #         C0=0.00, C1=0.002, K=0.0035, n=18)
        print('inc', int(inc/100))
        
    #     sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, inc = get_stress_strain(
    # sigma_arr, lamda=13972.2, mu=20958.3, alpha=2237.5, beta=-2216.5,
    # g=-10.0, C0=0.00, C1=0.00188, K=0.003345, PlotCalibration_n=10)
    
    
    #     # C120
    #     sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, inc = get_stress_strain(
    #         sigma_arr, lamda=12500, mu=18750, alpha=2237.5, beta=-2216.5, g=-10.0,
    #         C0=0.00, C1=0.0019, K=0.00485, PlotCalibration_n=10)
    
    #     # C80 - alliche paper
    #     sigma_arr, eps_1_arr, eps_2_arr, w_arr, f_arr, inc = get_stress_strain(
    # sigma_arr, lamda=10555.55, mu=15833.33, alpha=2237.5, beta=-2216.5,
    # g=-9.788, C0=0.00, C1=0.002033, K=0.003345, PlotCalibration_n=10)
#     
#         #-----------------------------------------------------------------------
#         # plot 1
#         #-----------------------------------------------------------------------
#         plt.subplot(231)
#         plt.plot(t_arr[0:inc], abs(sigma_arr[0:inc]), 'k', linewidth=1, alpha=1.0)
#         plt.title('loading history')
#
#         plt.xlabel('Time')
#         plt.ylabel('$\sigma_{1}$')
#         # plt.legend(loc=4)
#
#         #-----------------------------------------------------------------------
#         # plot 2
#         #-----------------------------------------------------------------------
#         plt.subplot(232)
#         plt.plot(abs(eps_1_arr[0:inc]), abs(
#             sigma_arr[0:inc]), 'k', linewidth=1, alpha=1.0)
#         plt.title('$ \epsilon_{11} - \sigma_{11}$')
#         plt.xlabel('$\epsilon_{11}$')
#         plt.ylabel('$\sigma_{11}$[MPa]')
#         # plt.legend(loc=4)
#
#         #-----------------------------------------------------------------------
#         # plot 2
#         #-----------------------------------------------------------------------
#         plt.subplot(236)
#         plt.plot(abs(sigma_arr[0:inc]), w_arr[0:inc], 'k', linewidth=1, alpha=1.0)
#         plt.title('$ \epsilon_{11} - \sigma_{11}$')
#         plt.xlabel('$\epsilon_{11}$')
#         plt.ylabel('$\sigma_{11}$[MPa]')
#         # plt.legend(loc=4)
#         #-----------------------------------------------------------------------
#         # plot 3
#         #-----------------------------------------------------------------------
#     #     plt.subplot(222)
#     #     plt.plot(abs(sigma_arr), w_arr, 'k', linewidth=1, alpha=1)
#     #     plt.xlabel('$\sigma_{11}$')
#     #     plt.ylabel('Damage')
#     #     #plt.ylim(0, 1)
#     #     # plt.legend()
#
#         #-----------------------------------------------------------------------
#         # plot 4
#         #-----------------------------------------------------------------------
#     #     plt.subplot(223)
#     #     plt.plot(t_arr, w_arr, 'b', linewidth=1, alpha=1)
#     #     plt.xlabel('Time')
#     #     plt.ylabel('Damage')
#     #     #plt.ylim(0, 1)
#     #     # plt.legend()
#
#         #-----------------------------------------------------------------------
#         # plot 5
#         #-----------------------------------------------------------------------
#     #     plt.subplot(225)
#     #
        n = (n1 + n2)
#     #     eps_1 = np.zeros(PlotCalibration_n)
#     #     cycle = np.zeros(PlotCalibration_n)
#     #     sig_1 = np.zeros(PlotCalibration_n)
#     #
#     #     for i in range(0, PlotCalibration_n, 1):
#     #         idx = m + 2 * i * m  # - 1
#     #         eps_1[i] = eps_1_arr[idx]
#     #         cycle[i] = i
#     #         #sig_1[i] = sigma_1_arr[idx]
#     #         # print sig_1
#     #     plt.plot(cycle, eps_1, 'b', linewidth=1, alpha=1)
#     #     plt.xlabel('number of cycles')
#     #     plt.ylabel('max $\epsilon_{11}$')
#         # plt.legend()
#
#     #     #-----------------------------------------------------------------------
#     #     # plot 5
#     #     #-----------------------------------------------------------------------
#         plt.subplot(234)
#
        eps_1_max = np.zeros(n)
        eps_1_min = np.zeros(n)
        cycle = np.zeros(n)
        for i in range(0, n, 1):

            idx_1 = m + 2 * i * m
            idx_2 = 2 * i * m
            if idx_1 <= len(eps_1_arr[0:inc]):
                idx_1 = idx_1
            else:
                idx_1 = m + 2 * (i - 1.0) * m
                break

            if idx_2 <= len(eps_1_arr[0:inc]):
                idx_2 = idx_2
            else:
                idx_2 = 1 * (i - 1.0) * m
                break

    #         print('idx_1', idx_1)

            eps_1_max[i] = eps_1_arr[int(idx_1)]
            eps_1_min[i] = eps_1_arr[int(idx_2)]
            cycle[i] = i + 1

#         plt.plot(cycle[0:i], abs(eps_1_max[0:i]), 'k', linewidth=1, alpha=1)
#
#         plt.xlabel('number of cycles')
#         plt.ylabel('max $\epsilon_{11}$')
#
#     # p =  int(max_s_level * 100)
#     # g = open(r"D:\IMB\DrBaktheer\reports\06102020\s2_damage_%s.txt"  % p , "w+")
#     # np.savetxt(g, D[0:i])
#
#
#         #     #-----------------------------------------------------------------------
#     #     # plot 6
#     #     #-----------------------------------------------------------------------
#         plt.subplot(235)
#         w = np.zeros(n)
#         cycle = np.zeros(n)
#         for i in range(0, n, 1):
#             idx = m + 2 * i * m
#             if idx <= len(w_arr[0:inc]):
#                 idx = idx
#             else:
#                 idx = m + 2 * (i - 1.0) * m
#                 break
#
#             w[i] = w_arr[idx]
#             cycle[i] = i + 1
#
#         plt.plot(cycle[0:i], w[0:i], 'k', linewidth=1, alpha=1)
#         plt.xlabel('number of cycles')
#         plt.ylabel('Damage')
#         #plt.ylim(0, 1)
#         # plt.legend()
#
#
#     #     #-----------------------------------------------------------------------
#     #     # plot 6
#     #     #-----------------------------------------------------------------------
#         plt.subplot(233)
#         D = np.zeros(n)
#         cycle = np.zeros(n)
#         for i in range(0, n, 1):
#             idx = m + 2 * i * m
#             if idx <= len(D_arr[0:inc]):
#                 idx = idx
#             else:
#                 idx = m + 2 * (i - 1.0) * m
#                 break
#
#             D[i] = D_arr[idx]
#             cycle[i] = i + 1
#
#         plt.plot(cycle[0:i], D[0:i], 'k', linewidth=1, alpha=1)
#         plt.xlabel('number of cycles')
#         plt.ylabel('$Dissipation$')
      
        #=========================================================================
        # saving results
        #=========================================================================
#         eps_max_record = np.zeros(1)
#         eps_min_record = np.zeros(1)
#         N_record = np.zeros(1)
#         w_record = np.zeros(1)
#         stiffness_record = np.zeros(1)
#        
#         for i in range(0, i, 1):
#             eps_max_record = np.vstack((eps_max_record, eps_1_max[i]))
#             eps_min_record = np.vstack((eps_min_record, eps_1_min[i]))
#             N_record = np.vstack((N_record, cycle[i]))
#             w_record = np.vstack((w_record, w[i]))
#             stiffness_record = np.vstack(
#                 (stiffness_record, stress_level_1_max / eps_1_max[i]))
#         p =  (str(int(max_s_level * 100)) + str(int(min_s_level * 100)))
#         p =  int(max_s_level * 100)
#         g = open(r"D:\IMB\DrBaktheer\reports\06102020\s2_damage_%s.txt"  % p , "w+")
#         np.savetxt(g, D[0:i])


        # plt.show()

    # sub_p =  (str(int(max_s_level * 100)) + str(int(min_s_level * 100))) + '_' + str(int(inc/100)) + '_' + str(num)
    #
    # p = pathlib.Path("09112020 report outputs") / "H_L"
    # p.mkdir(parents=True, exist_ok=True)
    #
    # np.savetxt(p / ("eps_1_max_%s.txt" % sub_p), abs(eps_1_max[0:i]), delimiter=" ", fmt="%s")
