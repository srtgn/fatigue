'''
Created on 14.06.2017  @author: abaktheer

Update on Jun 11, 2020 @author: Saeed

Implementation of the fatigue model for plain concrete [A.Alliche, 2004] under uniaxial compressive two-block loading loading
(stress driven algorithm)

'''


from array import array
import matplotlib.pyplot as plt
import numpy as np
import operator     
from Functions import*
from prompt_toolkit import input







H = [1 ,   1  ,  1  ,  4  ,  7  ,  15  ,  22 ,   30  ,  37]
L = [543  ,  1086   , 2172   , 5430   , 10860   , 21720   , 32580 ,   43440 ,   54300]

# Study = alliche_model(m = 50, n1 =300000 , n2=1, 
#                           sigma_u= -120, H_max_level = 0.95, 
#                           L_max_level = 0.725, min_level = 0.20,rep_num = 1)
#       
#     # sigma_arr,t_arr = Study.get_sigma_arr_2step()
# sigma_arr = Study.get_sigma_arr_2block()
#     # print(sigma_arr)
#       
# #     eps_1_i, eps_2_i, w_i, D_i, eps_1_arr, eps_2_arr, w_arr, f_arr, D_arr, inc, cycle, phi_arr = Study.get_stress_strain(sigma_arr)
# Study.get_stress_strain(sigma_arr)
#  
# input('break')
# rep_num is the number of repetition of the two block loading scenario
for i in range (6,9):
#         
    Study = alliche_model(m = 50, n1 =L[i] , n2=H[i], 
                          sigma_u= -120, H_max_level = 0.65, 
                          L_max_level = 0.95, min_level = 0.20,rep_num = 5)
    # sigma_arr,t_arr = Study.get_sigma_arr_2step()
    sigma_arr = Study.get_sigma_arr_2block()
    # print(sigma_arr)
     
#     eps_1_i, eps_2_i, w_i, D_i, eps_1_arr, eps_2_arr, w_arr, f_arr, D_arr, inc, cycle, phi_arr = Study.get_stress_strain(sigma_arr)
    Study.get_stress_strain(sigma_arr)
    # Study.plot(t_arr,sigma_arr,w_arr,D_arr,inc,eps_1_arr)
     
