import numpy as np
import modeling_of_uncertainty as moe
import normal_probability_distribution as npd
from scipy import stats as st 
import matplotlib.pyplot as plt

#Load and resistence normal variables - single load
#R = Resistence array variables
#S = Load array variables
#k_R = How much times the nominal resistance is below the mean resistance
#k_S = How much times the nominal load is above the mean load

def single_load_norm_beta_value(R_mean, S_mean, R_std, S_std): 
    beta = (R_mean - S_mean)/(R_std**2 + S_std**2)**0.5

    return beta

def p_failure_norm_single_load(R_mean, S_mean, R_std, S_std):
    beta = single_load_norm_beta_value(R_mean, S_mean, R_std, S_std)
    p_f = 1 - npd.phi(beta, 0, 1)

    print('The probability of failure is {}% with beta-value of {}.'.format(p_f*100, beta))

    return p_f, beta

def single_load_norm_factors(R_mean, S_mean, R_std, S_std, k_R, k_S):  
    e = 0.75
    delta_R = R_std/R_mean
    delta_S = S_std/S_mean
    beta = single_load_norm_beta_value(R_mean, S_mean, R_std, S_std)

    phi_factor = (1 - e*beta*delta_R)/(1 - k_R*delta_R)
    gamma_factor = (1 + e*beta*delta_S)/(1 - k_S*delta_S)

    print('\n\nNominal Capacity Reduction Factor = {}\nLoad Factor is {}\n\n'.format(phi_factor, gamma_factor))

    return(phi_factor, gamma_factor)

#Load and resistence normal variables - multiple loads
def multiple_load_norm_beta(R_mean, S_mean_array, R_std, S_std_array):
    S_sum_std_array_square = []
    for S_std in S_std_array:
        S_sum_std_array_square.append(S_std**2)
    beta = (R_mean - sum(S_mean_array))/(R_std**2 + sum(S_sum_std_array_square))**0.5

    return beta

def p_failure_norm_multiple_load(R_mean, S_mean_array, R_std, S_std_array):
    beta = multiple_load_norm_beta(R_mean, S_mean_array, R_std, S_std_array)

    p_f = 1 - npd.phi(beta, 0, 1)
    print('\n\nProbability of failure = {}%\nBeta-value = {}\n\n.'.format(p_f*100, beta))

    return p_f, beta


def multiple_load_norm_factors(R_mean, S_mean_array, R_std, S_std_array, k_R, k_S_array):
    beta = multiple_load_norm_beta(R_mean, S_mean_array, R_std, S_std_array)
    S_sum_std_array_square = []
    for S_std in S_std_array:
        S_sum_std_array_square.append(S_std**2)
    enn = sum(S_sum_std_array_square)**0.5/sum(S_std_array)
    e = 0.75
    delta_R = R_std/R_mean
    phi_factor = (1 - e*beta*delta_R)/(1 - k_R*delta_R)
    gamma_factor = []

    for i in range(0, len(S_mean_array)):
        delta_S_array = S_std_array[i]/S_mean_array[i]
        gamma_factor.append((1 + e*enn*delta_S_array*beta)/(1 + k_S_array[i]*delta_S_array))

    print('\n\nNominal Reduction Factor = {}.\nLoad Factors = {}\n\n'.format(phi_factor, gamma_factor))
    

multiple_load_norm_factors(4962.16, [750.06, 775.845], 645.08, [97.47, 287.01], 2, [2, 2])