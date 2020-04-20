import numpy as np
import modeling_of_uncertainty as moe
import normal_probability_distribution as npd
from scipy import stats as st 
import matplotlib.pyplot as plt
import scipy as sp


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
    
def afosm_l_ls(R_mean, S_mean, R_std, S_std): #reduced hasofer-lind method for linear limit state

    beta_HL = (R_mean - S_mean)/(R_std**2 + S_std**2)**0.5

    R_line = np.arange(-R_mean, R_mean + 2*R_mean/100, 2*R_mean/100)

    def S_line_calc(R_line):
        S_line = (R_std*R_line + R_mean - S_mean)/S_std
        return S_line

    S_line = S_line_calc(R_line)

    r_linha_ast = -(R_std/(R_std**2 + S_std**2)**0.5)*beta_HL
    s_linha_ast = (S_std/(R_std**2 + S_std**2)**0.5)*beta_HL

    plt.plot(R_line, S_line, label = 'Limit State Line', color = 'black')
    plt.scatter(r_linha_ast, s_linha_ast, label = 'Design Point', color = 'red')
    plt.scatter(0, 0, color = 'blue', label = 'Checking Point')
    plt.legend()
    plt.show()

    p_f = 1 - npd.phi(beta_HL, 0, 1)

    print('\nBy using the analytical model for the Hasofer-Lind Method with the linear limit state equation: ').\
    print('Reliability index: {}\nProbability of Failure = {}%.\n'.format(beta_HL, p_f*100))
    


def afosm_nl_ls(g, X_1_mean, X_2_mean):
    d = 2
    R = np.eye(d)
    marginals = np.array(R, ndmin=1)
    d_2 = len(marginals)
    x_init = np.tile([0.1],[d_2,1]) #initial search point

    # objective function from the book
    dist_fun = lambda u: np.linalg.norm(u) 

    # method for minimization
    alg = 'SLSQP'

    H = lambda u: g(u)
    cons = ({'type': 'ineq', 'fun': lambda u: -H(u.reshape(-1,1))})
    
    #constraints 
    result = sp.optimize.minimize(dist_fun, x0 = x_init, constraints = cons, method=alg)

    #geting main results
    beta_value = result.fun
    iterations = result.nit

    p_f = 1 - npd.phi(beta_value, 0, 1)

    print('\nBy using the', alg, 'Method for minimizing the limit state equation from scipy package:')
    print('Iterations: {}\nReliability index = {}\nProbability of failure = {}%\n\n'.format(iterations, beta_value, p_f*100))

    return p_f, beta_value

g     = lambda x : 18*x[0,:] - 12*x[1,:] + 120 - 50
beta = afosm_nl_ls(g, 1, 1)
