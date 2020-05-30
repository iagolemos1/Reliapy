import numpy as np
import modeling_of_uncertainty as moe
import normal_probability_distribution as npd
import lognormal_probability_distribution as lpd
from scipy import stats as st 
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib import colors
import math as mt
from scipy import stats as st
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from sympy.solvers import solve
from sympy import Symbol


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


def single_load_lognorm_beta_value(R_mean, S_mean, R_std, S_std): 
    delta_R = R_std/R_mean
    delta_S = S_std/S_mean
    beta = (mt.log(R_mean/S_mean * ((1 + delta_S**2)/(1 + delta_R**2))**0.5))/(mt.log((1 + delta_S**2)*(1 + delta_R**2)))**0.5

    return beta

def p_failure_lognorm_single_load(R_mean, S_mean, R_std, S_std):
    beta = single_load_lognorm_beta_value(R_mean, S_mean, R_std, S_std)
    p_f = 1 - npd.phi(beta, 0, 1)

    print('The probability of failure is {}% with beta-value of {}.'.format(p_f*100, beta))

    return p_f, beta

def single_load_lognorm_factors(R_mean, S_mean, R_std, S_std, k_R, k_S):
    beta = single_load_lognorm_beta_value(R_mean, S_mean, R_std, S_std)
    delta_R = R_std/R_mean
    delta_S = S_std/S_mean 

    e_l = (mt.log(1 + delta_R**2) + mt.log(1 + delta_S**2))**0.5/((mt.log(1 + delta_R**2))**0.5 + (mt.log(1 + delta_S**2))**0.5)
    phi_linha = np.exp(-beta*e_l*(mt.log(1 + delta_R**2))**0.5)/(1 + delta_R**2)**0.5
    gamma_linha = np.exp(-beta*e_l*(mt.log(1 + delta_S**2))**0.5)/(1 + delta_S**2)**0.5

    phi_factor = phi_linha*np.exp(k_R*delta_R)
    gamma_factor = gamma_linha*np.exp(k_S*delta_S)

    print('\n\nNominal Capacity Reduction Factor = {}\nLoad Factor is {}\n\n'.format(phi_factor, gamma_factor))

    return(phi_factor, gamma_factor)

def multiple_load_lognorm_beta(R_mean, S_mean_array, R_std, S_std_array):
    S_mean_total = sum(S_mean_array)
    S_std_total = np.linalg.norm(S_std_array)
    delta_R = R_std/R_mean
    delta_S = S_std_total/S_mean_total
    beta = (mt.log(R_mean/S_mean_total * ((1 + delta_S**2)/(1 + delta_R**2))**0.5))/(mt.log((1 + delta_S**2)*(1 + delta_R**2)))**0.5

    return beta

def p_failure_norm_multiple_load(R_mean, S_mean_array, R_std, S_std_array):
    beta = multiple_load_norm_beta(R_mean, S_mean_array, R_std, S_std_array)
    p_f = 1 - npd.phi(beta, 0, 1)
    print('\n\nProbability of failure = {}%\nBeta-value = {}\n\n.'.format(p_f*100, beta))

    return p_f, beta

def p_failure_lognorm_multiple_load(R_mean, S_mean_array, R_std, S_std_array):
    beta = multiple_load_lognorm_beta(R_mean, S_mean_array, R_std, S_std_array) 
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

def multiple_load_lognorm_factors(R_mean, S_mean_array, R_std, S_std_array, k_R, k_S_array):
    from scipy.optimize import fsolve

    beta = multiple_load_lognorm_beta(R_mean, S_mean_array, R_std, S_std_array)

    S_mean_total = sum(S_mean_array)
    S_std_total = np.linalg.norm(S_std_array)
    delta_R = R_std/R_mean
    delta_S = S_std_total/S_mean_total 

    delta_S_array = []
    for i in range(0, len(S_std_array)):
        delta_S_array.append(S_std_array[i]/S_mean_array[i])
    
    e_l = (mt.log(1 + delta_R**2) + mt.log(1 + delta_S**2))**0.5/((mt.log(1 + delta_R**2))**0.5 + (mt.log(1 + delta_S**2))**0.5)

    eq1 = S_mean_total*np.exp(beta*e_l*(mt.log(1 + delta_S**2))**0.5)/(1 + delta_S**2)**0.5

    func = lambda e_n_l: eq1 - sum(S_mean_array[i]*mt.e**(beta*e_l*e_n_l*(mt.log(1 + delta_S_array[i]**2))**0.5)/(1 + delta_S_array[i]**2)**0.5 for i in range(0, len(delta_S_array)))

    e_n_l_initial_guess = e_l
    e_n_l_solution = fsolve(func, e_n_l_initial_guess)[0]

    phi_linha = np.exp(-beta*e_l*(mt.log(1 + delta_R**2))**0.5)/(1 + delta_R**2)**0.5
    
    gamma_linha_array = []
    for i in range(0, len(S_std_array)):
        gamma_linha_array.append(np.exp(beta*e_l*e_n_l_solution*(mt.log(1 + delta_S_array[i]**2))**0.5)/(1 + delta_S_array[i]**2)**0.5)

    phi_factor = phi_linha*np.exp(k_R*delta_R)

    gamma_factor_array=[]
    for i in range(0, len(S_std_array)):
        gamma_factor_array.append(gamma_linha_array[i]*np.exp(-k_S_array[i]*delta_S_array[i]))

    print('\n\nNominal Reduction Factor = {}.\nLoad Factors = {}\n\n'.format(phi_factor, gamma_factor_array))
    
    return (phi_factor, gamma_factor_array)
  

def afosm_l_ls(R_mean, S_mean, R_std, S_std): #reduced hasofer-lind method for linear limit state equation

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

    print('\nBy using the analytical model for the Hasofer-Lind Method with the linear limit state equation: ')
    print('Reliability index: {}\nProbability of Failure = {}%.\n'.format(beta_HL, p_f*100))
    
def afosm(g, d_input):
    d = d_input
    R = np.eye(d)
    marginals = np.array(R, ndmin=1)
    d_2 = len(marginals)
    x_init = np.tile([0.1],[d_2,1]) #initial search point

    # objective function from the book
    dist_fun = lambda u: np.linalg.norm(u) 

    # method for minimization
    alg = 'SLSQP'

    H = lambda u: g(u)
    cons = ({'type': 'eq', 'fun': lambda u: -(H(u.reshape(-1,1)))})
    
    #constraints 
    result = sp.optimize.minimize(dist_fun, x0 = x_init, constraints = cons, method=alg)

    #geting main results
    beta_value = result.fun
    iterations = result.nit
    u_ast = result.x

    if d == 2:
        plt.figure()
        xx      = np.linspace(-10, 10, 200)
        [X1,X2] = np.meshgrid(xx,xx)
        xnod    = np.array([X1,X2])
        ZX      = g(xnod)
        plt.pcolor(X1,X2,ZX, cmap = 'RdBu')
        plt.contour(X1,X2,ZX,[0]) 
        plt.scatter(0,0, color = 'dimgray', label = 'Standard space origin')     
        plt.plot([0, u_ast[0]],[0, u_ast[1]],label = 'SLSQP solver')    # reliability index beta
        plt.scatter(u_ast[0], u_ast[1], color = 'black', label = 'Design Point')                       # design point in standard    
        plt.title('Standard space')
        plt.xlabel(r"$X_1^{*}$'")
        plt.ylabel(r"$X_2^{*}$'")
        plt.legend()
        plt.show()


    p_f = 1 - npd.phi(beta_value, 0, 1)
    print('------------------------')
    print('\nSolver:'algo)
    print('Iterations: {}\nReliability index = {}\nProbability of failure = {}%\n\n'.format(iterations, beta_value, p_f*100))
    print('------------------------')

    return p_f, beta_value, u_ast

def lognorm2norm(x, mean, std):
    delta = std/mean
    xi    = mt.log(1 + delta**2)**0.5
    lambda_est = mt.log(mean) - 0.5*(xi**2)

    lognorm_pdf = (1/((2*mt.pi)**0.5 * xi * x))*np.exp(-0.5*((mt.log(x) - lambda_est)/xi)**2)
    variable_lognorm = (mt.log(x) - lambda_est)/xi

    pdf_normal = st.norm.pdf(variable_lognorm)

    std_normal = pdf_normal/lognorm_pdf
    mean_normal = x - variable_lognorm*std_normal

    print('\n--- Variables transformation - Rackawtiz-Fiessler method ---')
    print('\nVariable value = {}\n'.format(x))
    print('- Transformation from lognormal to normal -\n\n--- Lognormal parameters ---\nLognormal mean = {}\nLognormal standard deviation = {}'.format(mean, std))
    print('\n--- Normal parameters ---\nNormal mean = {}\nNormal standard deviation = {}'.format(mean_normal, std_normal))
    print('---------------------------------------------------------------\n')

    return (mean_normal, std_normal)


def form_I(g, ln_mean, ln_std, n_mean, n_std):
    alg = 'Rackwitz-Fiessler-Ayyub-Haldar'
    # -- Defining some important functions -- #
    def lognorm2norm2(x, mean, std): #Transformation from log to norm without prints
        delta = std/mean
        xi    = mt.log(1 + delta**2)**0.5
        lambda_est = mt.log(mean) - 0.5*(xi**2)

        lognorm_pdf = (1/((2*mt.pi)**0.5 * xi * x))*np.exp(-0.5*((mt.log(x) - lambda_est)/xi)**2)
        variable_lognorm = (mt.log(x) - lambda_est)/xi

        pdf_normal = st.norm.pdf(variable_lognorm)

        std_normal = pdf_normal/lognorm_pdf
        mean_normal = x - variable_lognorm*std_normal
        return (mean_normal, std_normal)
    
    from scipy.misc import derivative
    def partial_derivative(func, var=0, point=[]): #partial derivate
        args = point[:]
        def wraps(x):
            args[var] = x
            return func(*args)
        return derivative(wraps, point[var], dx = 1e-6)

    variables_mean = ln_mean + n_mean
    variables_list = []
    for i in range(0, len(variables_mean)):
        variables_list.append(i)

    beta_init = 3
    new_beta = beta_init
    beta = 0
    iterations = 0
    while abs(new_beta - beta) > 0.001:
        beta = new_beta
        cosines_2 = np.ones(len(variables_list))
        cosines_1 = np.zeros(len(variables_list))
        new_search_points = ln_mean + n_mean

        while (abs(cosines_2[0] - cosines_1[0]) >= 0.005):
            iterations = iterations + 1
            cosines_1 = cosines_2
            initial_search_points = new_search_points
            normal_mean_trans = []
            normal_std_trans = []

            if ln_mean != []:
                for i in range(0, len(ln_mean)):
                    [mean_normal, std_normal] = lognorm2norm2(initial_search_points[i], ln_mean[i], ln_std[i])
                    normal_mean_trans.append(mean_normal)
                    normal_std_trans.append(std_normal)
                std = normal_std_trans + n_std 
                mean = normal_mean_trans + n_mean
            else:
                std = n_std
                mean = n_mean 

            diff = []
            for i in range(0, len(initial_search_points)):
                diff.append(partial_derivative(g, variables_list[i], initial_search_points))
        

            g_std = []
            g_std_square = []
            for i in range(0, len(initial_search_points)):
                g_std.append(diff[i]*std[i])
                g_std_square.append((diff[i]*std[i])**2)

            cosines = []    
            for i in range(0, len(initial_search_points)):
                cosines.append(g_std[i]/(sum(g_std_square))**0.5)
            cosines_2 = cosines

            new_search_points = []
            for i in range(0, len(initial_search_points)):
                new_search_points.append(mean[i] - cosines[i]*new_beta*std[i])
    
        if ln_mean != []:
            normal_mean_trans_test = []
            normal_std_trans_test = []
            for i in range(0, len(ln_mean)):
                [mean_normal_test, std_normal_test] = lognorm2norm2(initial_search_points[i], ln_mean[i], ln_std[i])
                normal_mean_trans_test.append(mean_normal_test)
                normal_std_trans_test.append(std_normal_test)
            std_test = normal_std_trans_test + n_std 
            mean_test = normal_mean_trans_test + n_mean 
        else: 
            std_test = std
            mean_test = mean  
     
        coordinates = []
        B = Symbol('B')
        for i in range(0, len(variables_list)):
            coordinates.append(mean_test[i] - cosines[i]*B*std_test[i]) 
        new_beta = solve(g(*coordinates), B)[0]
    
    p_f = 1 - st.norm.cdf(float(new_beta))

    print('-------------------------')
    print('Algorithm: FORM method I -', alg, 'solver')
    print('Iterations: {}\nReliability index = {}\nProbability of failure = {}%'.format(iterations - 1, new_beta, p_f*100))
    print('-------------------------')
