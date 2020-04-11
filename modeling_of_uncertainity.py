import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as st
import math as mt

sns.set(style = 'whitegrid')

def mean_sample(sample): #First moment, compute the mean from the observed data
    mean = sum(sample)/len(sample)

    return mean

def std_sample(sample):  #Second moment, compute the standard deviation from the sample (degree of freedom = 1)
    mean = mean_sample(sample)
    array = []
    for data in sample:
        array.append((data - mean)**2)
    std_sample = ((1/(len(sample) - 1))*sum(array))**0.5

    return std_sample

def std_pop(sample):  #compute the standard deviation from the population (degree of freedom = 0)
    mean = mean_sample(sample)
    array = []
    for data in sample:
        array.append((data - mean)**2)
    std_pop = ((1/(len(sample)))*sum(array))**0.5

    return std_pop

def var_sample(sample): #compute the variance from the sample
    var_sample = std_sample(sample)**2

    return var_sample

def var_pop(sample): #compute the variance from the population
    var_pop = std_pop(sample)**2
    
    return var_pop

def cov_sample(sample):  #compute the coeficient of variation from the sample
    #Smaller COV -> smaller amount of uncertainity
    #in many engineering problems, COV is between 0.1 and 0.3
    cov = std_sample(sample)/mean_sample(sample)

    return cov

def cov_pop(sample): #compute the coeficient of variation from the population
    cov = std_pop(sample)/mean_sample(sample)

    return cov

def skewness(sample): #Third moment, compute the skewness from the sample
    #Skewness will say if the observed data is symmetrical or unsymmetrical,
    #if it's zero, then the data is symmetric. 
    mean = mean_sample(sample)
    array = []
    for data in sample:
        array.append((data - mean)**3)
    skw = (1/len(sample))*sum(array)

    return skw

def skewness_coef(sample): #Skewness coeficient is a nondimensial measure of the skewness to avoid dimensional problems
    #If the skewness coeficient is positive, the dispersion is more above the mean than below the mean,
    #if the skewness coeficient is negative, the dispersion is more below the mean than above the mean.
    skw_coef = skewness(sample)/(std_sample(sample)**3)

    return skw_coef

def kurtosis(sample): #The fourth moment, it says how much flatten is the probability density curve
    mean = mean_sample(sample)
    array = []
    for data in sample:
        array.append((data - mean)**4)    
    kt = (1/len(sample))*sum(array)

    return kt

def kurtosis_coef(sample): #Kurtosis coeficiente is also a nondimensional measure of the kurtosis to avoid dimensional problems
    kurtosis_coef = kurtosis(sample)/(std_sample(sample)**4)
    
    return kurtosis_coef

def hist(sample, bins, dens_norm, cumulative_norm): #Simple plot of the histogram from the observed data
    plt.hist(sample, bins, density = dens_norm, cumulative = cumulative_norm)
    plt.xlabel('Observed Data')
    plt.ylabel('Frequency')
    plt.show()

def empirical_cdf(sample, alpha): #Plot empirical cfd with confidence interval
    n = len(sample)
    y = np.arange(1, n+1)/n

    #Computing confidence interval with the Dvoretzky–Kiefer–Wolfowitz method based on the empirical points
    F1 = []
    F2 = []
    for i in range(0, len(sample)):
        e = (((mt.log(2/alpha))/(2*n))**0.5)  
        F1.append(y[i] - e)
        F2.append(y[i] + e) 
    plt.scatter(sorted(sample), y, label='Empirical CDF')
    plt.plot(sorted(sample), F1, linestyle='--', color='red', alpha = 0.8, lw = 0.9, label = 'Dvoretzky–Kiefer–Wolfowitz Confidence Bands')
    plt.plot(sorted(sample), F2, linestyle='--', color='red', alpha = 0.8, lw = 0.9)
    plt.ylabel('Empirical Cumulative Distribution Function')
    plt.xlabel('Observed Data')
    plt.legend()
    plt.show()

    return y

def mode_sample(sample):
    md = st.mode(sample)

    return md

def median(sample):
    med = np.median(sample)

    return med

sample = [28900, 29200, 27400, 28700, 28400, 29900, 30200, 29500, 29600, 28400, 28300, 29300, 29300, 28100, 30200, 30200, 30300, 31200, 28800, 27600, 29600, 25900, 32000, 33400, 30600, 32700, 31300, 30500, 31300, 29000, 29400, 28300, 30500, 31100, 29300, 27400, 29300, 29300, 31300, 27500, 29400]

empirical_cdf(sample, 0.05)

#print(mean_sample(sample),'\n',var_sample(sample),'\n',std_sample(sample),'\n',cov_sample(sample),'\n',skewness(sample),'\n',skewness_coef(sample))
#hist(sample, 'sturges', True, False)