import numpy as np
from scipy import stats as st
import modeling_of_uncertainty as mou
import matplotlib.pyplot as plt
import math as mt
import seaborn as sns

def norm_pdf_plot(sample): #Plots the normal pdf with normalized histogram
    mean = mou.mean_sample(sample)
    std = mou.std_sample(sample)
    step = abs(min(sample) - max(sample))/100
    x = np.arange(min(sample), max(sample) + step, step)
    pdf = st.norm.pdf(x, mean, std)
    plt.plot(x, pdf, color = 'black', label = 'Theoretical PDF')
    plt.legend()
    mou.hist(sample, 'sturges', dens_norm = True, cumulative_norm = False)
    
    return pdf

def norm_cdf_plot(sample, alpha): #Plots the normal theoretical cdf compared to the empirical one
    mean = mou.mean_sample(sample)
    std = mou.std_sample(sample)
    step = abs(min(sample) - max(sample))/100
    x = np.arange(min(sample), max(sample) + step, step)
    cdf = st.norm.cdf(x, mean, std)

    n = len(sample)
    y = np.arange(1, n+1)/n

    F1 = []
    F2 = []
    for i in range(0, len(sample)):
        e = (((mt.log(2/alpha))/(2*n))**0.5)  
        F1.append(y[i] - e)
        F2.append(y[i] + e)
  
    plt.plot(x, cdf, color = 'black', label = 'Theoretical CDF') 
    plt.scatter(sorted(sample), y, label='Empirical CDF')
    plt.plot(sorted(sample), F1, linestyle='--', color='red', alpha = 0.8, lw = 0.9, label = 'Dvoretzky–Kiefer–Wolfowitz Confidence Bands')
    plt.plot(sorted(sample), F2, linestyle='--', color='red', alpha = 0.8, lw = 0.9)
    plt.ylabel('Cumulative Distribution Function')
    plt.xlabel('Observed Data')
    plt.legend()
    plt.show()

    return cdf


def phi(x, mean, std): #Computes the probability of getting a value of x or lesser in a random variable (cdf) 
    phi = st.norm.cdf(x, mean, std)

    return phi

def mean_norm(sample, alpha):
    n = len(sample)
    std = mou.std_sample(sample)
    k = st.norm.ppf(1 - alpha/2)
    e = k*std/(n**0.5)

    print("The population mean with a confidence level of {}% is {} \u00B1 {}.".format(int((1-alpha)*100), mou.mean_sample(sample), e))

def phi_inverse(prob): #Computes the inverse of the normal for a given probability
    phi_inv = st.norm.ppf(prob)

    return phi_inv

def norm_qq_plot(sample, alpha): #plots the quantile-quantie plot for the given data
    y = np.arange(1, len(sample)+1)/(len(sample)+1)
    mean = mou.mean_sample(sample)
    std = mou.std_sample(sample)
    theo_qq = phi_inverse(y)
    x = theo_qq*std + mean

    #Kolmogorov-Smirnov Test for getting the confidence interval
    K = (-0.5*mt.log(alpha/2))**0.5
    M = (len(sample)**2/(2*len(sample)))**0.5
    CI_qq_high = []
    CI_qq_low = []
    for prob in y:
        F1 = prob - K/M
        F2 = prob + K/M
        s_low = phi_inverse(F1)
        s_high = phi_inverse(F2)
        CI_qq_low.append(s_low*std + mean)
        CI_qq_high.append(s_high*std + mean)
    sns.regplot(x, sorted(sample), ci = None, line_kws={'color':'black','label':'Regression Line'})
    plt.plot(sorted(sample), CI_qq_low, linestyle='--', color='red', alpha = 1, lw = 0.8, label = 'Kolmogorov-Smirnov Confidence Bands')
    plt.legend()
    plt.plot(sorted(sample), CI_qq_high, linestyle='--', color='red', alpha = 1, lw = 0.8)
    plt.xlabel('Theoretical Normal Quantiles')
    plt.ylabel('Sample Quantiles')

    plt.show()






sample = [28900, 29200, 27400, 28700, 28400, 29900, 30200, 29500, 29600, 28400, 28300, 29300, 29300, 28100, 30200, 30200, 30300, 31200, 28800, 27600, 29600, 25900, 32000, 33400, 30600, 32700, 31300, 30500, 31300, 29000, 29400, 28300, 30500, 31100, 29300, 27400, 29300, 29300, 31300, 27500, 29400]
