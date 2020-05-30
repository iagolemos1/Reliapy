import numpy as np
from scipy import stats as st
import modeling_of_uncertainty as mou
import matplotlib.pyplot as plt
import math as mt
import seaborn as sns
from scipy import stats as st
sns.set_style("ticks")

def lognorm_pdf_plot(sample):
    ln_x = []
    for data in sample:
        ln_x.append(mt.log(data))
    mean = mou.mean_sample(sample)
    std = mou.std_sample(ln_x)
    step = abs(min(sample) - max(sample))/100
    x = np.arange(min(sample), max(sample) + step, step)
    pdf = st.lognorm.pdf(x, s = std, scale = mean, loc = 0)
    plt.plot(x, pdf, color = 'dimgray', label = 'Theoretical PDF')
    plt.legend()
    mou.hist(sample, 'sturges', dens_norm = True)
    print(pdf)
    
    return pdf

def lognorm_cdf_plot(sample, alpha): #Plots the normal theoretical cdf compared to the empirical one
    ln_x = []
    for data in sample:
        ln_x.append(mt.log(data))
    mean = mou.mean_sample(sample)
    std = mou.std_sample(ln_x)
    step = abs(min(sample) - max(sample))/100
    x = np.arange(min(sample), max(sample) + step, step)
    cdf = st.lognorm.cdf(x, s = std, scale = mean, loc = 0)

    n = len(sample)
    y = np.arange(1, n+1)/n

    F1 = []
    F2 = []
    for i in range(0, len(sample)):
        e = (((mt.log(2/alpha))/(2*n))**0.5)  
        F1.append(y[i] - e)
        F2.append(y[i] + e)
  
    plt.plot(x, cdf, color = 'dimgray', label = 'Theoretical CDF') 
    plt.scatter(sorted(sample), y, label='Empirical CDF')
    plt.plot(sorted(sample), F1, linestyle='--', color='red', alpha = 0.8, lw = 0.9, label = 'Dvoretzky–Kiefer–Wolfowitz Confidence Bands')
    plt.plot(sorted(sample), F2, linestyle='--', color='red', alpha = 0.8, lw = 0.9)
    plt.ylabel('Cumulative Distribution Function')
    plt.xlabel('Observed Data')
    plt.legend()
    plt.show()

    return cdf
    
def lognorm_cdf(x, mean, std): #Computes the probability of getting a value of x or lesser in a random variable (cdf) 
    delta = std/mean
    xi    = mt.log(1 + delta**2)**0.5
    lambda_est = mt.log(mean) - 0.5*(xi**2)

    cdf = st.norm.cdf(mt.log(x), lambda_est, xi)

    return cdf

def lognorm_ppf(prob, mean, std): #Computes the inverse of the normal for a given probability
    delta = std/mean
    xi    = mt.log(1 + delta**2)**0.5
    lambda_est = mt.log(mean) - 0.5*(xi**2)
    ppf = st.norm.ppf(prob, lambda_est, xi)

    return ppf

def lognorm_qq_plot(sample, alpha): #plots the quantile-quantie plot for the given data
    y = np.arange(1, len(sample)+1)/(len(sample)+1)
    ln_x = []
    for data in sample:
        ln_x.append(mt.log(data))
    mean = mou.mean_sample(sample)
    std = mou.std_sample(ln_x)
    theo_qq = st.lognorm.ppf(y, s = std, loc = 0, scale = mean)
    x = theo_qq

    #Kolmogorov-Smirnov Test for getting the confidence interval
    K = (-0.5*mt.log(alpha/2))**0.5
    M = (len(sample)**2/(2*len(sample)))**0.5
    CI_qq_high = []
    CI_qq_low = []
    for prob in y:
        F1 = prob - K/M
        F2 = prob + K/M
        s_low = st.lognorm.ppf(F1, s = std, loc = 0, scale = mean)
        s_high = st.lognorm.ppf(F2, s = std, loc = 0, scale = mean)
        CI_qq_low.append(s_low)
        CI_qq_high.append(s_high)
    sns.regplot(x, sorted(sample), ci = None, line_kws={'color':'dimgray','label':'Regression Line'})
    plt.plot(sorted(sample), CI_qq_low, linestyle='--', color='red', alpha = 1, lw = 0.8, label = 'Kolmogorov-Smirnov Confidence Bands')
    plt.legend()
    plt.plot(sorted(sample), CI_qq_high, linestyle='--', color='red', alpha = 1, lw = 0.8)
    plt.xlabel('Theoretical Lognormal Quantiles')
    plt.ylabel('Sample Quantiles')

    plt.show()
