import numpy as np 
import lognormal_probability_distribution as lnpd 
import modeling_of_uncertainty as mou
import math as mt
from scipy import stats as st
from scipy.stats import chi2

def chi_square_test(sample, alpha): 
    m =  mt.ceil(1 + 3.322*mt.log(len(sample),10))
    k = 2 
    f = m - 1 - k
    c = chi2.ppf(1-alpha, f)

    h = (max(sample)-min(sample))/m

    ln_x = []
    for data in sample:
        ln_x.append(mt.log(data))
    mean = mou.mean_sample(sample)
    std = mou.std_sample(ln_x)
  
    class_ = [min(sample)]
    for i in range(0, m):
        class_.append(class_[i] + h)

    e = []
    for i in range(1, len(class_)):
        if i == 1:
            e.append(lnpd.phi(class_[i], mean, std)*len(sample))
        elif i == (len(class_)):
            e.append((1 - lnpd.phi(class_[i-1], mean, std))*len(sample))
        else:
            e.append((lnpd.phi(class_[i], mean, std)-lnpd.phi(class_[i-1], mean, std))*len(sample))

    n = []
    t = 0
    i = 1
    sample = sorted(sample)
    sample.append(max(sample)+1)

    for j in range(0, len(sample)-1):
        if sample[j+1] > class_[i]:
            n.append((j+1)-t)
            t = j + 1
            i = i + 1

    test_1 = []
    for i in range(0, len(n)):
        test_1.append(((e[i]-n[i])**2)/e[i])
    
    statistics = sum(test_1)

    if statistics < c:
        print('Chi-Square Statistics = {}\nP-Value = {}\nConsidering that Chi-Square Statistics < P-Value, with a confidence level of {}%, the lognormal distribution is acceptable.'.format(statistics, c, (1-alpha)*100))
    
    elif statistics > c:
        print('Chi-Square Statistics = {}\nP-Value = {}\nConsidering that Chi-Square Statistics > P-Value, with a confidence level of {}%, the lognormal distribution is not acceptable.'.format(statistics, c, (1-alpha)*100))


def ks_test(sample, alpha):
    sort_sample = sorted(sample)
    n = len(sort_sample)
    S = np.arange(1, n + 1)/n
    ln_x = []
    for data in sample:
        ln_x.append(mt.log(data))
    mean = mou.mean_sample(sample)
    std = mou.std_sample(ln_x)
    
    F = lnpd.phi(sort_sample, mean, std) 
    D = max(abs(F-S))
    D_ks = st.ksone.ppf(1 - alpha/2, n)

    if D < D_ks:
        print('Kolmogorov-Smirnov Statistics = {}\nP-Value = {}\nConsidering that K-S Test Statistics < P-Value, with a confidence level of {}%, the normal distribution is acceptable.'.format(D, D_ks, (1-alpha)*100))
    
    elif D > D_ks:
        print('Kolmogorov-Smirnov Statistics = {}\nP-Value = {}\nConsidering that K-S Test Statistics > P-Value, with a confidence level of {}%, the normal distribution is not acceptable.'.format(D, D_ks, (1-alpha)*100))

