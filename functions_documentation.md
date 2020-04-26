## Functions Documentations
This file presents a documentation of the functions presented in each module of the package.

## Modeling of Uncertainty
* **``mean_sample(x)``**: This will compute the mean of the observed data, where ``x`` is an 1-D array with the observed data.
* **``std_pop(x)``**: Using zero degree of freedom, this will compute the standard deviation of the population, where ``x`` is an 1-D array with the observed data. 
* **``std_sample(x)``**: With one degree of freedom, this will compute the standard deviation of the mean, where ``x`` is an 1-D array with the observed data.
* **``var_pop(x)``**: Using a zero degree of freedom, this will compute the variance for a population, where ``x`` is an 1-D array with the observed data.
* **``var_sample(x)``**: Using one degree of freedom, this will compute the variance for a sample, where ``x`` is an 1-D array with the observed data.
* **``cov_sample(x)``**: This will compute the coefficient of variation of the sample, where ``x`` is an 1-D array with the observed data. A small COV implies in a smaller amount of uncertainty. Often, in engineering, the COV is between 0.1 and 0.3.
* **``cov_pop(x)``**: This will compute the coefficient of variation of the population, where ``x`` is an 1-D array with the observed data.
* **``skewness(x)``**:  This is going to compute the skewness of the sample, where ``x`` is a 1-D array with the observed data. 
* **``skewness_coef(x)``**: This is going to compute the skewness coefficient in order to avoid dimensional problems with the skewness value, where ``x`` is an 1-D array with the observed data. A positive coefficient means that the dispersion is more above the mean than below the mean, while a negative one means the opposite. 
* **``kurtosis(x)``**: This is going to compute the kurtosis of the sample, where ``x`` is a 1-D array with the observed data. 
* **``kurtosis_coef(x)``**: his is going to compute the kurtosis coefficient in order to avoid dimensional problems with the kurtosis value, where ``x`` is an 1-D array with the observed data. It shows how much flatten is the curve of the density of probability. 
* **``hist(x, bins, dens_norm)``**: This will plot the histogram of the sample. ``x`` is 1-D array with the observed data, ``bins`` may be an ``int`` or one of the following methods to compute the number of bins of a histogram: 'sturges', 'doane', 'scott', 'fd' (Freedman-Diaconis estimator), 'stone', 'rice' and 'sqrt' (for more information, see [numpy.hist](https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html)) and ``dens_norm`` may be a boolean (True or False), if it’s true, it will normalize the counts to construct a density plot.
* **`` empirical_cdf(x, alpha)``**: This will plot empirical cdf with confidence interval based on the Dvoretzky–Kiefer–Wolfowitz method. ``x`` is an 1-D array with the observed data and ``alpha`` is the significance level. 

##Normal Probability Distribution

* **``norm_pdf_plot(x)``**: This will plot probability density function, where ``x`` is an 1-D array.
* **``norm_cdf_plot(x, alpha)``**: This will plot the normal theoretical  cumulative distribution function compared to the observed data, where ``x`` is an 1-D array with the observed data and ``alpha`` is the significance level. 
* **``phi(x, mean, std)* **``: This will compute the probability of getting a value of x or lesser in a random variable in a  cumulative distribution function. Here, ``x`` is a observation vector or any integer, ``mean`` is a given mean, and ``std`` stands for the standard deviation for the given mean.
* **``mean_norm(x, alpha)``**: This will return the mean of the observed data, where ``x`` is an 1-D array with the observed data and ``alpha`` is the significance level.
* **``phi_inverse(prob)``**: This will return the inverse of the normal for a given probability.
* **``norm_qq_plot(x, alpha)``**:This will plot the quantile-quantie for the observed data, where ``x`` is an 1-D array with the observed data and ``alpha`` is the significance level.

##Normal Statistical Tests
* **``chi_square_test(x, alpha))``**: This computes the chi square test for, where ``x`` is an 1-D array with the observed data and ``alpha`` is the significance level. 
* **``ks_test(x, alpha))``**: This computes the Kolmogorov-Smirnov test,  where ``x`` is an 1-D array with the observed data and ``alpha`` is the significance level.

##Lognormal Statistical Tests
* **``chi_square_test(x, alpha))``**: This computes the chi square test for, where ``x`` is an 1-D array with the observed data and ``alpha`` is the significance level. 
* **``ks_test(x, alpha))``**: This computes the Kolmogorov-Smirnov test,  where ``x`` is an 1-D array with the observed data and ``alpha`` is the significance level.

##Lognormal Probability Distribution
* **``lognorm_pdf_plot(x)``**: This will plot probability density function, where ``x`` is an 1-D array.
* **``lognorm_cdf_plot(x, alpha)``**: This will plot the normal theoretical  cumulative distribution function compared to the observed data, where ``x`` is an 1-D array with the observed data and ``alpha`` is the significance level. 
* **``phi(x, mean, s)* **``: This will compute the probability of getting a value of x or lesser in a random variable in a  cumulative distribution function. Here, ``x`` is a observation vector or any integer, ``mean`` is a given mean, and ``s`` stands for the standard deviation for the given mean.
* **``phi_inverse(prob, s, mean)``**: This will return the inverse of the normal for a given probability, where ``prob`` is the probability, ``s`` stands for the standard deviation for the given mean , and mean is the given mean.
* **``lognorm_qq_plot(x, alpha)``**:This will plot the quantile-quantie for the observed data, where ``x`` is an 1-D array with the observed data and ``alpha`` is the significance level.
