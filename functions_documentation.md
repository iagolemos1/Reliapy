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
* **``hist(x, bins, dens_norm)``**: This will plot the histogram of the sample. ``x`` is 1-D array with the observed data, ``bins`` may be an ``int`` or one of the following methods to compute the number of bins of a histogram: 'sturges', 'doane', 'scott', 'fd' (Freedman-Diaconis estimator), 'stone', 'rice' and 'sqrt' (for more information, see [numpy.hist](https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html)) and ``dens_norm`` may be a True or False, if it’s true, it will normalize the counts to construct a density plot.
* **`` empirical_cdf(x, alpha)``**: This will plot empirical cdf with confidence interval based on the Dvoretzky–Kiefer–Wolfowitz method. ``x`` is an 1-D array with the observed data and ``alpha`` is the significance level. 

