# matrixdist
Statistics for Matrix Distributions

This package implements tools which are useful for the statistical analysis of homogeneous and in-homogeneous phase--type distributions. These distributions are absorption times of Markov jump processes, and thus the maximization of their likelihood for statistical estimation is best dealt with using the EM algorithm. 

The main method for fitting in-homogeneous phase--type is *fit*, which allows data to be right-censored or aggregated. Several in-homogeneity transforms are considered, as well as the most common sub-intensity matrix structures. 

Simulation can be efficiently performed using the method *sim*. Various other methods for computing functionals of in-homogeneous phase--type objects are the following: *dens* (density), *cdf* (cumulative distribution function), *haz* (hazard rate), and *quan* (quantile). 

For phase--type distributions, the additional methods *minimum*, *maximum*, *"+"* obtain the matrix representations of the resulting phase--type distribution obtained through their application to two random variables. The method *moment* computes exact moments of phase--type distributions of any order. 

Finally, *show*, *coef*, *logLik* are auxiliary methods for visualizing, extracting coefficients and log-likelihood of phase--type objects, respectively.