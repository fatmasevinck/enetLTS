# enetLTS: Robust and Sparse Methods for High Dimensional Linear and Binary and Multinomial Regression

## Summary

`enetLTS` is an `R` package that provides a fully robust version of 
elastic net estimator for high dimensional linear and logistic (including 
binary and multinomial) regression. The elastic net penalization provides 
intrinsic variable selection and coefficient estimates for highly correlated 
variables in particular for high-dimensional low sample size 
data sets, and it has been extended to generalized linear regression models 
([Friedman et al. 2010](https://www.jstatsoft.org/article/download/v033i01/361)). 
Combining these advantages with trimming idea yields the robust solutions.
The main idea of the algorithm is to search for outlier-free subsets on which the classical elastic 
net estimators can be applied. Outlier-free subsets are determined by trimming 
the penalized log-likelihood function belonging to the regression model. 
The algorithm starts with 500 elemental subsets
only for one combination of $\alpha$ and $\lambda$, and takes the *warm start* strategy
for subsequent combinations in order to save the computation time.
The final reweighting step is added to improve the statistical 
efficiency of the proposed methods. 
From this point of view, the enet-LTS estimator can be seen as trimmed version 
of the elastic net regression estimator for linear, binary and multinomial 
regression ([Friedman et al. 2010](https://www.jstatsoft.org/article/download/v033i01/361)). 
Selecting optimal model with optimal tuning parameters is done via cross-validation, 
and various plots are available to illustrate model selection and to evaluate the 
final model estimates. 

## Implemented Methods 

- `enetLTS()`: elastic net trimmed squared regression for families:

   1- `gaussian`

   2- `binomial`
   
   3- `multinomial`
                                                                  

## Installation

Package `enetLTS` is on CRAN (The Comprehensive R Archive Network), hence the latest release can be easily installed from the R command as follows

```R
> install.packages("enetLTS")
```

## Building from source

To install the latest (possibly unstable) version from GitHub, you can pull this repository and install it from the R command line as follows

```R
> install.packages("devtools")
> devtools::install_github("fatmasevinck/enetLTS")
```

If you already have package devtools installed, the first line can be skipped.


# Example: Robust and Sparse Linear Regression

Like many other packages, the easy way to use the package `enetLTS` is to install it directly from `CRAN`. 

```R
> # install and load package
> install.packages("enetLTS")
> library(enetLTS)
> # fit the model for family="gaussian"
> fit.gaussian <- enetLTS(X,y)
```

Several plots are available for the results: plotCoef.enetLTS() visualizes the coefficients, 
plotResid.enetLTS() plots the values of residuals vs fitted values, 
and plotDiagnostic.enetLTS() allows to produce various diagnostic
plots for the final model fit. 
Examples of these plots are shown in Figure \ref{fig:plotexamples}.

# Example: Robust and Sparse Binary Regression 

# Example: Robust and Sparse Multinomial Regression

# References 

Friedman J., Hastie T. and Tibshirani R. (2010) Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software, 33(1), 1-22. DOI
[10.1163/ej.9789004178922.i-328.7](https://www.jstatsoft.org/article/download/v033i01/361)

