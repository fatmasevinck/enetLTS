---
title: '`enetLTS`: Robust and Sparse Methods for High Dimensional Linear and binary and Multinomial Regression'
tags:
  - R
  - Robust regression
  - Elastic net
  - outlier detection
authors:
  - name: Fatma Sevinc KURNAZ
    orcid: 0000-0002-5958-7366 
    affiliation: 1
  - name: Peter FILZMOSER
    orcid: 0000-0002-8014-4682
    affiliation: 2
affiliations:
 - name: Department of Statistics, Yildiz Technical University, Istanbul, Turkey
   index: 1
 - name: Institute of Statistics and Mathematical Methods in Economics, TU Wien, Vienna, Austria
   index: 2



date: 24 May 2022
bibliography: paper.bib
---

# Summary


`enetLTS` is an R [@R21] package that provides a fully robust version of 
elastic net estimator for high dimensional linear and logistic (including 
binary and multinomial) regression. The elastic net penalization provides 
intrinsic variable selection and coefficient estimates for highly correlated 
variables in particular for high-dimensional low sample size 
data sets, and it has been extended to generalized linear regression models 
[@Friedman10]. 
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
regression [@Friedman10]. 
Selecting optimal model with optimal tuning parameters is done via cross-validation, 
and various plots are available to illustrate model selection and to evaluate the 
final model estimates. 


# State of the field

# Statement of need

A number of new robust linear regression methods have been developed during the last 
decade to improve the calculation for high dimensional linear regression, such as 
[@Alfons21R; @Keplinger21R]. 
However, to the best of our knowledge, the robust logistic (both binary and multinomial) 
regression for high dimensional data is not available in elsewhere.
Package `enetLTS` therefore provides
researchers with access to robust solutions and variable selection at the same time
with high-dimensional linear and logistic regression data. 
It has been used in many benchmarking studies in the statistical
literature e.g. [@Insolia21b; @Insolia21a; @Monti21], 
as well as in empirical research e.g. [@Segaert18; @Jensch22].


# Installation and basic usage

Sample code 

```R
> install.packages("enesLTS")
```


# Mathematics

Package `enetLTS` fits the linear regression, binary and multinomial regression via penalized maximum likelihood solving the following problem:

\begin{equation}
\label{eq:mlproblem}
\operatorname*{arg\,min}_{\pmb{\beta}}
\left\{\sum_{i=1}^h \ell_i(y_i,\beta_0+x_i^T\pmb{\beta}) + h\lambda\left[(1-\alpha)\frac{1}{2} \sum_{j=1}^p\pmb{\beta}_j^2 + \alpha \sum_{j=1}^p|\pmb{\beta}_j|\right]
\right\}
\end{equation}

In the equation \autoref{eq:mlproblem}, $\ell_i(y_i,\beta_0+x_i^T\pmb{\beta})$ corresponds to the negative log-likelihood contribution belonging to the `family` for observation $i$. The non-negative tuning parameter $\lambda$ controls the entire strength of the
penalty. The tuning parameter $\alpha \in [0,1]$ allows to mix the
proportion of the ridge ($L_2$) and the lasso ($L_1$) penalty.

# Example: Robust and Sparse Linear Regression

Like many other packages, the esy way to use the package `enetLTS` is to install it directly from `CRAN`. 

```{r, echo = FALSE, eval = FALSE}
# install and load package
install.packages("enetLTS")
library(enetLTS)
# fit the model for family="gaussian"
fit.gaussian <- enetLTS(X,y)
```

Several plots are available for the results: plotCoef.enetLTS() visualizes the coefficients, 
plotResid.enetLTS() plots the values of residuals vs fitted values, 
and plotDiagnostic.enetLTS() allows to produce various diagnostic
plots for the final model fit. 
Examples of these plots are shown in Figure \ref{fig:plotexamples}.


# Example: Robust and Sparse Binary Regression




# Example: Robust and Sparse Multinomial Regression



# References
