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


`enetLTS` is an R ([Development Core Team, 2021](https://www.R-project.org/)) package that provides a fully robust version of 
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
regression [@Friedman10]. 
Selecting optimal model with optimal tuning parameters is done via cross-validation, 
and various plots are available to illustrate model selection and to evaluate the 
final model estimates. 


# State of the field

# Statement of need

A number of new robust linear regression methods have been developed during the last 
decade to improve the calculation for high dimensional linear regression, such as 
`[@Alfons21R:2021; @Keplinger21R:2021]`. 
However, to the best of our knowledge, the robust logistic (both binary and multinomial) 
regression for high dimensional data is not available in elsewhere.
Package `enetLTS` therefore provides
researchers with access to robust solutions and variable selection at the same time
with high-dimensional linear and logistic regression data. 
It has been used in many benchmarking studies in the statistical
literature e.g. `[@Insolia21a:2021; @Insolia21b:2021; Monti21:2021]`, 
as well as in empirical research e.g. `[@Segaert18:2018; @Jensch22:2022]`.



Some text.... 
Here is the sample reference [@Insolia21a]

[@Alfons21R]


# Statement of need 

Yet another section


# Installation and basic usage

Sample code 

```R
> install.packages("enesLTS")
```

Sample equation 

$$
y = \beta_0 + \beta_1 x + \varepsilon
$$

References will automatically be placed in the references part. 

# References
