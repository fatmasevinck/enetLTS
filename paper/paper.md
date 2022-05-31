---
title: '`enetLTS`: Robust and Sparse Methods for High Dimensional Linear and Binary and Multinomial Regression'
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


`enetLTS` is an R package [@R21] that provides a fully robust version of 
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


# Example: Robust and Sparse Linear Regression

We have considered the [NCI-60 cancer cell panel](https://discover.nci.nih.gov/cellminer/) data [@Reinhold12] in order to illustrate the functionality of the `enetLTS` model for linear regression. As in [@Alfons21] the response variable is determined by the protein expressions for a specific protein, which is 92th protein, and
the explanatory variable is determined by the gene expressions of the 100 genes that have the highest (robustly estimated) correlations with the response variable. This data set is available in package `robustHD`.

```R
> # load data
> library("robustHD")
> data("nci60")  # contains matrices 'protein' and 'gene'

> # define response variable
> y <- protein[, 92]
> # screen most correlated predictor variables
> correlations <- apply(gene, 2, corHuber, y)
> keep <- partialOrder(abs(correlations), 100, decreasing = TRUE)
> X <- gene[, keep]
```

Like many other packages, the easy way to use the package `enetLTS` is to install it directly from `CRAN`. 

```R
> # install and load package
> install.packages("enetLTS")
> library(enetLTS)
> # fit the model for family="gaussian"
> fit.gaussian <- enetLTS(X,y)
> [1] "optimal model: lambda = 0.1391 alpha = 0.6"
>
> fit.gaussian
enetLTS estimator 

Call:  enetLTS(xx = X, yy = y, family = "gaussian", alphas = alphas, 
 lambdas = lambdas, lambdaw = NULL, intercept = TRUE, scal = TRUE,
 hsize = 0.75, nsamp = c(500, 10), nCsteps = 20, nfold = 5,     
 repl = 1, ncores = 1, tol = -1e+06, seed = NULL, crit.plot = TRUE) 

Coefficients:
          1           2           3           4           5           6 
-2.52390658  0.31135509  0.00000000  0.00000000  0.12285091  0.00000000 
          7           8           9          10          11          12 
 0.00000000  0.00000000  0.00000000  0.07457828  0.00000000  0.00000000 
         13          14          15          16          17          18 
 0.09740240  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000 
         19          20          21          22          23          24 
 0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.05946501 
         25          26          27          28          29          30 
-0.27371935  0.00000000  0.00000000  0.00000000 -0.07528489  0.00000000 
         31          32          33          34          35          36 
 0.00000000  0.00000000  0.13028362  0.00000000  0.00000000  0.00000000 
         37          38          39          40          41          42 
 0.00000000  0.00000000  0.12902065 -0.02704781  0.00000000  0.00000000 
         43          44          45          46          47          48 
 0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000 
         49          50          51          52          53          54 
 0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000 
         55          56          57          58          59          60 
 0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000 
         61          62          63          64          65          66 
 0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000 
         67          68          69          70          71          72 
 0.06071918  0.00000000  0.10203048  0.00000000  0.00000000  0.00000000 
         73          74          75          76          77          78 
 0.09639756  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000 
         79          80          81          82          83          84 
 0.00000000  0.00000000 -0.03343054  0.00000000  0.00000000  0.00000000 
         85          86          87          88          89          90 
 0.00000000 -0.08272298  0.00000000 -0.09226326  0.00000000  0.00000000 
         91          92          93          94          95          96 
 0.00000000  0.17930052  0.00000000  0.00000000  0.09609275 -0.10894526 
         97          98          99         100         101 
 0.00000000  0.00000000  0.04865659  0.00000000  0.00000000 

 number of the nonzero coefficients:
[1] 21

 alpha: 0.725
 lambda: 0.1391
 lambdaw: 0.07936752
```

Several plots are available for the results: `plotCoef.enetLTS()` visualizes the coefficients, 
`plotResid.enetLTS()` plots the values of residuals vs fitted values, 
and `plotDiagnostic.enetLTS()` allows to produce various diagnostic
plots for the final model fit. 
Examples of these plots are shown in Figure \ref{fig:plotexamplesGuas}.

![Examples of plot functions of residuals (left); diagnostic (right)\label{fig:plotexamplesGuas}](JOSSgausNCI60.png)

# Example: Robust and Sparse Binary Regression




# Example: Robust and Sparse Multinomial Regression

The fuit data set has been well-known in the context of robust discrimination. 
It contains spectral information with 256 wavelengths,
thus is high-dimensional, for observations from 3 different cultivars of the same fruit, named
D, M, and HA, with group sizes 490, 106, and 500. 

```R
> # load data
> library(rrcov)
> data(fruit)
> 
> d <- fruit[,-1]  # first column includes the fruid names 
> xx <- as.matrix(d)
> # define response variable
> grp <- c(rep(1,490),rep(2,106),rep(3,500)) 
> yy <- factor(grp-1)
>
> fit.binomial <- enetLTS(xx, yy, family="multinomial",
                    alphas=seq(from=0.01,to=0.1,by=0.01), 
                    lambdas=seq(from=0.01,to=0.1,by=0.01),
                    lambdaw=NULL, intercept=TRUE, hsize=0.75, 
                    nsamp=c(500,10), nCsteps=20, nfold=5, repl=1, ncores=1, 
                    tol=-1e6, scal=TRUE, seed=NULL, crit.plot=TRUE)
 ```                


# Related Software

Package `robustHD` provides the sparseLTS estimator for linear regression based on the trimming idea for high dimensional linear regression [@Alfons21R]. Package `pense` provides implementations of robust S- and MM-type estimators using elastic net
regularization for linear regression [@Keplinger21R]. Package `glmnet` implements the elastic net estimator for generalized linear regression models [@Friedman21R]. Moreover, the procedure of the R package `enetLTS` [@Kurnaz22Rcran]. is implemented using internally the R package `glmnet` [@Friedman21R].  


# Acknowledgements

Fatma Sevinc KURNAZ is supported by grant TUBITAK 2219 from Scientific and Technological Research Council of Turkey (TUBITAK). 


# References
