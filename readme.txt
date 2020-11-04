Description 
  Tuning paramameters selection for variable selection and estimation in  reproducing kernel Hilbert space. 
Usage
  generalGaussianKRR(Xtr,Ytr,fold,Xte)
Arguments
  Xtr, Xte   Xtr and Xte are the training and testing covariates
  Ytr           Ytr is the training reponse
  fold         fold > 1 is an interger with default value 5
Value
  A list with components:
  fit.value    The estimator with Xte
  selection   The index  of  selected variables
  gamma      The chosen tuning parameters of general Gaussian kernel
Examples
  require(MASS)
  source("generalGaussianKRR.R")
  n <- 100
  p <- 10
  library(MASS)
  set.seed(123)
  X <- mvrnorm(n,mu=rep(0,p),Sigma = diag(1,p))
  Y <- (2*X[,1]-1)*(2*X[,2]-1)+rnorm(n)
  generalGaussianKRR(X,Y)
