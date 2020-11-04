generalGaussianKernel<- function(D,gamma){
  p=length(gamma)
  Dk=0                        ##Dk_ij=sum(gamma*(x_j^k-x_i^2)^2)
  for (k in 1:p) {
    Dk <- Dk+gamma[k]^2*D[,,k]
  }
  kernel=exp(-Dk) 
}
fitKRR <- function(X,y,X_new,gamma)
{
  n=nrow(X)
  alambda= seq(1e-5,30,len=50)*log(n)/n
  
  n_new=nrow(X_new)
  p=ncol(X)
  
  D <- array(0,c(n,n,p))        
  for(k in 1:p){
    D[,,k]<- as.matrix(dist(X[,k],method = "euclidean")^2)
  }
  I=diag(n)
  D_new <- array(0,c(n_new,n,p))        
  for(k in 1:p){
    xx=matrix(rep(X[,k],n_new),nrow = n,ncol = n_new,byrow = FALSE)
    xx_new=matrix(rep(X_new[,k],n),nrow = n_new,ncol=n,byrow = FALSE)
    D_new[,,k]<- (xx_new-t(xx))^2
  }
  
  K=generalGaussianKernel(D,gamma)
  K_new=generalGaussianKernel(D_new,gamma)
  
  GCV=rep(NA,length(alambda))
  for (lamk in 1:length(alambda)) {
    S=K%*%solve(K+alambda[lamk]*I)
    GCV[lamk]=n*sum(((I-S)%*%y)^2)/((n-sum(diag(S)))^2)
  }
  opt.lambda=alambda[which.min(GCV)]
  alpha=solve(K+opt.lambda*I)%*%y
  
  pred=K_new%*%alpha
  
  return(pred)
}
choosegamma0TMP <- function(X,y,lambda=1/length(y),tau=1/length(y), 
                            gamma0=-1, lambda0 = 1/sqrt(n),    fold)
{
  n=nrow(X)
  p=ncol(X)
  
  D <- array(0,c(n,n,p))       
  for(k in 1:p){
    D[,,k]<- as.matrix(dist(X[,k],method = "euclidean")^2)
  }
  
  ni=floor(n/fold)
  Iv <- list()   # sub-samples 
  Im <- list()   # Identity matrix
  for (i in 1:(fold-1))
  {
    Iv[[i]] = setdiff(1:n,((i-1)*ni+1):(i*ni))
    Im[[i]]  = diag(1, ni)
    
  }
  Iv[[fold]] = setdiff(1:n,((fold-1)*ni+1):n)
  Im[[fold]] = diag(1, n-(fold-1)*ni)
  
  Ialpha <- list()   # alpha
  E <- list()
  Minv <- list()
  
  s_gamma <- function(gamma, lambda=1, tau=0){
    K <<- generalGaussianKernel(D,gamma)
    s_gammaV = tau*sum(abs(gamma))   #LASSO penalty
    for (i in 1:fold)
    {
      Minv[[i]] <<- solve(K[-Iv[[i]],-Iv[[i]]]+lambda*Im[[i]])
      Ialpha[[i]] <<- Minv[[i]]%*%y[-Iv[[i]]]
      E[[i]] <<- y[Iv[[i]]]-K[Iv[[i]],-Iv[[i]]]%*%Ialpha[[i]]
      s_gammaV= s_gammaV + (mean(E[[i]]^2))
    }
    return(s_gammaV)
  }
  
  gradient <- function(gamma, lambda=1, tau=0){
    gradientV = rep(0,p)
    
    for (k in Imin) {
      DK=D[,,k]
      gradientV[k] = tau*sign(gamma[k])  
      
      for (i in 1:fold)
      {
        Dk_1=DK[Iv[[i]],Iv[[i]]]
        Dk_2=DK[-Iv[[i]],-Iv[[i]]]
        Dk_12=DK[Iv[[i]],-Iv[[i]]]
        
        gradientV[k] =gradientV[k] + 
          4*gamma[k]*mean(E[[i]]*((Dk_12*K[Iv[[i]],-Iv[[i]]])%*%Ialpha[[i]] - 
                                    K[Iv[[i]],-Iv[[i]]]%*%Minv[[i]]%*%(Dk_2*K[-Iv[[i]],-Iv[[i]]])%*%Minv[[i]]%*%y[-Iv[[i]]]))
      }
    }
    return(gradientV)
  }
  
  f0 <- function(coef0, gamma_0, lambda=1)
  {
    gamma = gamma_0*coef0
    f = s_gamma(gamma, lambda=lambda, tau=0)
    return(f)
  }
  
  f1 <- function(coef0=lambda, gamma_0)
  {
    f = s_gamma(gamma_0, lambda = coef0, tau=0)
    return(f)
  }
  
  if  (max(gamma0) < 0)   
  {
    gamma_0 = rep(1/sqrt(p),p)
    lambda0 = 1/sqrt(n)
    
    for (iter in 1:20)
    {
      coef0 = optimize(f0, c(0.05, 5), tol=0.0001, gamma_0=gamma_0, lambda=lambda0)
      gamma_0 = gamma_0*coef0$minimum
      gamma0 = gamma_0
      
      coef0 = optimize(f1, c(0.01, 10), tol=0.0001, gamma_0=gamma_0)
      lambda_0 = coef0$minimum
      
      if (abs(lambda0-lambda_0) < 1.0e-5) break
      
      lambda0 = lambda_0
      
    }  
  }
  
  gamma_1 = gamma0
  
  m = 500
  for (k in 1:ceiling(p/m))
  {
    Imin  = ((k-1)*m+1):min((k*m),p)
    optgamma=optim(gamma_1,fn=s_gamma,gr=gradient,method="BFGS", lambda=lambda0, tau=tau)
    gamma_1 = abs(optgamma$par)
  }
  
  opt.gamma=abs(optgamma$par)
  opt.gamma = opt.gamma*(opt.gamma/max(opt.gamma) > 1.0e-3)
  coef1 = optimize(f0, c(0.05, 20), tol=0.0001, gamma_0=opt.gamma)
  opt.gamma = opt.gamma*coef1$minimum
  
  err = s_gamma(opt.gamma, lambda = lambda, tau=0)
  criterion = log(err) + log(n)*sum(abs(opt.gamma)>0)/n
  
  
  return(list(gamma= abs(opt.gamma),gamma0 =gamma0, 
              criterion=criterion, lambda0 = lambda0 ))
}

chooseGamma = function(X, y, fold)
{
  n=nrow(X)
  p=ncol(X)
  
  stdX = apply(X, 2, sd)+1.0e-10
  meanX = colMeans(X)
  X = (X-matrix(meanX, n, p, byrow = TRUE))/matrix(stdX, n, p, byrow = TRUE)
  
  stdy=sd(y)
  meany=mean(y)
  y=(y-meany)/stdy
  
  criterion0 = 1.0e10
  criterionALL = c()
  
  lambda0 = -1
  gamma0  = -1
  
  for (tau in (1:10/10)^2*10)
  {
    A = choosegamma0TMP(X, y, tau=tau, lambda0 = lambda0, gamma0=gamma0, fold=fold)
    lambda0 = A$lambda0 
    gamma0 =  A$gamma0
    criterionALL = c(criterionALL, A$criterion)
    if (A$criterion < criterion0)
    {
      BIC0 = A$criterion
      A0 = A
    }
  }
  
#  print(criterionALL)
  
  return(list(gamma= A0$gamma, gamma0 = A0$gamma0,lambda = A0$lambda,
              criterion=A0$criterion, criterionALL=criterionALL))
}

generalGaussianKRR <- function(X,y,fold=5,X_new=X){
  p=ncol(X)
  n=nrow(X)
  n_new=nrow(X_new)
  
  stdX = apply(X, 2, sd)+1.0e-10
  meanX = colMeans(X)
  zX = (X-matrix(meanX, n, p, byrow = TRUE))/matrix(stdX, n, p, byrow = TRUE)
  zX_new=((X_new-matrix(meanX, n_new, p, byrow = TRUE))/matrix(stdX, n_new, p, byrow = TRUE))
  
  stdy=sd(y)
  meany=mean(y)
  zy=(y-meany)/stdy
  
  choose=chooseGamma(zX,zy,fold)
  optgamma=choose$gamma
  fit=fitKRR(zX,zy,zX_new,optgamma)
  
  selectindex=which(optgamma!=rep(0,p))
  
  return(list(fit.value=fit,selection=selectindex,gamma=optgamma))
}