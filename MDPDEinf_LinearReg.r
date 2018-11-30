## This function 'lmdpd' compute the MDPDE of (beta, sigma) and test for H0: beta_i = beta_i0 
## under a Linear regression Model y=Xbeta + eps, eps ~ N(0, sigma^2)
##
## References:
##  1. Ghosh, A., & Basu, A. (2013). Robust estimation for independent non-homogeneous observations using density power divergence with applications to linear regression. Electronic Journal of statistics, 7, 2420-2456.
##  2. Basu, A., Ghosh, A., Martin, N., & Pardo, L. (2018). Robust Wald-type tests for non-homogeneous observations based on the minimum density power divergence estimator. Metrika, 81(5), 493-522.
##
## Inputs: 
##  y       = Response vector   [n X 1] 
##  X       = Covariate Matrix  [n X p]
##  alpha   = DPD tuning paramter
##  p       = number of covariates (optional)
##  Initial = Initial values of the parameter for estimation process. Default: 1 for all parameters (initial of sigma needs to be given via log(sigma))
##  beta0   = Value of beta under null to be tested componenetwise> Default: all zero (significance of variable testing)
##
## Outputs: a list containing the folwloing 4 elements. 
## (1) Paramter estimates of beta and sigma [(p+1) X 1]
## (2) Asymptotic variance estimates of the MDPDEs in (1) [(p+1) X 1]
## (3) P-values of the test for individual component of beta=beta0 [p X 1]
## (4) An indicator if R-iterations in the estimation process has converged: 0 == convergence


lmdpd <- function(y,X,alpha,p=ncol(X),Initial=matrix(c(rep(1,p),0)),beta0=matrix(rep(1,p))) {
  
  n <- length(y)

#Objective function for the computation of MDPDE (Minimum DPD estimators)
objfunction_Linear_reg <- function(t) {
    beta<-t[1:p]  
    sigma<-exp(t[p+1])
    r<- (-((y-X%*%beta)^2)/(2*sigma*sigma))
    c<-1/(sigma^alpha)
  
    if ( alpha == 0) 
      { f = -mean(r)+t[p+1]}
    else
      { f =c*((1/sqrt(1+alpha))-((1+alpha)*mean(exp(r))/alpha))}
    return(f)
}


## Obtain the MDPDE in 'est' 
## Default method used in optim is "L-BFGS-B", in case of non-convergence, this method can be changed
  result<-optim(Initial, objfunction_Linear_reg , gr = NULL, 
              method =  "L-BFGS-B", 
              lower = -Inf, upper = Inf, control = list(), hessian = FALSE)
  est=result$par
  s<-exp(est[p+1])  ###sigma estimate
  est[p+1]<-s
  
  
## Record convergence indicator in 'conv'
  conv=result$convergence 
  
  
## Compute asymptotic variance estimate in 'AV'
  varb=s*s*((1+(alpha*alpha/(1+2*alpha)))^(3/2))
  Sb<-varb*diag(T)     #var of beta
  vars=4*(s^4)*((2*(1+2*alpha*alpha)*((1+(alpha*alpha/(1+2*alpha)))^(5/2)))-(alpha*alpha*(1+alpha)*(1+alpha)))/((2+alpha*alpha)^2) #var of sigma
  AV=cbind(t(Sb),vars)

    
## Perform Wald-type test for beta=beta0 (individually) & record P-values in 'PV'
  eb=est[1:p]   #beta estimate
  W<-((eb-beta0)^2)/Sb   #Test statistics
  PV=1-pchisq(t(W),1) #p-values


## Output  
  output<-list(est, AV, PV, conv)
  return(output)
}