rm(list = ls())
source("MDPDEinf_LinearReg.r")


###########################################################################
###Data Specification (Example)  
# my.data<-read.delim("Aus_AIDS_data.csv",header=T,sep=",") 
my.data<-salinity
y <- my.data[,4]
n=length(y)
p=4
X<-cbind(matrix(rep(1,n)), matrix(as.matrix(my.data[,c(1,2,3)]),ncol=3))
T<-solve(t(X)%*%X)
Initial<-matrix(c(rep(1,p),0))  #Initial value for estimation interations
beta0<-c(18, 0.77, 0, -0.6)     #the value of \beta to be tested (H_0)


# Compute the estimation and perform test for different alpha values in alpha1
alpha1<-c(seq(0,1,length.out=11))
L<-length(alpha1)
OUT<-matrix(data=seq(0, 0, length.out = 4*L*(p+1)),nrow=p+1,ncol=4*L)
H <-cbind(t(alpha1), t(alpha1), t(alpha1), t(alpha1))

for (j in 1:L){
  
  alpha<-alpha1[j]
  result<-lmdpd(y,X,alpha,p,Initial,beta0)
  
  #Estimates
  OUT[,j]<-result[[1]]
  #Asymptotic variances
  OUT[,L+j]<-result[[2]]
  # p-values
  OUT[,2*L+j]<-rbind(t(result[[3]]), matrix(-1))
  #Record convergence indicator
  OUT[,3*L+j]=result[[4]]
}

#Output has 4 block of columns each of size L (columns corr to alpha values in each block)
#First block contains the MDPDEs of beta and sigma (nrow = p +1)
#Second block contains their asymptotic variance estimates
#3rd block contains p-value of the test for beta = beta0 componentwise
#4th block contains the convergence indicator 

#Write output in a csv file
OUT<-rbind(H, OUT)
write.csv(OUT, file = "Result.csv")
