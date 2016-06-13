k=30
x<-c(0,0,0,1,2,1,0,2,0,0,0,1,3,0,0,3,0,2,1,2,0,0,0,1,5,0,1,0,0,1,0,1,0,1,0,1)
true.lambda<-c(1.17,.83,.5,1,.83,.83,1.17,.83,.67,.17,0,.33,1.5,.5,1.17,1.33,.5,1.17,
               .5,.5,1.33,.83,.33,1.5,1.33,.67,.67,.33,.33,.33,.5,.83,.67,.33,0,.5)
em.lambda<-29/65*x
# function to make A-matrix
make_big_A<-function(x,lambda)
{
  t(apply(as.matrix(x),1,function(s) {exp(-lambda)*lambda^s/factorial(s) } ))
  
}


# EM

EM <- function(y_init,A,eps=1e-6,max.iter=1e6){
  y_old <- y_init
  n <- nrow(A)
  J <- 1
  err <- Inf
  while(err > eps & J <= max.iter){
    y_new <- -grad(y_old,A)*y_old/n
    err <- abs(obj(y_new,A) - obj(y_old,A))/abs(obj(y_old,A))
    y_old <- y_new
    J <- J+1  
  }
  list(y=y_new,iterations=J)
}


#
# Auxiliary functions (for estimating mixing distribution)
#

# compute neg. log-likelhood

obj <- function(y,A) -sum(log(A%*%y))

# compute gradient

grad <- function(y,A) -t(A)%*%(1/(A%*%y))
#######################
#######################


lambda_MLE<-x
max.lambda_MLE<-max(lambda_MLE)
lambda_NP<-seq(1/6,max.lambda_MLE,length=k)
A_NP<-make_big_A(x,lambda_NP)
fit_NP <- EM(y_init=rep(1/k,k),A=A_NP) 
y_NP<-fit_NP$y
plot(lambda_NP,y_NP,type="l")
est.lambda_NP <- t(A_NP%*%(y_NP*lambda_NP)/(A_NP%*%y_NP))
est.try<-est.lambda_NP-est.lambda_NP[1]

est.lambda_NP[x==0]<-0
est.lambda_NP
est.exc0<-est.lambda_NP[true.lambda!=0]
est.exc0
true.lambda.exc0<-true.lambda[true.lambda!=0]
sse.lambdaLoss<-sum((est.exc0-true.lambda.exc0)^2/true.lambda.exc0)
sse.lambdaLoss

est.try
est.try.exc0<-est.try[true.lambda!=0]
sse<-sum((est.try.exc0-true.lambda.exc0)^2/true.lambda.exc0)
sse
  

# exclude x=0 at first
A_NP.exc0<-make_big_A(x[x!=0],lambda_NP)
fit_NP.exc0 <- EM(y_init=rep(1/k,k),A=A_NP.exc0) 
y_NP.exc0<-fit_NP.exc0$y
plot(lambda_NP,y_NP.exc0,type="l")
est.lambda_NP.exc0 <- t(A_NP.exc0%*%(y_NP.exc0*lambda_NP)/(A_NP.exc0%*%y_NP.exc0))
est.lambda_NP.all<-rep(0,length(x))
est.lambda_NP.all[x!=0]<-est.lambda_NP.exc0
est.lambda_NP.all.exc0<-est.lambda_NP.all[true.lambda!=0]
sse.lambdaLoss.exc0<-sum((est.lambda_NP.all.exc0-true.lambda.exc0)^2/true.lambda.exc0)
sse.lambdaLoss.exc0

