
ols <- lm(y.net~.,data = barro)
par(mfrow = c(2,2))
plot(ols,pch = 19,which = 5)
par(mfrow = c(1,1))


# -------------------------------------------------------------------------



#Improved Testimator

library(quantreg)
library(mvtnorm)
library(smoothmest)

beta<-function(p1,p2){
  p = p1 + p2
  set.seed(88888)
  rbinom(p,1,p1/p)
}

LAD.LASSO<-function(X,Y){
  n=nrow(X)
  p=ncol(X)
  lam<-numeric(p)
  e.mat<-diag(p)
  Y = Y - median(Y)
  lad<-coef(rq(Y~ X-1,0.5))
  for(j in 1:p){
    lam[j]= log(n)/(n*abs(lad[j]))
  }
  X.new<-rbind(X, matrix(lam*e.mat,nr=p))
  Y.new<-rbind(Y,matrix(rep(0,p),nr=p))
  lad.lasso<-coef(rq(Y.new ~ X.new-1,0.5))
  return(lad.lasso)
}

uLAD<-function(X,Y,w){
  n=nrow(X)
  p=ncol(X)
  data=as.data.frame(cbind(X,Y))
  UE<-coef(rq(Y ~ X-1,0.5,data,weights = w))
  return(UE)
}

uLAD.1<-function(X,Y,w,index){
  UE.1<- uLAD(X,Y,w)[index]
  return(UE.1)
}

ULAD1 <- function(X,Y,w,index){
  p = ncol(X)
  UE.1<-rep(0,p)
  UE.1[index] <- uLAD(X,Y,w)[index]
  return(UE.1)
}

uLAD.2<-function(X,Y,w){
  UE.2<-uLAD(X,Y,w)[-index]
  return(UE.2)
}

rLAD<-function(X,Y,w,index){
  n=nrow(X)
  p=ncol(X)
  X1<-X[,index]
  data1=as.data.frame(cbind(X1,Y))
  R<-coef(rq(Y ~ X1-1,0.5,data1, weights = w))
  return(R)
}

RLAD <- function(X,Y,w,index){
  p = ncol(X)
  UE.1<-rep(0,p)
  UE.1[index] <- rLAD(X,Y,w,index)
  return(UE.1)
}

Sum.L1<-function(X,Y,b,w){
  n=nrow(X)
  p=ncol(X)
  data=as.data.frame(cbind(X,Y))
  s<-sum(w*abs(X%*%b-Y))
  return(drop(s))
}

L.test<-function(X,Y,index){
  n=nrow(X)
  p=ncol(X)
  X1<-X[,index]
  data=as.data.frame(cbind(X,Y))
  data1=as.data.frame(cbind(X1,Y))
  drop(Sum.L1(X1,Y,rLAD(X,Y,rep(1,n),index),rep(1,n))-Sum.L1(X,Y,uLAD(X,Y,rep(1,n)),rep(1,n)))
}

L.values<-function(X,Y,index){
  n=nrow(X)
  p=ncol(X)
  w<-matrix(rexp(n,rate=1),nr=n)
  bc.s<-rLAD(X,Y,w,index)
  b.s<-uLAD(X,Y,w)
  bc<-rLAD(X,Y,c(rep(1,n)),index)
  b<-uLAD(X,Y,c(rep(1,n)))
  X1<-X[,index]
  Sum.L1(X,Y,b.s,w)-Sum.L1(X1,Y,bc.s,w)-Sum.L1(X,Y,b,w)+Sum.L1(X1,Y,bc,w)
}

alpha=0.05
ptLAD<-function(X,Y,alpha,index){  #Testimator
  L.vec=numeric(length = 1000)
  for(ii in 1:1000){
    L.vec[ii]<-L.values(X,Y,index)
  }
  uLAD.1(X,Y,c(rep(1,n)),index)-as.integer(L.test(X,Y,index)<=quantile(L.vec,1-alpha))*(uLAD.1(X,Y,c(rep(1,n)),index)-rLAD(X,Y,c(rep(1,n)),index))
}

PTLAD <- function(X,Y,alpha,index){
  p = ncol(X)
  UE.1<-rep(0,p)
  UE.1[index] <- ptLAD(X,Y,alpha,index)
  return(UE.1)
}

sLAD<-function(X,Y,index){
  k <- ncol(X) - length(index) - 2
  n <- nrow(X)
  w <- rep(1,n)
  uLAD.1(X,Y,w,index)-(k/L.test(X,Y,index))*(uLAD.1(X,Y,w,index)-rLAD(X,Y,w,index))
}


SLAD <- function(X,Y,index){
  UE.1<-rep(0,p)
  UE.1[index] <- sLAD(X,Y,index)
  return(UE.1)
}

prLAD<-function(X,Y,index){#Improved testimator
  k <- ncol(X) - length(index) - 2
  n <- nrow(X)
  w <- rep(1,n)
  rLAD(X,Y,w,index)+max(1-(k/L.test(X,Y,index)),0)*(uLAD.1(X,Y,w,index)-rLAD(X,Y,w,index))
}

PRLAD <- function(X,Y,index){
  p <- ncol(X)
  UE.1<-rep(0,p)
  UE.1[index] <- prLAD(X,Y,index)
  return(UE.1)
}


p1=3 #Number of nonzero Covariates
p2=5 #Number of Zero Prameter
p=p1+p2
N.sim=1000
MSE.mat<-matrix(0,nr=N.sim,nc=4)

set.seed(88888888)

colnames(MSE.mat) <- c('PTR','PR','LAD-LASSO')
b=as.matrix(beta(p1,p2))



# n = 30
# n = 50
n = 100 #Number of Observations

X<-matrix(rmvnorm(n,rep(0,p),diag(1,p,p)),nr=n,nc=p)

j= 1
while( j < N.sim){
  cat("j = ",j, "\n")
  e<-matrix(rdoublex(n,mu=0,lambda=1),nr=n,nc=1)
  Y<-X%*%b+e
  ladlasso <- LAD.LASSO(X,Y)
  index <- which(abs(ladlasso) > 1e-2)
  ladlasso[-index] <- 0
  index_zero <- which(abs(ladlasso) <= 1e-2)
  if (length(index_zero) < 4)
    next
  print(ladlasso)
  
  MSE.mat[j,1]<-t(PTLAD(X,Y,alpha,index)-b)%*%(PTLAD(X,Y,alpha,index)-b)
  MSE.mat[j,2]<-t(PRLAD(X,Y,index)-b)%*%(PRLAD(X,Y,index)-b)
  MSE.mat[j,3]<-t(ladlasso-b)%*%(ladlasso-b)
  
  j = j + 1
}
res <- apply(MSE.mat,2,mean)
res[1]/res

save(MSE.mat, file = "n100-p1-3-p2-5")
