
###################################################
##Generate matrices
#################################################################

genp <- function(p,sparsity,sigma){
  Delta = matrix(0,p,p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      if(runif(1)<sparsity){
        Delta[i,j] = sign(rnorm(1))*sigma
      }
    }
  }
  Delta = Delta+t(Delta)
  meig = eigen(Delta,only.values=T)$values[p]
  Delta = Delta-(meig-1)*diag(p)
  return(Delta)
}

genp1 <- function(p,n,sigma){
  A = matrix(1,p,p)
  l = which(upper.tri(A)==TRUE)
  m = sample(l,n)
  theta = matrix(0,p,p)
  theta[m]=sign(rnorm(n))*sigma
  theta = theta+t(theta)
  return(theta)
}



#################################################################
## d-trace loss algorithm
#################################################################

L1_dts <- function(SigmaX,SigmaY,rho,lambda)
{
  Delta0<-solve(SigmaY+diag(nrow(SigmaY)))-solve(SigmaX+diag(nrow(SigmaX)))
  ##Lambda0<-Delta0
  ##Delta0<-diag(nrow(SigmaX))
  Lambda0<-matrix(1,nrow(SigmaX),ncol(SigmaX))
  tol = 1e-3
  p = dim(Delta0)[1]
  k = 0
  M1 = eigen(SigmaX)
  M2 = eigen(SigmaY)
  Ux = M1$vectors
  Dx = M1$values
  Uy = M2$vectors
  Dy = M2$values
  C1 = matrix(0,p,p)
  C2 = matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      C1[i,j] = 1/(Dy[i]*Dx[j]+2*rho)
      C2[i,j] = 1/(Dy[j]*Dx[i]+2*rho)
    }
  }
  Delta1 = Delta0
  Delta2 = Delta0
  Delta3 = Delta0
  Lambda1 = Lambda0
  Lambda2 = Lambda0
  Lambda3 = Lambda0


  #while (k<10000){
  #z1 = kro(SigmaX-SigmaY+rho*Delta2+rho*Delta3+Lambda3-Lambda1,C1,Ux,Uy)
  #z2 = kro(SigmaX-SigmaY+rho*z1+rho*Delta3+Lambda1-Lambda2,C2,Uy,Ux)
  #z3 = S((Lambda2/rho-Lambda3/rho+z1+z2)/2,(lambda/rho)/2)
  #if (F_dis(Delta1,z1)<tol&&F_dis(Delta2,z2)<tol&&F_dis(z1,z2)<tol)
  #{
  #   break
  #}
  #Delta1 = z1
  #Delta2 = z2
  #Delta3 = z3
  #Lambda1 = Lambda1+rho*(Delta1-Delta2)
  #Lambda2 = Lambda2+rho*(Delta2-Delta3)
  #Lambda3 = Lambda3+rho*(Delta3-Delta1)
  #k = k+1
  #}
  l=p^2
  Delta3 = .C("dtl", as.double(Delta0),Delta3=double(l),as.double(Lambda0),as.double(SigmaX),
              as.double(SigmaY), as.integer(rho), as.double(C1),as.double(C2),
              as.double(Ux),as.double(Uy),as.double(tol),as.integer(p),as.double(lambda))$Delta3;
  Delta3<-matrix(Delta3,nrow=p,ncol=p)
  return (Delta3)
}

#########################################################################


#########################################################################
##lambda selection
#########################################################################
Dpmdtl<- function(X1,X0,
                  lambda=NULL,nlambda=10,lambda.min.ratio=NULL,
                  rho=NULL,shrink=NULL,prec=0.001,
                  correlation=FALSE,
                  tuning=c("none","aic","bic","nbic"))
{
  ## ==========================================================
  ## calculate covariance matrices, calculate lambdas, etc
  ## ==========================================================
  if(ncol(X1)!=ncol(X0))
  {
    cat("X1 and X0 need to have the same number of columns.\n");
    return(NULL);
  }
  n1 <- nrow(X1); n0 <- nrow(X0);
  ## the number of parameters is p(p+1)/2
  p <- ncol(X1);


  if(correlation){ S1 <- cor(X1); S0 <- cor(X0); } else
  { S1 <- cov(X1)*(1-1/n1); S0 <- cov(X0)*(1-1/n0); }
  e <- S1-S0;
  if(!is.null(lambda)){ nlambda <- length(lambda); }
  if(is.null(lambda))
  {
    if(is.null(lambda.min.ratio)){ lambda.min.ratio <- 0.04; }
    lambda.max <- 2*max(abs(e));
    lambda.min <- lambda.min.ratio*lambda.max;
    lambda <- exp(seq(log(lambda.max),log(lambda.min),
                      length=nlambda));
  }


  if(is.null(rho)){ rho <- 1; }
  if(is.null(shrink)){ shrink <- 1.5; }
  lambda <- lambda - shrink*prec;

  ret = vector("list", nlambda);
  for(i in 1:nlambda)
  {
    ret[[i]] <- matrix(NA,nrow=p,ncol=p);
    ret[[i]]<- L1_dts(S1,S0,rho,lambda[i]);

  }


  ## ==============================================================
  ## run tuning
  ## ==============================================================
  opt <- switch(tuning[1], ## default is "none"
                none=NA,
                aic=dpmdtl.ic(S1,S0,ret,n1+n0,2),
                bic=dpmdtl.ic(S1,S0,ret,n1+n0,log(n1+n0)),
                nbic=dpmdtl.ic(S1,S0,ret,n1+n0,2*log(n1+n0)/(n1+n0)))

  if(!is.na(opt[1])){ names(opt) <- c("max","1","L1","sp","F","nc"); }
  return(list(Dpmdtl=ret,lambda=lambda/2,nlambda=nlambda,opt=opt));

}
## **************************************************************
## function to calculate loss
## return a vector of different types of losses
## D=estimated difference matrix
## S1,S0=true, validation set, or training set covariances
## **************************************************************
loss <- function(D,S1,S0)
{
  ##err <- sum(diag(D%*%D%*%(S1%*%S0+S0%*%S1)))/4-sum(diag(D%*%(S1-S0)));
  err <- (S1%*%D%*%S0+S0%*%D%*%S1)/2-S1+S0;
  return(c(max(abs(err)), ## max
           sum(abs(err)), ## l1, element-wise
           max(apply(err,1,function(r){ sum(abs(r)) })), ## matrix L1
           svd(err,nu=0,nv=0)$d[1], ## spectral
           sqrt(sum(err^2)), ## frobenius
           sum(svd(err,nu=0,nv=0)$d))) ## nuclear

}


## **************************************************************
## tuning methods
## **************************************************************
## ==============================================================

dpmdtl.ic <- function(S1,S0,ret,n,penalty)
{
  lowtri <- which(lower.tri(ret[[1]],diag=TRUE));
  df <- sapply(ret,function(x){ sum(x[c(lowtri)]!=0); });
  ic <- scale(n*sapply(ret,loss,S1,S0),
              center=-penalty*df,scale=FALSE);
  return(apply(ic,1,which.min));
}


