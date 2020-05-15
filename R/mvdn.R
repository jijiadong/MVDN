mvdn <- function(data1,data2,lambda=NULL,nlambda=20,method=c("none","aic","bic","nbic"),...){

  if(class(data1)!="array" | class(data2)!="array"){
    stop("The class of data1 and data2 must be 'array' ")
  }

  if(dim(data1)[1]!=dim(data2)[1]){
    stop("data1 and data2 must have the the same number of rows")
  }

  if(dim(data1)[2]!=dim(data2)[2]){
    stop("data1 and data2 must have the the same number of columns")
  }

  p <- dim(data1)[1]
  q <- dim(data1)[2]
  n1 <- dim(data1)[3]
  n2 <- dim(data2)[3]

  data1 <- arr2list(data1)
  data2 <- arr2list(data2)

  X1 <- scale(t(sapply(data1,unlist)),center=T,scale=F)
  X2 <- scale(t(sapply(data2,unlist)),center=T,scale=F)

  X1 <- toList(X1,dims=c(p,q))
  X2 <- toList(X2,dims=c(p,q))

  covs <- function(xx){
    xx %*% t(xx)
  }

  gamma1 <- cov2cor(array(rowMeans(sapply(X1,covs),na.rm=T),dim=c(p,p)))
  gamma2 <- cov2cor(array(rowMeans(sapply(X2,covs),na.rm=T),dim=c(p,p)))

  optmvdn <- mvdn.cov(S1=gamma1,S2=gamma2,n1=n1,n2=n2,lambda=lambda,nlambda=nlambda,method=method,...)
  optmvdn
}



mvdn.cov <- function(S1,S2,n1,n2,lambda=NULL,nlambda=10,method=c("none","aic","bic","nbic"),
                     lambda.min.ratio=NULL,rho=NULL,shrink=NULL,prec=0.001){

  # if (!isSymmetric(S1) | !isSymmetric(S2) ){
  # stop("S1 and S2 must be symmetric matrix")
  # }

  p <- ncol(S1)

  if(!is.null(lambda)){ nlambda <- length(lambda) }

  if(is.null(lambda)){

    if(is.null(lambda.min.ratio)){ lambda.min.ratio <- 0.04 }

    lambda.max <- 2*max(abs(S2-S1))
    lambda.min <- lambda.min.ratio*lambda.max
    lambda <- exp(seq(log(lambda.max),log(lambda.min),length=nlambda))
  }


  if(is.null(rho)){ rho <- 1 }
  if(is.null(shrink)){ shrink <- 1.5 }

  lambda <- lambda - shrink*prec

  ret <- vector("list", nlambda)
  for(i in 1:nlambda){
    ret[[i]] <- matrix(NA,nrow=p,ncol=p)
    ret[[i]]<- L1_dts(S1,S2,rho,lambda[i])
  }

  ## run tuning
  opt <- switch(method[1], ## default is "none"
                none=NA,
                aic=dpmdtl.ic(S1,S2,ret,n1+n2,2),
                bic=dpmdtl.ic(S1,S2,ret,n1+n2,log(n1+n2)),
                nbic=dpmdtl.ic(S1,S2,ret,n1+n2,2*log(n1+n2)/(n1+n2)))

  if(!is.na(opt[1])){ names(opt) <- c("max","1","L1","Spectral","Frobenius","Nuclear") }

  return(list(mvdn=ret,lambda=lambda/2,nlambda=nlambda,opt=opt))
}



arr2list <- function(x){
  # x is a 3 dims array
  xlist <- lapply(1:dim(x)[3],function(k){x[,,k]})
  names(xlist) <- dimnames(x)[[3]]
  xlist
}



toList <- function(xx,dims){
  MM <- list()
  for(i in 1:nrow(xx)){
    MM[[i]] <- array(xx[i,],dim=dims)
  }
  MM
}

