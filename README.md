# Matrix-variate Differential Network Analysis

## Install
```{r}
# install.packages(“devtools”)
library(devtools)
install_github('jijiadong/MVDN')
```

## Usage

```
mvdn(data1,data2,lambda=NULL,nlambda=20,method=c("none","aic","bic","nbic"),...)
```
### Arguments

  - data1: a three-dimensional *array* (p×q×n1) containing n1 samples.
  - data2:  a three-dimensional *array* (p×q×n2) containing n2 samples.
  - lambda: the tuning parameter of lasso penalty. user-supplied lambda sequence; default is NULL. 
  - nlambda: the number of lambda values, default is 20.
  - method: the method used in the lambda selection.
  - ... :    other arguments that can be passed to *mvdn.cov*.


### Value

  - mvdn: differential network.
  - lambda: the actual sequence of lambda values used.
  - nlambda: the number of lambda values used.
  - opt: indicating which one is the optimized lambda based on different norms ("max" max norm, "1" element-wise, "L1" matrix L1 max norm, "Spectral", "Frobenius", "Nuclear").

### Example

```{r}
##-------------------------------------------------------------------------------
##generate data
p <- 6; q=20; n <- 15
set.seed(1212)
Theta <- graph.generate(p, graph="hub",m=2)
omega1 <- Theta * sample(c(-1, 1), p * p, replace = TRUE) * runif(p * p, 0.1, 0.5)
omega1[lower.tri(omega1, diag = FALSE)] <- 0
omega1 <- omega1 + t(omega1)
diag(omega1) <- abs(min(eigen(omega1)$values)) + 1
sigma1 <- cov2cor(solve(omega1))
omega1 <- solve(sigma1)
omega1[abs(omega1)<10^-4] <- 0
omega2 <- omega1
omega2[1:(p/2),1:(p/2)] <- -1*omega2[1:(p/2),1:(p/2)]
diag(omega2) <- diag(omega1)
sigma2 <- solve(omega2)
delta <- omega2-omega1
sigmaT1 <- 0.4^abs(outer(1:q,1:q,"-"))
sigmaT2 <- 0.5^abs(outer(1:q,1:q,"-"))

X1 <- rmatnorm(n,mean=matrix(0,p,q),sigmaS=sigma1,sigmaT=sigmaT1,method="chol")
X2 <- rmatnorm(n,mean=matrix(0,p,q),sigmaS=sigma2,sigmaT=sigmaT2,method="chol")
##generate data end
##-------------------------------------------------------------------------------

## ------ not run ------
obj <- mvdn(data1=X1,data2=X2,nlambda=10,method="bic")

obj$opt

# use the Frobenius norm
deltahat <- obj$mvdn[[obj$opt["Frobenius"]]]
deltahat

```

