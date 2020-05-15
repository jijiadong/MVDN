\name{rmatnorm}
\alias{rmatnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix normal distribution
}
\description{
The function provide a random number generator for the Matrix normal distribution \code{vec(X) ~ N(vec(M), kron(sigmaT, sigmaS))}.
}
\usage{
rmatnorm(n, mean = matrix(0, nrow = nrow(sigmaS), ncol = ncol(sigmaT)), sigmaS, sigmaT, method = c("eigen", "svd", "chol"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of observations.
}
  \item{mean}{
mean matrix, default is \code{matrix(0, nrow=nrow(sigmaS), ncol=ncol(sigmaT))}.
}
  \item{sigmaS}{
covariance matrix between rows.
}
  \item{sigmaT}{
covariance matrix between columns.
}
  \item{method}{
   string specifying the matrix decomposition used todetermine the matrix root of \code{sigmaS} and \code{sigmaT}. Possible methods are eigenvalue decomposition (\code{"eigen"}, default), singular value decomposition (\code{"svd"}), and Cholesky decomposition (\code{"chol"}).
}
}

\value{
returns a p*q*n array if n > 1. returns a p*q matrix if n = 1.
}


\seealso{
 \code{\link{rmatdata}}
}
\examples{

p <- 5; q <- 6; n <- 10
set.seed(1212)
s1 <- 0.5^abs(outer(1:p,1:p,"-"))
s2 <- 0.4^abs(outer(1:q,1:q,"-"))
x <- rmatnorm(n,mean=matrix(1,p,q),sigmaS=s1,sigmaT=s2)


}

