\name{mvdn.cov}
\alias{mvdn.cov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MVDN: Matrix-variate Differential Network Analysis.
}
\description{
 Matrix-Variate Differential Network (MVDN) model which is particularly useful in modelling
    connectivity alteration for matrix-variate data such as in brain connectivity analysis.
}
\usage{
mvdn.cov(S1, S2, n1, n2, lambda = NULL, nlambda = 10, method = c("none", "aic", "bic", "nbic"), lambda.min.ratio = NULL, rho = NULL, shrink = NULL, prec = 0.001)
}
\arguments{
  \item{S1}{
A p*p matrix.
}
  \item{S2}{
A p*p matrix.
}
  \item{n1}{
Sample size of group1.
}
  \item{n2}{
Sample size of group2.
}
  \item{lambda}{
The tuning parameter of lasso penalty. user-supplied lambda sequence; default is NULL.
}
  \item{nlambda}{
The number of lambda values, default is 20.
}
  \item{method}{
The method used in the lambda selection.
}
  \item{lambda.min.ratio}{
The ratio used to get the min lambda if \code{lambda} is NULL.
}
  \item{rho}{
The parameter in augmented Lagrange method.
}
  \item{shrink}{
The parameter shrink the \code{lambda}.
}
  \item{prec}{
The parameter shrink the \code{lambda}
}
}

\value{
 \item{mvdn}{differential network.}
 \item{lambda}{the actual sequence of lambda values used.}
 \item{nlambda}{the number of lambda values used.}
 \item{opt}{indicating which one is the optimized lambda based on different norms ("max": max norm, "1": element-wise, "L1": matrix   L1 max norm, "Spectral", "Frobenius", "Nuclear").}
}


\seealso{
    \code{\link{mvdn}}
}
