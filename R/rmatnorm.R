rmatnorm <- function (n, mean = matrix(0, nrow=nrow(sigmaS),ncol=ncol(sigmaT)),
                      sigmaS, sigmaT, method = c("eigen", "svd", "chol")){

  # vec(X) ~ N(vec(M), kron(sigmaT, sigmaS))
  # sigmaS Covariance matrix between rows
  # sigmaT Covariance matrix between columns

  nr <- nrow(sigmaS)
  nc <- ncol(sigmaT)

  if (!isSymmetric(sigmaS, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("sigmaS must be a symmetric matrix")
  }
  if (nrow(mean) != nr)
    stop("nrow of mean and nrow of sigmaS have non-conforming size")

  if (!isSymmetric(sigmaT, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("sigmaT must be a symmetric matrix")
  }
  if (ncol(mean) != nc)
    stop("ncol of mean and ncol of sigmaT have non-conforming size")

  mat = array( dim = c(nr, nc, n) )
  method <- match.arg(method)
  if (method == "eigen") {
    ev <- eigen(sigmaT, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigmaT is numerically not positive semidefinite")
    }
    R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values,0))))
  } else if (method == "svd") {
    s. <- svd(sigmaT)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
      warning("sigmaT is numerically not positive semidefinite")
    }
    R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  } else if (method == "chol") {
    R <- chol(sigmaT, pivot = TRUE)
    R <- R[, order(attr(R, "pivot"))]
  }

  ##
  if (method == "eigen") {
    ev <- eigen(sigmaS, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigmaS is numerically not positive semidefinite")
    }
    Q <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values,0))))
  } else if (method == "svd") {
    s. <- svd(sigmaS)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
      warning("sigmaS is numerically not positive semidefinite")
    }
    Q <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  } else if (method == "chol") {
    Q <- chol(sigmaS, pivot = TRUE)
    Q <- Q[, order(attr(Q, "pivot"))]
  }

  for(i in 1:n){
    mat[,,i] = mean + t(Q) %*% matrix(rnorm(nr*nc), nrow=nr) %*% R
  }

  if(n==1){
    mat <- mat[,,1]
  }
  mat
}


###############

rmatdata <- function (n, mean = matrix(0, nrow=nrow(sigmaS),ncol = ncol(sigmaT)),
                      sigmaS, sigmaT, method = c("eigen", "svd", "chol"),distribution = c("uniform","t"),df){

  # vec(X) ~ distribution(vec(M), kron(sigmaT, sigmaS))
  # sigmaS Covariance matrix between rows
  # sigmaT Covariance matrix between columns

  nr <- nrow(sigmaS)
  nc <- ncol(sigmaT)

  if (!isSymmetric(sigmaS, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("sigmaS must be a symmetric matrix")
  }
  if (nrow(mean) != nr)
    stop("nrow of mean and nrow of sigmaS have non-conforming size")

  if (!isSymmetric(sigmaT, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("sigmaT must be a symmetric matrix")
  }
  if (ncol(mean) != nc)
    stop("ncol of mean and ncol of sigmaT have non-conforming size")

  mat = array( dim = c(nr, nc, n) )
  method <- match.arg(method)
  if (method == "eigen") {
    ev <- eigen(sigmaT, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigmaT is numerically not positive semidefinite")
    }
    R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values,0))))
  } else if (method == "svd") {
    s. <- svd(sigmaT)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
      warning("sigmaT is numerically not positive semidefinite")
    }
    R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  } else if (method == "chol") {
    R <- chol(sigmaT, pivot = TRUE)
    R <- R[, order(attr(R, "pivot"))]
  }

  ##
  if (method == "eigen") {
    ev <- eigen(sigmaS, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigmaS is numerically not positive semidefinite")
    }
    Q <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values,0))))
  } else if (method == "svd") {
    s. <- svd(sigmaS)
    if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
      warning("sigmaS is numerically not positive semidefinite")
    }
    Q <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  } else if (method == "chol") {
    Q <- chol(sigmaS, pivot = TRUE)
    Q <- Q[, order(attr(Q, "pivot"))]
  }

  if(distribution == "uniform"){
    for(i in 1:n){
      mat[,,i] = mean + t(Q) %*% matrix(runif(nr*nc,min=-df,max=df), nrow=nr) %*% R
    }
  } else if(distribution == "t"){
    for(i in 1:n){
      mat[,,i] = mean + t(Q) %*% matrix(rt(nr*nc,df=df), nrow=nr) %*% R
    }
  }

  if(n==1){
    mat <- mat[,,1]
  }
  mat
}
