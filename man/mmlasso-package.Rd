\name{mmlasso-package}
\alias{mmlasso-package}
\docType{package}

\title{Robust and Sparse Estimators for Linear Regression Models}

\description{
Functions to calculate the MM-Lasso and adaptive MM-Lasso estimators proposed in Smucler and Yohai (2015). The S-Ridge estimator of Maronna (2011) is used as the initial estimator.
}
\details{
\tabular{ll}{
Package: \tab mmlasso\cr
Type: \tab Package\cr
Version: \tab 1.3.4\cr
Date: \tab 2016-2-26\cr
License: \tab GPL (>= 2)\cr
Imports: \tab Rcpp, robustHD, robustbase, parallel, doParallel, foreach, MASS\cr
LinkingTo: \tab Rcpp, RcppArmadillo\cr
NeedsCompilation: \tab yes\cr
}

Index:
\preformatted{
mmlasso        Function to calculate the adaptive MM-Lasso
               and the MM-Lasso
}
\preformatted{
sridge         Function to calculate the S-Ridge
}
}
\author{
Ezequiel Smucler <ezequiels.90@gmail.com>

Maintainer: Ezequiel Smucler <ezequiels.90@gmail.com>
}

\references{
Ezequiel Smucler and Victor J. Yohai. Robust and sparse estimators for linear regression models (2015). Available at \url{http://arxiv.org/abs/1508.01967}.

Maronna, R.A. (2011). Robust Ridge Regression for High-Dimensional Data. Technometrics 53 44-53.
}

\keyword{ package }

\examples{
require(MASS)
p <- 8
n <- 60
rho <- 0.5
desv <- 1
beta.true <- c(rep(0,p+1))
beta.true[2] <- 3
beta.true[3] <- 1.5
beta.true[7] <- 2
mu <- rep(0,p)
sigma <- rho^t(sapply(1:p, function(i, j) abs(i-j), 1:p))
set.seed(1234)
x <- mvrnorm(n,mu,sigma)
u <- rnorm(n)*desv
y <- x\%*\%beta.true[2:(p+1)]+beta.true[1]+u

###Calculate estimators

set.seed(1234)
RobSparse <- mmlasso(x,y)
}
