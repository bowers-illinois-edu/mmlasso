\name{mmlasso}
\alias{mmlasso}
\title{Robust and Sparse Estimators for Linear Regression Models}
\description{
Function to calculate the MM-Lasso and adaptive MM-Lasso estimators proposed in Smucler and Yohai (2015). The S-Ridge estimator of Maronna (2011) is used as the initial estimator.
}


\usage{mmlasso(x, y, varsigma=1, cualcv.mm=5, cualcv.S=5,
        numlam.mm=30, numlam.S=30, niter.S=50, niter.mm=50, ncores=1)}
\arguments{
  \item{x}{A matrix of carriers. Intercept is added automatically.}
  \item{y}{A vector of response variables.}
  \item{varsigma}{Power to elevate the weights for the adaptive MM-Lasso. Default is 1.}
  \item{cualcv.mm}{A natural number greater than 2. Method for estimating prediction error of MM and adaptive MM-Lasso: cualcv-fold cross validation. Default is 5.}
  \item{cualcv.S}{A natural number greater than 2. Method for estimating prediction error of S-Ridge: cualcv-fold cross validation. Default is 5.}
  \item{numlam.mm}{Number of candidate penalization parameter values for MM and adaptive MM-Lasso. Default is 30.}
  \item{numlam.S}{Number of candidate penalization parameter values for S-Ridge. Default is 30.}
  \item{niter.mm}{Maximum number of weighted Lasso iterations for MM and adaptive MM-Lasso. Default is 50.}
  \item{niter.S}{Maximum number of iterations of IWLS for S-Ridge. Default is 50.}
  \item{ncores}{Number of cores to use for parallel computations. Default is one core.}
}

\value{
 \item{coef.SE}{Initial S-Ridge estimate. First coordinate is the intercept.}
 \item{coef.MMLasso}{MM-Lasso estimate. First coordinate is the intercept.}
 \item{coef.MMLasso.ad}{adaptive MM-Lasso estimate. First coordinate is the intercept.}

}
\author{
Ezequiel Smucler, \email{ezequiels.90@gmail.com}.
}
\references{
Ezequiel Smucler and Victor J. Yohai. Robust and sparse estimators for linear regression models (2015). Available at \url{http://arxiv.org/abs/1508.01967}.

Maronna, R.A. (2011). Robust Ridge Regression for High-Dimensional Data. Technometrics 53 44-53.
}

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

\keyword{models}
\keyword{regression}
\keyword{robust}


