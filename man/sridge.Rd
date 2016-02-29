\name{sridge}
\alias{sridge}
\title{Robust and Sparse Estimators for Linear Regression Models}
\description{
Function to calculate the S-Ridge estimator of Maronna (2011), adapted from Maronna's original MATLAB code. It is not intended to be used on its own, but rather as an initial estimator for the MM-Lasso.
}


\usage{sridge(x,y,cualcv.S=5,numlam.S=30,niter.S=50,normin=0,
          denormout=0,alone=0,ncores=1)}
\arguments{
  \item{x}{A matrix of carriers. Intercept is added automatically.}
  \item{y}{A vector of response variables.}
    \item{cualcv.S}{A natural number greater than 2. Method for estimating prediction error of S-Ridge: cualcv-fold cross validation. Default is 5.}
  \item{numlam.S}{Number of candidate penalization parameter values for S-Ridge. Default is 30.}
  \item{niter.S}{Maximum number of iterations of IWLS for S-Ridge. Default is 50.}
  \item{normin}{Center and scale input data? 0=no, default ; 1=yes.}
  \item{denormout}{Return final estimate in the original coordinates? 0=no, default ; 1=yes.}
  \item{alone}{Are you calculating the estimator for its sake only? 0=no, default ; 1=yes.}
  \item{ncores}{Number of cores to use for parallel computations. Default is one core.}

}

\value{
 \item{coef}{S-Ridge estimate. First coordinate is the intercept.}
 \item{scale}{M-estimate of scale of the residuals of the final regression estimate.}
 \item{edf}{Final equivalent degrees of freedom.}
\item{lamda}{Optimal lambda.}
\item{delta}{Optimal delta.}
}
\author{
Ezequiel Smucler, \email{ezequiels.90@gmail.com}.

Ricardo Maronna.
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
SRidge <- sridge(x,y,normin=1,denormout=1,alone=1)
}

\keyword{models}
\keyword{regression}
\keyword{robust}


