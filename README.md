# mmlasso
This package provides the functions to calculate the MM-Lasso and adaptive MM-Lasso estimators proposed in Smucler and Yohai (2015).

## Development

If you add new code to `src/fnsauxarma.cpp` you will need to do the following two steps before installing the package.

```
library(Rcpp)
library(RcppArmadillo)
compileAttributes()
```
