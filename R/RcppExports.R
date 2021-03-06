# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rho_marina <- function(x, cw) {
    .Call('mmlasso_rho_marina', PACKAGE = 'mmlasso', x, cw)
}

rho_bisq <- function(x, cw) {
    .Call('mmlasso_rho_bisq', PACKAGE = 'mmlasso', x, cw)
}

psi_marina <- function(x, cw) {
    .Call('mmlasso_psi_marina', PACKAGE = 'mmlasso', x, cw)
}

weight_marina <- function(x, cw) {
    .Call('mmlasso_weight_marina', PACKAGE = 'mmlasso', x, cw)
}

Mscale_mar <- function(x, b, cc) {
    .Call('mmlasso_Mscale_mar', PACKAGE = 'mmlasso', x, b, cc)
}

Mscale_bisq <- function(x, b, cc, cuad) {
    .Call('mmlasso_Mscale_bisq', PACKAGE = 'mmlasso', x, b, cc, cuad)
}

scale_tau <- function(x) {
    .Call('mmlasso_scale_tau', PACKAGE = 'mmlasso', x)
}

spa_med <- function(x) {
    .Call('mmlasso_spa_med', PACKAGE = 'mmlasso', x)
}

rob_sq <- function(z) {
    .Call('mmlasso_rob_sq', PACKAGE = 'mmlasso', z)
}

my_svdecon <- function(x) {
    .Call('mmlasso_my_svdecon', PACKAGE = 'mmlasso', x)
}

SPCC <- function(x) {
    .Call('mmlasso_SPCC', PACKAGE = 'mmlasso', x)
}

MMLassoCpp_ini <- function(xx, y, beta_ini) {
    .Call('mmlasso_MMLassoCpp_ini', PACKAGE = 'mmlasso', xx, y, beta_ini)
}

MMLassoCpp1 <- function(x, y, beta_ini, scale_ini, c1) {
    .Call('mmlasso_MMLassoCpp1', PACKAGE = 'mmlasso', x, y, beta_ini, scale_ini, c1)
}

MMLassoCpp2 <- function(xjota, yast, beta_lars, beta_o, alpha) {
    .Call('mmlasso_MMLassoCpp2', PACKAGE = 'mmlasso', xjota, yast, beta_lars, beta_o, alpha)
}

desrobrid <- function(x, y, niter, lam, betinte, betinslo, cc, delsca, epsilon) {
    .Call('mmlasso_desrobrid', PACKAGE = 'mmlasso', x, y, niter, lam, betinte, betinslo, cc, delsca, epsilon)
}

regrid <- function(x, y, lambda) {
    .Call('mmlasso_regrid', PACKAGE = 'mmlasso', x, y, lambda)
}

rr_se <- function(X, y, lambda2, deltaesc, cc_scale, nkeep, niter, epsilon) {
    .Call('mmlasso_rr_se', PACKAGE = 'mmlasso', X, y, lambda2, deltaesc, cc_scale, nkeep, niter, epsilon)
}

rr_se_vec <- function(X, y, lambda2, deltaesc, cc_scale, nkeep, niter, epsilon) {
    .Call('mmlasso_rr_se_vec', PACKAGE = 'mmlasso', X, y, lambda2, deltaesc, cc_scale, nkeep, niter, epsilon)
}

