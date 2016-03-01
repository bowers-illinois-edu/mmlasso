sridge <- function(x, y, cualcv.S=5, numlam.S=30, niter.S=50, normin=0, denormout=0, alone=0, ncores=1){
    # Solves n*s_n^2 +lam*||beta1||^2} = min. Adapted from Ricardo Maronna's original MATLAB code.
    # INPUT
    # cualcv.S: method for estimating prediction error. cualcv-fold cross-validation
    # normin: normalize input data?. 0=no,  default ; 1=yes
    # denormout: denormalize output?. 0=no,  default ; 1=yes
    # alone: are you calculating the estimator for its sake only? 0=no,  default ; 1=yes
    # numlam.S: number of lambda values,  default 30
    # niter.S : number of maximum iterations of IWLS
    # ncores : number of cores to use for parallel computations. Default is one core.
    # OUTPUT
    # coef: (p+1)-vector of regression parameters,  beta(1)=intercept
    # scale: M-estimate of scale of the final regression estimate
    # edf: final equivalent degrees of freedom
    # lamda: optimal lambda
    # delta: optimal delta


    ###Center and scale covariates and response using median and mad
    if (normin == 1){
        prep <- prepara(x,  y)
        Xnor <- prep$Xnor
        ynor <- prep$ynor
        mux <- prep$mux
        sigx <- prep$sigx
        muy <- prep$muy
    }else{
        Xnor <- x
        ynor <- y
    }

    # Spherical Principal Components (no centering) privar,  Beig= vector of robust
    # "eigenvalues" and matrix of eigenvectors Xnor is now = PCA scores
    pca <- SPCC(Xnor)
    privar <- pca$lamda
    Beig <- pca$b
    Xnor <- pca$scores
    n <- nrow(Xnor)
    p <- ncol(Xnor)
    nkeep <- 5 # number of candidates to be kept for full iteration in the Pena-Yohai procedure
    privar <- privar*n # Makes the robust eigenvalues of the same order as those of classical PCA used for LS
    nlam <- min(c(p, numlam.S))
    pmax <- min(c(p,  n/2))   # edf<=n/2 to keep BDP >=0.25
    pp <- seq(1, pmax, length=nlam)
    lamdas <- findlam(privar, pp) # find lambdas corresponding to the edf's
    deltas <- 0.5*(1-(pp)/n)  # for the M-scale used with Penia-Yohai

    ### Reorder data for CV
    srt <- sort.int(sample(1:n), index.return=TRUE)  # Random permutation
    tt <- srt$x
    orden <- srt$ix
    Xnord <- Xnor[orden, ]
    ynord <- ynor[orden]
    ###

    if (alone==1){
        if(length(find("cl")) != 1){
            newclus <- TRUE
            ### Set-up cluster for parallel computations
            cores <- min(detectCores(), ncores)
            try(cl <- makeCluster(cores))
        }
        newclus <- FALSE
        try(registerDoParallel(cl))
        locallib <- .libPaths()[[1]]
        clusterExport(cl,  "locallib",  envir=environment())
        clusterEvalQ(cl, .libPaths(locallib))
        ###
    }
    if(alone == 0){
        if(length(find("cl") != 1)){
            newclus <- TRUE
        }
        newclus <- FALSE
    }

    ### Parallel CV
    exp.se <- c("CVSE")
    klam <- NULL
    mse <- foreach(klam = 1:length(lamdas),  .combine = c,  .packages = c("mmlasso"),
                   .export = exp.se) %dopar% {
        CVSE(Xnord, ynord, nfold=cualcv.S, lamdas[klam], pp[klam], nkeep, niter.S)
    }

    if(any(is.infinite(mse))){
        warning("IWLS for S-Ridge failed when lambda equaled:")
        print(lamdas[which(is.infinite(mse))])
    }
    indmin <- which.min(mse)
    lamin <- lamdas[indmin]
    delmin <- deltas[indmin]
    ###

    if (alone == 1 & newclus){
        stopCluster(cl)
    }

    fin <- rr_se(Xnor, ynor, lamin, deltaesc=delmin, cc_scale=1, nkeep, niter.S, epsilon=1e-4)
    beta <- fin$coef
    betaslo <- beta[2:(length(beta))]
    bint <- beta[1]
    res <- ynor-Xnor%*%betaslo-as.vector(bint)
    edf <- fin$edf
    deltult <- 0.5 * (1 - (edf + 1)/n)  ##  'delta' for final scale
    deltult <- max(c(deltult,  0.25))
    # c0: constant for consistency of final scale
    c0 <- const_marina(deltult)
    sigma <- Mscale_mar(res, deltult, c0)
    a_cor <- mean(psi_marina(res/sigma, c0)^2)
    b_cor <- mean(Mchi(res/sigma,  c0,  "bisquare",  2))
    c_cor <- mean(psi_marina(res/sigma, c0)*(res/sigma))
    corr <- 1+(edf+1)/(2*n)*(a_cor/(b_cor*c_cor))
    fscale <- sigma*corr  # bias correction for final scale
    # Back from PC to ordinary coordinates
    betaslo <- Beig%*%betaslo

    # De-normalize beta if required by user
    if (normin==1 & denormout==1){
        betaslo <- betaslo/sigx
        bint <- muy+bint-mux%*%betaslo
    }
    beta <- c(bint, betaslo)
    re <- list(coef=beta, scale=fscale, edf=edf, lamda=lamin, delta=delmin)
    return(re)
}


CVSE <- function(X, y, nfold, lam, gradlib, nkeep, niter.S){
    # Performs nfold-CV for S-Ridge
    # INPUT
    # beta.ini,  scale.ini: initial estimates of regression and scale
    # X, y: data
    # lam: penalization parameter
    # gradlib: degrees of freedom,  used to calculate the delta for the M-scale
    # nkeep: number of candidates for full iterations of IWLS
    # niter.S : number of maximum iterations of IWLS
    # OUTPUT
    # mse: resulting MSE (estimated using a tau-scale)

    ### Segment data
    n <- nrow(X)
    p <- ncol(X)
    indin <- 1:n
    nestim <- n*(1-1/nfold)
    lamcv <- lam
    deltaesc <- 0.5*(1-(gradlib)/nestim)
    inint <- floor(seq(0, n, length.out=nfold+1))
    resid <- vector(mode = "numeric",  length = n)
    ###

    for (kk in 1:nfold){
        testk <- (inint[kk]+1):inint[kk+1]
        estik <- setdiff(indin,  testk)
        Xtest <- X[testk, ]
        Xesti <- X[estik, ]
        ytest <- y[testk]
        yesti <- y[estik]
        se <- try(rr_se(Xesti, yesti, lamcv, deltaesc, cc_scale=1, nkeep, niter.S, epsilon=1e-4))
        if (class(se)=="try-error"){
            return(Inf)
        }
        bet <- se$coef
        bint <- bet[1]
        beta <- bet[2:(p+1)]
        fitk <- Xtest%*%beta+bint
        resid[testk] <- ytest-fitk
    }
    mse <- scale_tau(resid)
    return(mse)
}

findlam <- function(vals, r){
    # Finds lamdas which yield edf=r
    p <- length(vals)
    nr <- length(r)
    lamr <- vector(mode = "numeric",  nr)
    lam1 <- 0
    lam2 <- max(vals)*(p/r-1)
    for (i in 1:nr){
        lam0 <- c(lam1,  lam2[i]+0.5)   # the value 0.5 is for lam=0
        lamr[i] <- uniroot(sumrat, lam0, vals, r[i])$root
    }
    return(lamr)
}

sumrat <- function(lam, vals, r){
    susu <- sum(vals/(vals+lam))-r
    return(susu)
}


const_marina <- function(delta){
    integrand <-  function(x, c){dnorm(x)*rho_marina(x, c)}
    expe <- function(c, delta){integrate(integrand, lower=-Inf, upper=Inf, c)$value-delta}
    init <- c(0.1, 100)
    try(cw <- uniroot(expe, init, delta)$root, silent=TRUE)
    if (class(cw)=="try-error"){
        warning("Something's wrong,  could not find tuning constant for the scale")
        return(NULL)
    }
    return(cw)
}


#' A function to fit a ridge penalized MM-estimator
#' This version does not do any cross-validation. It merely returns the fits at the chosen lambdas.
#' @param lambdas are a vector of lambda values or a scalar lambda value.
#' @param deltas should be the delta values associated with the lambda values.
#' @export
sridge2 <- function(x,  y,  cualcv.S = 5,  nkeep = 5,  numlam.S = 30,  niter.S = 50,  normin = 0,
                    denormout = 0,  alone = 0,  ncores = 1,  lamdas = NULL,  deltas = NULL,  adaptive = FALSE) {
    #  Solves n*s_n^2 +lam*||beta1||^2 = min. Adapted from Ricardo Maronna's original
    #  MATLAB code. INPUT cualcv.S: method for estimating prediction error.
    #  cualcv-fold cross-validation nkeep: number of candidates to be kept for full
    #  iteration in the Pena-Yohai procedure (default=5) normin: normalize input
    #  data?. 0=no,  default ; 1=yes denormout: denormalize output?. 0=no,  default ;
    #  1=yes alone: are you calculating the estimator for its sake only? 0=no,  default
    #  ; 1=yes numlam.S: number of lambda values,  default 30 niter.S : number of
    #  maximum iterations of IWLS ncores : number of cores to use for parallel
    #  computations. Default is all available cores OUTPUT beta: (p+1)-vector of
    #  regression parameters,  beta(1)=intercept fscale: M-estimate of scale of the
    #  final regression estimate edf: final equivalent degrees of freedom lamin:
    #  optimal lambda delmin: optimal delta

    ###  Center and scale covariates and response using median and mad
    if (normin == 1) {
        stop("Not able to scale by mad when most of variables are binary. Maybe will scale in some other way later.")
        prep <- prepara(x,  y)
        Xnor <- prep$Xnor
        ynor <- prep$ynor
        mux <- prep$mux
        sigx <- prep$sigx
        muy <- prep$muy
    } else {
        Xnor <- x
        ynor <- y
        sigx <- rep(1, ncol(x))
    }

    #  Spherical Principal Components (no centering) privar,  Beig= vector of robust
    #  "eigenvalues" and matrix of eigenvectors Xnor is now = PCA scores
    pca <- SPCC(Xnor)
    privar <- pca$lamda
    Beig <- pca$b
    Xnor <- pca$scores
    n <- nrow(Xnor)
    p <- ncol(Xnor)
    if(is.null(lamdas)){ ##  Create lambda vector if one was not provided
        privar <- privar * n  ##  Makes the robust eigenvalues of the same order as those of classical PCA used for LS
        nlam <- min(c(p,  numlam.S))
        pmax <- min(c(p,  n/2))  # edf<=n/2 to keep BDP >=0.25
        pp <- seq(1,  pmax,  length = nlam)
        lamdas <- findlam(privar,  pp)  # find lambdas corresponding to the edf's
        deltas <- 0.5 * (1 - (pp)/n)  # for the M-scale used with Penia-Yohai
    }

    fin <- rr_se_vec(X = Xnor,  y = ynor,  lambda2 = lamdas,
                     deltaesc = deltas,  cc_scale = 1,  nkeep,  niter.S,  epsilon = 1e-04)
    betas <- fin$coef
    res <- apply(betas,2,function(b){
                     ynor - Xnor %*% b[-1] - as.vector(b[1])
                     })
    edf <- fin$edf
    deltult <- 0.5 * (1 - (edf + 1)/n)  ##  'delta' for final scale
    ## deltult <- max(c(deltult,  0.25))
    ## deltult <- pmax(deltult,  rep(0.25,length(deltult)))
    deltult <- deltult * (deltult > .25) + .25 * ( 1 - (deltult > .25) )
    #  c0: constant for consistency of final scale
    getfscale<-function(thedeltult,theres,theedf){
        c0 <- const_marina(thedeltult)
        sigma <- Mscale_mar(theres,  thedeltult,  c0)
        a_cor <- mean(psi_marina(theres/sigma,  c0)^2)
        b_cor <- mean(Mchi(theres/sigma,  c0,  "bisquare",  2))
        c_cor <- mean(psi_marina(theres/sigma,  c0) * (theres/sigma))
        corr <- 1 + (theedf + 1)/(2 * n) * (a_cor/(b_cor * c_cor))
        fscale <- sigma * corr  # bias correction for final scale
        return(fscale)
    }
    fscale <- sapply(1:ncol(res),function(i){
                         getfscale(thedeltult=deltult[i],
                                   theres=res[,i],theedf=edf[i]) })
    #  Back from PC to ordinary coordinates
    betaslo <- apply(betas[-1,],2,function(b){ Beig %*% b})

    if(!adaptive){
        lamdaad <- t(as.matrix(c(lamda=NA,delta=NA)))
        betares <- rbind(betas[1,],  betaslo)
        results <- list(coef = betares,  scale = fscale,  edf = edf,  lamda = lamdas,  delta = deltas,
                        lamdaad = lamdaad[,"lamda"], deltaad = lamdaad[,"delta"])
        return(results)
    } else {
        resmaker <- function( Xnor,  ynor,  nkeep,  niter.S,  Beigi, numlam.S){
            force(Xnor); force(ynor); force(nkeep); force(niter.S); force(Beig); force(numlam.S) ##  ; force(mux)
            function(thelamda, thedelta, thebetas){
                activ <- which(abs(thebetas[2:length(thebetas)]) > .Machine$double.eps)
                ###Take out covariates not chosen by MMLasso and scale the remaining ones
                wad <- abs(thebetas[2:length(thebetas)][activ])
                Xnorw <- as.matrix(scale(Xnor[, activ], center=FALSE,  scale = 1/wad))
                ###Calculate candidate lambdas for adaptive MM-ridge
                lamdaad <- genlambdas(Xnorw, ynor, numlam.S=numlam.S)
                lamdaad <- lamdaad[lamdaad[,"lamda"]!=0,]
                fin <- rr_se_vec(X = Xnorw,  y = ynor,  lambda2 = lamdaad[,"lamda",drop=FALSE],
                                 deltaesc = lamdaad[,"delta",drop=FALSE],  cc_scale = 1,  nkeep,  niter.S,  epsilon = 1e-04)
                betas <- fin$coef
                betasbig <- matrix(0,ncol=p,nrow=nrow(fin$coef))
                betasbig[,activ] <- fin$coef
                betasbig[-1,] <- betasbig[-1,]*wad
                res <- apply(betasbig,2,function(b){
                                 ynor - Xnorw %*% b[-1] - as.vector(b[1])
                                 })
                edfbig<-rep(0,p)
                edfbig[activ] <- fin$edf
                deltult <- 0.5 * (1 - (edfbig + 1)/n)  ##  'delta' for final scale
                deltult <- deltult * (deltult > .25) + .25 * ( 1 - (deltult > .25) )
                fscale <- sapply(1:ncol(res),function(i){
                                     getfscale(thedeltult=deltult[i],
                                               theres=res[,i],theedf=edf[i]) })
                #  Back from PC to ordinary coordinates
                betasload <- apply(betasbig[-1,],2,function(b){ Beig %*% b})
                beta <- rbind(betas[1,],  betasload)
                results <- list(coef = betasload,  scale = fscale,  edf = edf,  lamda = thelamda,  delta = thedelta,
                                lamdaad = lamdaad[,"lamda"], deltaad = lamdaad[,"delta"])
                return(results)

                ##  De-normalize beta if required by user
                ##  if (normin == 1 & denormout == 1) {
                ##      betaslo <- betaslo/sigx
                ##      bint <- muy + bint - mux %*% betaslo
                ##  }
            }
        }

        getresults <- resmaker(Xnor,  ynor,  nkeep,  niter.S,  Beig, numlam.S)

        if(ncores > 1){
            ##  assuming name of cluster is cl for now
            clusterExport(cl, c("getresults","getfscale"), envir=environment())
            results <- parLapply(cl, 1:ncol(betas),  function(i){
                                     getresults(thelamda=lamdas[[i]], thedelta=deltas[[i]], thebetas=betas[,i])
                                })
        } else {
            results <- lapply(1:ncol(betas), function(i){
                                  message(i)
                                  getresults(thelamda=lamdas[[i]], thedelta=deltas[[i]], thebetas=betas[,i])
                                })
        }

    }
}

#' Generate lambda sequences
genlambdas <- function(x,  y,  normin = 0,  numlam.S = 30) {

    ###  Center and scale covariates and response using median and mad
    if (normin == 1) {
        stop("Not able to scale by mad when most of variables are binary. Maybe will scale in some other way later.")
        prep <- prepara(x,  y)
        Xnor <- prep$Xnor
        ynor <- prep$ynor
        mux <- prep$mux
        sigx <- prep$sigx
        muy <- prep$muy
    } else {
        Xnor <- x
        ynor <- y
    }

    #  Spherical Principal Components (no centering) privar,  Beig= vector of robust
    #  'eigenvalues' and matrix of eigenvectors Xnor is now = PCA scores
    pca <- SPCC(Xnor)
    privar <- pca$lamda
    Beig <- pca$b
    Xnor <- pca$scores
    n <- nrow(Xnor)
    p <- ncol(Xnor)
    privar <- privar * n  # Makes the robust eigenvalues of the same order as those of classical PCA used for LS
    nlam <- min(c(p,  numlam.S))
    pmax <- min(c(p,  n/2))  # edf<=n/2 to keep BDP >=0.25
    pp <- seq(1,  pmax,  length = nlam)
    lamdas <- findlam(privar,  pp)  # find lambdas corresponding to the edf's
    deltas <- 0.5 * (1 - (pp)/n)  # for the M-scale used with Penia-Yohai

    return(cbind(lamda=lamdas, delta=deltas))
}
