mmlasso<-function(x,y,varsigma=1,cualcv.mm=5,cualcv.S=5,numlam.mm=30,numlam.S=30,niter.S=50,niter.mm=50,ncores=Inf){
  #Main function to compute MM-Lasso estimates.
  #The initial estimate is the S-Ridge of Maronna (2012).
  #INPUT
  #X, y data set
  #cualcv.mm: method for estimating prediction error of MM and adaptive MM-Lasso, cualcv-fold CV
  #cualcv.S: method for estimating prediction error of S-Ridge, cualcv-fold CV
  #numlam.mm: number of candidate lambda values for MM and adaptive MM-Lasso
  #numlam.S: number of candidate lambda values for S-Ridge
  #niter.mm : number of maximum iterations of weighted Lasso iterations for MM and adaptive MM-Lasso
  #niter.S : number of maximum iterations of IWLS for S-Ridge
  #ncores: number of cores to use for parallel computations
  #varsigma: power to elevate the weights for the adaptive MM-Lasso
  #OUTPUT
  #coef.SE: Initial S-Ridge (p+1)-vector of estimated regression parameters, beta(1)=intercept
  #coef.MMLasso: Final MM-Lasso (p+1)-vector of estimated regression parameters, beta(1)=intercept
  #coef.MMLasso.ad: Final adaptive MM-Lasso (p+1)-vector of estimated regression parameters, beta(1)=intercept

  n<-nrow(x)
  p<-ncol(x)

  ###Center and scale covariates and response using median and mad
  prep<-prepara(x,y)
  xnor<-prep$Xnor
  ynor<-prep$ynor
  mux<-prep$mux
  sigx<-prep$sigx
  muy<-prep$muy
  ###

  ###Set-up cluster for parallel computations
  cores<-min(detectCores(),ncores)
  try(cl<-makeCluster(cores))
  try(registerDoParallel(cl))
  locallib<-.libPaths()[[1]]
  parallel::clusterExport(cl, "locallib", envir=environment())
  parallel::clusterEvalQ(cl,.libPaths(locallib))
  ###

  ###Calculate initial estimate and scale
  fit.SE<-sridge(xnor,ynor,normin=0,denormout=0,cualcv.S,numlam.S,niter.S) #Input data is already normalized, output is in normalized data scale
  beta.SE<-fit.SE$coef
  edf.SE<-fit.SE$edf+1
  beta.SE.slo<-beta.SE[2:length(beta.SE)]
  beta.SE.int<-beta.SE[1]
  scale.SE<-fit.SE$scale
  ###

  ###Adjust c1 as recommended by Maronna and Yohai (2010)
  if (edf.SE/n>0.1){
    c1<-3.89 #.90 de efi
  } else{
    c1<-3.44 #.85 de efi
  }
  ###

  ###Calculate candidate lambdas for MMLasso
  lambdamax0<-lambda0(xnor,ynor) #Initial candidate
  lambdamax1<-try(optim.lam(xnor,ynor,beta.SE,scale.SE,lambdamax0,c1,niter.mm))#Improvement
  if(class(lambdamax1)=='try-error'){
    lambdamax1<-lambdamax0
    warning('Using approximate lambdamax for MM-Lasso')
  }
  lambdas<-seq(0,lambdamax1,length.out=numlam.mm)
  if(p>=n){
    lambdas<-lambdas[2:numlam.mm]
  }
  ###

  ####Random permutation and order data for CV
  srt<-sort.int(sample(1:n),index.return=TRUE)
  tt<-srt$x
  orden<-srt$ix
  xnord<-xnor[orden,]
  ynord<-ynor[orden]
  ###

  ###Parallel CV
  exp.mm<-c('MMLasso','CVLasso')
  klam<-NULL
  mse<-foreach(klam=1:length(lambdas),.combine=c,.packages=c('mmlasso','robustHD'),.export=exp.mm)%dopar%{
    CVLasso(xnord,ynord,beta.SE,scale.SE,nfold=cualcv.mm,lambdas[klam],c1,niter.mm)
  }
  if(any(is.infinite(mse))){
    warning('Iteratively weighted Lasso failed when lambda equaled:')
    print(lambdas[which(is.infinite(mse))])
  }
  indmin<-which.min(mse)
  lamin<-lambdas[indmin]
  ###

  #Calculated final MMLasso estimate
  fit.MMLasso<-try(MMLasso(xnor,ynor,beta.SE,scale.SE,lamin,c1,niter.mm))
  if(class(fit.MMLasso)=='try-error'){
    beta.MMLasso.slo<-beta.SE[2:length(beta.SE)]
    beta.MMLasso.int<-beta.SE[1]
    warning('Calculation of the final MM-Lasso estimate failed, returning S-Ridge instead')
  }else{
  beta.MMLasso<-fit.MMLasso$coef
  beta.MMLasso.slo<-beta.MMLasso[2:length(beta.MMLasso)]
  beta.MMLasso.int<-beta.MMLasso[1]}

  activ <- which(beta.MMLasso.slo!=0)
  #If no variables were selected by the MMLasso, return the MMLasso
  if(length(activ)<=1){
    beta.MMLasso.slo<-beta.MMLasso.slo/sigx
    beta.MMLasso.int<-muy+beta.MMLasso.int-mux%*%beta.MMLasso.slo
    beta.SE.slo<-beta.SE.slo/sigx
    beta.SE.int<-muy+beta.SE.int-mux%*%beta.SE.slo
    beta.SE<-c(beta.SE.int,beta.SE.slo)
    beta.MMLasso<-c(beta.MMLasso.int,beta.MMLasso.slo)
    beta.MMLasso.ad<-beta.MMLasso
    re<-list(coef.MMLasso.ad=beta.MMLasso.ad,coef.MMLasso=beta.MMLasso,coef.SE=beta.SE)
    return(re)
  }


  ###Take out covariates not chosen by MMLasso and scale the remaining ones
  w.ad <- (abs(beta.MMLasso.slo[activ]))^varsigma
  xnor.w = as.matrix(scale(xnor[,activ],center=FALSE, scale = 1/w.ad))
  xnord.w <- as.matrix(xnor.w[orden,])
  beta.SE.w <- c(beta.SE[1],beta.SE[activ+1]/w.ad)

  ###Calculate candidate lambdas for adaptive MMLasso
  lambdamax0.ad<-lambda0(xnor.w,ynor)
  lambdamax1.ad<-try(optim.lam(xnor.w,ynor,beta.SE.w,scale.SE,lambdamax0.ad,c1,niter.mm))
  if (class(lambdamax1.ad)=='try-error'){
    lambdamax1.ad<-lambdamax0.ad
    warning('Using approximate lambdamax for adaptive MM-Lasso')
  }
  lambdas.ad<-seq(0,lambdamax1.ad,length.out=numlam.mm)
    if (p>=n){
      lambdas.ad<-lambdas.ad[2:numlam.mm]
    }
  ##

  ###Parallel CV for the adaptive MMLasso
  mse.ad<-foreach(klam=1:length(lambdas.ad),.combine=c,.packages=c('mmlasso','robustHD'),.export=exp.mm)%dopar%{
    CVLasso(xnord.w,ynord,beta.SE.w,scale.SE,nfold=cualcv.mm,lambdas.ad[klam],c1,niter.mm)
  }
  if(any(is.infinite(mse.ad))){
    warning('Iteratively weighted Lasso for adaptative MM-Lasso failed when lambda equaled:')
    print(lambdas.ad[which(is.infinite(mse.ad))])
  }
  indmin.ad<-which.min(mse.ad)
  lamin.ad<-lambdas.ad[indmin.ad]
  ###

  try(stopCluster(cl))

  ###Calculate final estimates and return to original coordinates

  fit.MMLasso.ad<-try(MMLasso(xnor.w,ynor,beta.SE.w,scale.SE,lamin.ad,c1,niter.mm)$coef)
  if(class(fit.MMLasso.ad)=='try-error'){
    beta.MMLasso.slo.ad<-beta.MMLasso.slo
    beta.MMLasso.int.ad<-beta.MMLasso.int
    warning('Calculation of the final adaptive MM-Lasso estimate failed, returning MM-Lasso instead')
  } else{
    beta.MMLasso.slo.ad <- rep(0,p)
    beta.MMLasso.slo.ad[activ] <- fit.MMLasso.ad[2:length(fit.MMLasso.ad)]*w.ad
    beta.MMLasso.int.ad <- fit.MMLasso.ad[1]
  }

  beta.MMLasso.slo.ad<-beta.MMLasso.slo.ad/sigx
  beta.MMLasso.int.ad<-muy+beta.MMLasso.int.ad-mux%*%beta.MMLasso.slo.ad
  beta.MMLasso.ad <- c(beta.MMLasso.int.ad,beta.MMLasso.slo.ad)

  beta.MMLasso.slo<-beta.MMLasso.slo/sigx
  beta.MMLasso.int<-muy+beta.MMLasso.int-mux%*%beta.MMLasso.slo
  beta.MMLasso<-c(beta.MMLasso.int,beta.MMLasso.slo)

  beta.SE.slo<-beta.SE.slo/sigx
  beta.SE.int<-muy+beta.SE.int-mux%*%beta.SE.slo
  beta.SE<-c(beta.SE.int,beta.SE.slo)

  ###

  re<-list(coef.MMLasso.ad=beta.MMLasso.ad,coef.MMLasso=beta.MMLasso,coef.SE=beta.SE)

  return(re)

}


CVLasso<-function(X,y,beta.ini,scale.ini,nfold,lam,c1,niter.mm){
  #Performs nfold-CV for MM-Lasso
  #INPUT
  #X,y: data
  #beta.ini, scale.ini: initial estimate of regression and scale
  #nfold: nfold-CV
  #lam: penalization parameter
  #c1: tuning constant for the rho function
  #niter.mm: maximum number of weighted Lasso iterations for MM-Lasso

  #OUTPUT
  #mse: resulting MSE (estimated using a tau-scale)

  ###Segment data
  n<-nrow(X)
  p<-ncol(X)
  indin<-1:n
  nestim<-n*(1-1/nfold) # Actual number of observations in an estim sample
  lamcv<-lam
  inint<-floor(seq(0,n,length.out=nfold+1))
  resid<-vector(mode='numeric',length=n)
  ###

  for (kk in 1:nfold){
    ###Get test and estimating samples
    testk<-(inint[kk]+1):inint[kk+1]
    estik<-setdiff(indin,testk);
    Xtest<-X[testk,]
    Xesti<-X[estik,]
    ytest<-y[testk]
    yesti<-y[estik]
    fit.MMBR<-try(MMLasso(Xesti,yesti,beta.ini,scale.ini,lambda=lamcv,c1=c1,niter.mm))
    if (class(fit.MMBR)=="try-error"){
      return(Inf)
    }
    beta.MMBR<-fit.MMBR$coef
    beta.MMBR.slo<-beta.MMBR[2:length(beta.MMBR)]
    beta.MMBR.int<-beta.MMBR[1]
    fitk<-Xtest%*%beta.MMBR.slo+beta.MMBR.int
    resid[testk]<-ytest-fitk
  }
  mse<-scale_tau(resid)
  return(mse)
}


MMLasso<-function(xx,y,beta.ini,scale.ini,lambda,c1,niter.mm){
  #Performs iteratively re-weighted Lasso
  #INPUT
  #xx,yy: data
  #beta.ini, scale.ini: initial estimate of regression and scale
  #lambda: penalization parameter
  #c1: tuning constant for the rho function
  #niter.mm: maximum number of iterations
  #OUTPUT:
  #coef: final estimate

  MMcpp_ini<-MMLassoCpp_ini(xx,y,beta.ini)
  x<-MMcpp_ini$x
  resid.n<-MMcpp_ini$resid.n
  tol<-1
  m<-0
  beta.n<-beta.ini
  p<-length(beta.ini)

  while (tol>= 1e-4){
    beta.o<-beta.n
    if(all(beta.o==0)){
      return(list(coef=beta.o))
    }
    MMcpp1<-MMLassoCpp1(x,y,beta.o,scale.ini,c1)
    xort<-MMcpp1$xort
    xjota<-MMcpp1$xjota
    yast<-MMcpp1$yast
    alpha<-MMcpp1$alpha
    u.Gram<- !(p>=500)
    try(res.Lasso <- robustHD:::fastLasso(xort, yast,lambda,intercept=FALSE, normalize=FALSE, use.Gram=u.Gram),silent=TRUE)
    if (class(res.Lasso)=="try-error"){
      warning("fastLasso failed")
      return(list(coef=beta.o))
    }
    beta.Lasso <- coef(res.Lasso)
    MMcpp2<-MMLassoCpp2(xjota,yast,beta.Lasso,beta.o,alpha)
    beta.n<-MMcpp2$beta.n
    tol<-MMcpp2$tol
    m<-m+1
    if (m >= niter.mm) {tol<-0}
  }
  return(list(coef=beta.n))
}



optim.lam<-function(x,y,beta.ini,scale.ini,lambdamax,c1,niter.mm){
  #Finds (approximately) smallest lambda such that the slope of the MM-Lasso estimate is zero
  #INPUT
  #X, y: data set
  #beta.ini, scale.ini: initial estimates and scale
  #lambdamax: initial guess for lambda
  #c1: tuning constant for the rho function
  #niter.mm: maximum number of iterations
  #OUTPUT
  #lambda: smallest lambda such that the slope of the MM-Lasso estimate is zero

  n<-nrow(x)
  p<-ncol(x)
  lambda2<-lambdamax
  beta.n<-MMLasso(x,y,beta.ini,scale.ini,lambda2,c1,niter.mm)$coef
  zeros.n<-sum(beta.n[2:(p+1)]==0)
  while (zeros.n<p & lambda2<max(n,1e4)){
    lambda2<-2*lambda2
    beta.n<-MMLasso(x,y,beta.ini,scale.ini,lambda2,c1,niter.mm)$coef
    zeros.n<-sum(beta.n[2:(p+1)]==0)
  }
  if(lambda2>=max(n,1e4))
  {
    return(lambda2)
  }
  lambda1<-lambdamax/2
  beta.o<-MMLasso(x,y,beta.ini,scale.ini,lambda1,c1,niter.mm)$coef
  zeros.o<-sum(beta.o[2:(p+1)]==0)
  while(zeros.o==p){
    lambda1<-lambda1/2
    beta.o<-MMLasso(x,y,beta.ini,scale.ini,lambda1,c1,niter.mm)$coef
    zeros.o<-sum(beta.o[2:(p+1)]==0)
  }
  for (j in 1:5){
    lambda<-0.5*(lambda1+lambda2)
    beta.n<-MMLasso(x,y,beta.ini,scale.ini,lambda,c1,niter.mm)$coef
    zeros.n<-sum(beta.n[2:(p+1)]==0)
    if (zeros.n<p){
      lambda1<-lambda
    }else{
      lambda2<-lambda
    }
  }
  return(lambda2)
}



prepara<-function(X,y){
  #Centers y and the columns of X to zero median, and normalizes X to unit mad
  #INPUT
  #X, y: data set
  #OUTPUT
  #Xnor: centered normalized X
  #ycen: centered y
  #mux: median vector of X
  #scalex: mad vector of X
  #muy: median of y
  #sigy:scale of y
  muy<-median(y)
  ycen<-y-muy
  sigy<-mad(y)
  mux<-apply(X,2,'median')
  sigx<-apply(X,2,'mad')
  Xnor<-scale(X,center=mux,scale=sigx)
  res<-list(Xnor=Xnor,ynor=ycen,mux=mux,sigx=sigx,muy=muy,sigy=sigy)
  return(res)
}
