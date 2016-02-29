#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec rho_marina(arma::vec x,double cw) {
  int n = x.size();
  arma::vec out(n);
  for(int i = 0; i < n; ++i) {
    if (fabs(x[i])>=cw)
      out[i]=1;
    else 
      out[i]=1-pow(1-pow(x[i]/cw,2),3);
  }
  return out;
}

// [[Rcpp::export]]
arma::vec rho_bisq(arma::vec x,double cw) {
  int n = x.size();
  arma::vec out(n);
  for(int i = 0; i < n; ++i) {
    if (fabs(x[i])>=cw)
      out[i]=1;
    else 
      out[i]=1-pow(1-pow(x[i]/cw,2),3);
  }
  return out;
}


// [[Rcpp::export]]
arma::vec psi_marina(arma::vec x, double cw){
  int n = x.size();
  arma::vec out(n);
  for(int i = 0; i < n; ++i) {
    if (fabs(x[i])>=cw)
      out[i]=0;
    else 
      out[i]=6*(x[i]/pow(cw,2))*pow(1-pow(x[i]/cw,2),2);
  }
  return out;  
}

// [[Rcpp::export]]
arma::vec weight_marina(arma::vec x,double cw){
  int n=x.size();
  arma::vec out(n);
  for(int i = 0; i < n; ++i) {
    if (x[i]==0)
      out[i]=1;
    else if(fabs(x[i])<=cw){ 
      out[i]=pow(1-pow(x[i]/cw,2),2);}
    else if (fabs(x[i])>cw){
      out[i]=0;}
  }
  return out;
}


// [[Rcpp::export]]
double Mscale_mar(arma::vec x,double b, double cc) {
  int n = x.size();
  double sc = median(abs(x)) / 0.6745;
  double sc2 = 0;
  int maxit = 100;
  int i = 0;
  double eps = 1e-6;
  double err = 1;
  while (i<maxit and err>eps){
    i+=1;
    sc2 = sc * sqrt( (sum( rho_marina( x/sc , cc) ) / b   )/ n );
    err = fabs(sc2/sc - 1);
    sc = sc2;
  }
  return sc;
}


double Mscale_bisq(arma::vec x,double b, double cc,int cuad) {
  int n = x.size();
  double sc = median(abs(x)) / 0.6745;
  double sc2 = 0;
  int maxit = 100;
  int i = 0;
  double eps = 1e-6;
  double err = 1;
  while (i<maxit and err>eps){
    i+=1;
    sc2 = sc * sqrt( (sum( rho_bisq( x/sc , cc) ) / b   )/ n );
    err = fabs(sc2/sc - 1);
    sc = sc2;
  }
  if (cuad==1) {return pow(sc,2);}
  return sc;
}

// [[Rcpp::export]]
double scale_tau(arma::vec x){
    int n = x.n_elem;
    double c2 = 3;
    arma::vec x_dot = abs(x);
    double sigma0 = median(x_dot);
    x = x_dot/(sigma0);
    arma::vec rho = pow(x,2);
    double c22 = pow(c2,2);
    for (int j=0; j<n; j++){
      if(rho[j]>c22){rho[j] = c22;}
    }
    double consis = 1;
    double out = sigma0*sqrt(sum(rho)/(n*consis));
    return out;
}


arma::vec spa_med(arma::mat x){
  //Spatial M-median , w=weights 
  int n = x.n_rows;
  double del0 = 1.e-6;
  double deldis = 0;
  arma::vec w = vec(n);
  w = sqrt(sum(pow(x,2),1));
  deldis = del0*median(w);
  for (int j=0; j<n; j++){
    if(w[j]<deldis){w[j]=deldis;}
  }
  w = 1/w;
  w = w/sum(w);
  return w;
}


arma::vec rob_sq(arma::mat z){
  //Squared dispersion
  int p = z.n_cols;
  arma::rowvec mume = median(z);
  arma::vec out = vec(p);
  for (int j=0; j<p ; j++){
    out[j] = Mscale_bisq(z.col(j)-mume[j],0.5,1.54764,1);
  }
  return out;
}

// [[Rcpp::export]]
List my_svdecon(arma::mat x){
  int n = x.n_rows;
  int p = x.n_cols;
  arma::mat svd_Ux = mat(n,n);
  arma::mat svd_Vx = mat(p,p);
  arma::vec svd_sigmax = vec(p);
  arma::vec svd_sigmaxa = vec(n);
  uvec tes;
  double eps;
  svd(svd_Ux,svd_sigmax,svd_Vx,x);
  if (p>n){
    svd_Vx.shed_cols(n,p-1);
    svd_sigmaxa = svd_sigmax.rows(0,n-1);
    eps = 1e-8*svd_sigmaxa(1)*p;
    tes = find(svd_sigmaxa>eps);
  } else{
    eps = 1e-8*svd_sigmax(1)*n;
    tes = find(svd_sigmax>eps);
  }
  svd_Ux = svd_Ux.cols(tes);
  svd_Vx = svd_Vx.cols(tes);
  List ret;
  ret["U"] = svd_Ux;
  ret["V"] = svd_Vx;
  return ret;
}

// [[Rcpp::export]]
List SPCC(arma::mat x){
//  SPC Spherical Principal Components (Locantore et al., 1999)
//  lamda= Robust "eigenvalues" (increasing); b=Matrix of eigenvectors
//  mu=spatial median; scores=projection of x (centered) on eigenvectors
  arma::vec w = spa_med(x);
  arma::mat y = diagmat(w)*x;
  List svd_y = my_svdecon(y);
  arma::mat svd_Vy = svd_y["V"];
  arma::mat scores = x*svd_Vy;
  arma::vec lamda = rob_sq(scores);
  arma::uvec index = sort_index(lamda);
  arma::vec lamda_ord = sort(lamda);
  svd_Vy = svd_Vy.cols(index);
  scores = scores.cols(index);
  List ret;
  ret["lamda"] = lamda_ord;
  ret["b"] = svd_Vy;
  ret["scores"] = scores;
  return ret;
}



// [[Rcpp::export]]
List MMLassoCpp_ini(arma::mat xx,arma::vec y, arma::vec beta_ini){
  int n=xx.n_rows;
  int p=xx.n_cols;
  arma::mat x = mat(n,0);
  List ret;
  x.insert_cols(0, ones(n,1)); 
  x.insert_cols(1,xx) ;
  arma::vec resid_n = y - x*beta_ini;
  ret["x"] = x;
  ret["resid.n"] = resid_n;
  ret["nrow"] = n;
  ret["ncol"] = p;
  return ret;
}

// [[Rcpp::export]]
List MMLassoCpp1(arma::mat x,arma::vec y, arma::vec beta_ini, double scale_ini,double c1){
  int n=x.n_rows;
  int p=x.n_cols-1;
  arma::mat xast = mat(n,p+1);
  arma::mat xastrun = mat(n,p+1);
  arma::mat xort = mat(n,p);
  arma::mat xjotamat = mat(n,p);
  arma::vec xjota = vec(n);
  arma::vec yast = vec(n);
  arma::vec alpha = vec(p);
  arma::vec w=vec(n);
  arma::vec resid_n = y - x*beta_ini;
  List ret;
  w = sqrt(weight_marina(resid_n/scale_ini,c1));
  xast = diagmat(w)*x;
  yast = diagmat(w)*y;
  xjota = xast.col(0);
  xastrun = xast;
  xastrun.shed_col(0);
  if (sum(abs(xjota))>0){
    alpha = xastrun.t()*xjota/dot(xjota,xjota);
  } else{
    alpha= zeros(p);
  }
  xjotamat = repmat(xjota,1,p);
  xort = xastrun - xjotamat*diagmat(alpha);
  ret["xort"] = xort;
  ret["xjota"] = xjota;
  ret["yast"] = yast;
  ret["alpha"] = alpha;
  return ret;
}


// [[Rcpp::export]]
List MMLassoCpp2(arma::vec xjota,arma::vec yast, arma::vec beta_lars, arma::vec beta_o,arma::vec alpha){
  int p = beta_lars.n_elem;
  double beta_n_int = 0;
  if (sum(abs(xjota))>0){
    beta_n_int = dot(xjota,yast)/dot(xjota,xjota)-dot(beta_lars,alpha);
  } else{
    beta_n_int= 0;
  }
  
  arma::vec beta_n = vec(p+1);
  arma::vec u = vec(p+1);
  double tol;
  beta_n[0] = beta_n_int;
  for(int i = 1; i <= p; ++i) {
    beta_n[i] = beta_lars[i-1];
  }
  u = beta_n-beta_o;
  tol = sqrt(dot(u,u)/dot(beta_o,beta_o));
  List ret;
  ret["beta.n"] = beta_n;
  ret["tol"] = tol;
  return ret;
}


List desrobrid(arma::mat x, arma::vec y,int niter,double lam,double betinte,arma::vec betinslo,double cc,double delsca,double epsilon){
//  IWLS descent for S-Ridge
//  INPUT
//  X, y: data
//  lam: penalization parameter
//  betinte: betinslo: initial intercept and slope estimates
//  cc: tuning constant for the M-scale
//  delsca: delta for the M-scale
//  epsilon: tolerance parameter for converge
//  OUTPUT:
//  beta: final estimate
//  edf: estimated EDF
//  obj: value of the objective function at the optimum
  int n = x.n_rows;
  int p = x.n_cols;
  arma::vec res0 = y-x*betinslo;
  res0 = res0-betinte;
  arma::vec tt = vec(n);
  arma::vec w = vec(n);
  arma::vec rw = vec(n);
  arma::vec ycen = vec(n);
  arma::vec yw = vec(n);
  arma::vec yau = vec(n+p);
  arma::vec beta = vec(p);
  arma::vec beta_fin = vec(p+1);
  arma::vec resin = vec(n);
  arma::vec res = vec(n);
  arma::vec edf_pre = vec(n+p);
  arma::mat xw = mat(n,p);
  arma::mat xau = mat(n+p,p);
  arma::mat rxau = mat(p,p);
  arma::mat Q_xau = mat(n+p,p);
  arma::mat R_xau = mat(p,p);
  int iter=0;
  double lala ;
  double delta = 1;
  double conve = 1;
  double sig = 0;
  double sig0 = Mscale_mar(res0,delsca,cc);
  double crit0 = n*pow(sig0,2)+lam*dot(betinslo,betinslo);
  double crit = 0;
  double binter = betinte;
  double tolcrit = epsilon;
  double tolres = epsilon;
  double edf = 0;
  List ret;
  while (iter<niter and (delta>tolcrit or conve>tolres)){
    iter+=1;
    tt = res0/sig0;
    w = weight_marina(tt,cc);
    rw = sqrt(w);
    ycen = y-binter;
    xw = diagmat(rw)*x;
    yw = diagmat(rw)*ycen;
    lala = mean(w%pow(tt,2))*lam;
    xau = xw;
    rxau = sqrt(lala)*eye(p,p);
    xau.insert_rows(n,rxau);
    yau = yw ;
    yau.insert_rows(n,zeros(p));
    qr_econ(Q_xau,R_xau,xau);
    beta = solve(R_xau,Q_xau.t()*yau);
    resin = y-x*beta;
    if (sum(w)>0){
      binter = sum(resin%w)/sum(w);
      }else{
        binter = median(resin);
      }
    res =resin-binter ;
    sig = Mscale_mar(res,delsca,cc);
    crit = n*pow(sig,2)+lam*dot(beta,beta);
    delta = 1-crit/crit0;
    conve = max(abs(res-res0))/sig;
    res0 = res;
    sig0 = sig;
    crit0 = crit;
 }
  beta_fin[0] = binter;
  for(int i = 1; i <= p; ++i) {
    beta_fin[i] = beta[i-1];
  } 
  edf_pre = sum (pow(Q_xau,2),1);
  for(int j = 1; j <= n; ++j) {
    edf += edf_pre[j-1];
  }

  ret["edf"] = edf;
  ret["beta"] = beta_fin;
  ret["objF"] = crit;
  return ret;
}


List regrid(arma::mat x,arma:: vec y,double lambda){
  int n = x.n_rows;
  int p = x.n_cols;
  arma::mat xau = mat(n+p,p+1);
  arma::mat xext = mat(p,p+1);
  arma::mat rxau = mat(p,p);
  arma::vec yau = vec(n+p);
  arma::vec beta = vec(p+1);
  List ret;
  x.insert_cols(0, ones(n,1)); 
  xext = sqrt(lambda)*eye(p,p);
  xext.insert_cols(0,zeros(p));
  xau = join_cols(x,xext);
  yau = join_cols(y, zeros(p));
  beta = solve(xau,yau);
  ret["beta"] = beta;
  ret["Xau"] = xau;
  ret["yau"] = yau;
  return ret;
}

// [[Rcpp::export]]
List rr_se(arma::mat X, arma::vec y,double lambda2, double deltaesc, double cc_scale, int nkeep, int niter, double epsilon){
//  Initial estimator: adapted Pena-Yohai (P&Y) for Ridge Regression (details in Maronna, Technometrics 2011)
//  3*(p+1)+1 initial LS-Ridge estimators based on original and clean samples are computed
//  The estimate (or nkeep estimates) that minimizes the S-Ridge objective function is used as the initial estimator(s). 
//  INPUT
//  X,y: data
//  lambda2:l2-penalty
//  deltaesc:delta for initial M-scale estimator
//  nkeep:number of initial candidates selected, Obs: hard-coded nkeep=5 for now
//  niter:number of IWLS iterations
//  prosac: proportion of observations removed based on PSC(=deltaesc by default)
//  epsilon: effective zero
//  OUTPUT
//  coef: coefficients of the estimated regression vector
//  edf: estimated EDF
//  obj: value of the objective function at the optimum
  int n = X.n_rows;
  int p = X.n_cols;
  int m = floor(deltaesc*n);
  arma::mat Xuno = mat(n,p+1); 
  Xuno = join_rows(ones(n,1),X);
  arma::vec critkeep = vec(3*p+4);
  arma::mat betakeep = mat(p+1,3*p+4);
  double lamLS_2 = lambda2;  
  //corrects lambda to match RR-LS. Appendix Maronna (2012)
  //double factr = 12.68451+(-46.97558)*(deltaesc-0.25)+91.72468*(deltaesc-0.25)*(deltaesc-0.35);
  double factr = 23.9716 -73.4391*deltaesc+ 64.9480*pow(deltaesc,2);
  double lam0_2 = lambda2*factr;
  //Computes RR-LS finding the OLS on the extended X matrix and extended Y vector
  List ext_model = regrid(X,y,lamLS_2);
  arma::mat Xau = ext_model["Xau"];
  arma::vec yau = ext_model["yau"];
  arma::vec betaj = ext_model["beta"];
  arma::vec betaslo = vec(p);
  arma::vec resj = yau-Xau*betaj;
  arma::vec restrun = vec(n);
  
  //Singular Value Decomposition to obtain the H matrix
  List svd_Xau = my_svdecon(Xau);
  arma::mat svd_U = svd_Xau["U"];
  arma::mat svd_V = svd_Xau["V"];
  arma::mat UV = mat(n+p,p+1);
  arma::mat Q = mat(p+1,p+1);
  arma::mat eigvec_Q = mat(p+1,p+1);
  arma::mat eigvec_flip_Q = mat(p+1,p+1);
  arma::mat Z = mat(n+p ,p+1);
  arma::vec h = vec(n+p);
  arma::vec w = vec(n+p);
  arma::vec eigval_Q = vec(p+1);
  double sigj = 0;
  h = sum(pow(svd_U,2),1);
  w = pow(resj / (1-h),2);
  UV = svd_U*svd_V.t();
  Q = UV.t() * diagmat(w) * UV;
  eig_sym(eigval_Q,eigvec_Q,Q);
  eigvec_flip_Q = fliplr(eigvec_Q);
  Z = UV * eigvec_flip_Q;
  Z.shed_rows(n,n+p-1);
  //each columns of Z represent the forecast change for each observation in the direction of u_j (col of U)
  restrun = resj.rows(0,n-1);
  betaslo = betaj.rows(1,p);
  //Compute the M-scale for the LS-Ridge and value of objective function to be minimized by RR-SE 
  sigj = Mscale_mar(restrun,deltaesc,cc_scale);
  critkeep[0] = n*pow(sigj,2)+lam0_2*dot(betaslo,betaslo);
  betakeep.col(0) = betaj;

  
  int n1 = n-m; //new sample size
  double lam2 = lamLS_2*n1/n; //corrects lambda for LS-Ridge to acknowledge new sample size
  arma::mat Xord = mat(n,p);
  arma::mat Xj = mat(n1,p);
  arma::vec yord = vec(n);
  arma::vec yj = vec(n1);
  arma::vec betaj2 = vec(p+1);
  arma::vec betaj3 = vec(p+1);
  arma::uvec ordt = uvec(n);
  int l = 1;
  List listaj;
  
  for(int j = 0; j <=p ; ++j) {
      //Remove m observations with the highest forcast change (i.e., m-highest z_ij, order(zj) the sort is increasing) 
      ordt = sort_index(Z.col(j));    
      Xord = X.rows(ordt);
      yord = y.rows(ordt);
      Xj = Xord.rows(0,n1-1);
      yj = yord.rows(0,n1-1);
      listaj = regrid(Xj,yj,lam2);
      arma::vec betaj = listaj["beta"];
      betaslo = betaj.rows(1,p);
      resj = y-Xuno*betaj; //residuals of ALL observations using betaj based on clean sample
      sigj = Mscale_mar(resj,deltaesc,cc_scale);
      critkeep[l] = n*pow(sigj,2)+lam0_2*dot(betaslo,betaslo);
      betakeep.col(l) = betaj;
      
      //Remove m observations with the lowest forcast change (i.e., m-highest z_ij, order(zj) the sort is decreasing) 
      ordt = sort_index(Z.col(j),"descend");    
      Xord = X.rows(ordt);
      yord = y.rows(ordt);
      Xj = Xord.rows(0,n1-1);
      yj = yord.rows(0,n1-1);
      listaj = regrid(Xj,yj,lam2);
      arma::vec betaj2 = listaj["beta"];
      betaslo = betaj2.rows(1,p);
      resj = y-Xuno*betaj2; //residuals of ALL observations using betaj based on clean sample
      sigj = Mscale_mar(resj,deltaesc,cc_scale);
      critkeep[l+1] = n*pow(sigj,2)+lam0_2*dot(betaslo,betaslo);
      betakeep.col(l+1) = betaj2;
      
      //Remove m observations with the highest absolute forcast change (large absolute values) 
      ordt = sort_index(abs(Z.col(j)));    
      Xord = X.rows(ordt);
      yord = y.rows(ordt);
      Xj = Xord.rows(0,n1-1);
      yj = yord.rows(0,n1-1);
      listaj = regrid(Xj,yj,lam2);
      arma::vec betaj3 = listaj["beta"];
      betaslo = betaj3.rows(1,p);
      resj = y-Xuno*betaj3; //residuals of ALL observations using betaj based on clean sample
      sigj = Mscale_mar(resj,deltaesc,cc_scale);
      critkeep[l+2] = n*pow(sigj,2)+lam0_2*dot(betaslo,betaslo);
      betakeep.col(l+2) = betaj3;
      l += 3;
  }
  
  //Find the nkeep best initial estimators based to perform full iterations of IWLS #Esto es un poco como Croux
  
  arma::uvec ii = sort_index(critkeep);
  betakeep = betakeep.cols(ii);
  

  
  
//For each candidate initial estimator, we iterate weigthed normal equations niter times until a local minimum of RR-SE objective function is found
  arma::vec crit_fin = vec(5);
  arma::mat beta_fin = mat(p+1,5);
  arma::vec edf_fin = vec(5);
  List des_se;
  for (int k=0 ; k<=4; k++){
    des_se = desrobrid(X,y,niter,lam0_2,betakeep(0,k),betakeep(span(1,p),k),cc_scale,deltaesc,epsilon);
    crit_fin[k] = des_se["objF"];
    arma::vec betaj = des_se["beta"];
    beta_fin.col(k) = betaj;
    edf_fin[k] = des_se["edf"];
  }

//Output: the coefficients, critical values and estimated EDF
  uword opt_index;
  crit_fin.min(opt_index);
  double edf = edf_fin[opt_index];
  double objF = crit_fin[opt_index];
  arma::vec coef = beta_fin.col(opt_index);
  List ret;
  ret["coef"] = coef;
  ret["edf"] = edf;
  ret["objF"] = objF;
  return ret;
}
