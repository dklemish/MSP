#include <R.h>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <string>
#include <cmath>

// #define PI 3.141592653589793238462643383279502884197169399

// RcppArmadillo.h also includes Rcpp.h

// Declare dependency on RcppArmadillo so Rcpp knows to link libraries
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/***** Function headers *****/
void dmvnorm(const mat &, const vec &, const mat &, vec &, int, int);
void rmvnorm(vec&, mat &, vec &, int);
void rgev(const double &, const int &, vec &);
void rpstable(const double &, const int &, vec &);
void calcDistance(mat &, mat &, int);
void storeCalc(mat &, vec &, vec &, double &, int &, int &, int &);
void storeCalcSmith(mat &, mat &, vec &, vec &, vec &, double &, int &, int &, int &, int &);

void covMatern(int, mat &, double, double, double, mat &);
void covPowExp(int, mat &, double, double, double, mat &);

void mspBrownResnick(bool, int, int &, int &, double, vec &, const vec &, const mat &, mat &, vec &);
void mspExtremeT(bool, int, int &, int &, double cstar, vec &, const vec &, const mat &, const double, mat &, vec &);

// [[Rcpp::export]]
SEXP Sim_SP(
    SEXP x_, 
    SEXP cstar_,
    SEXP model_,
    SEXP corrFn_,
    SEXP keepPsi_,
    SEXP covPar_, 
    SEXP alpha_){
  /***** Initialization *****/
  // Variable declarations
  bool keepPsi;
  int n;                // Total # of grid points (n in 1D, n^2 in 2D)
  int d;                // # of dimensions
  int count = 0;
  int curr_storage = 100;
  double cstar;         // Approx max value of stochastic process
  double sigma2;        // Partial sill
  double phi;           // 1 / range covariance parameter
  double nu;            // Matern parameter
  double gamma;         // Exponential covariance parameter
  double alpha;         // Extreme T process parameter
  std::string corrFn;
  std::string model;
  mat x;                // Location matrix
  mat Sigma;            // Covariance matrix of underlying process
  mat distMat;          // Matrix of pairwise distances between spatial locations
  mat psiStore;         // Matrix to store psi functions
  vec Z;                // Output
  vec zero_mean;        // Zero vector, size n
  vec zetaStore;        // Matrix to store zeta values
  List covPar;          // List of covariance function parameters
  
  // Convert SEXP objects to Rcpp objects
  NumericMatrix X = as<NumericMatrix>(x_);
  cstar = as<double>(cstar_);
  corrFn = as<std::string>(corrFn_);
  model = as<std::string>(model_);
  keepPsi = as<bool>(keepPsi_);
  covPar = as<List>(covPar_);
  alpha = as<double>(alpha_);
  
  // Initialize other variables
  x   = mat(X.begin(), X.nrow(), X.ncol(), TRUE);
  n   = X.nrow();
  d   = X.ncol();
  
  sigma2 = as<double>(covPar[0]);
  phi    = as<double>(covPar[1]);
  nu     = as<double>(covPar[2]);
  gamma  = as<double>(covPar[3]);
  
  Sigma   = mat(n, n);
  distMat = mat(n, n);
  
  zero_mean = vec(n, fill::zeros);
  
  // Allocate output storage
  Z = vec(n, fill::zeros);
  psiStore = mat(curr_storage, n);
  zetaStore = vec(curr_storage);
  
  /**** Calculations ****/
  /* Calculate distance matrix between all points */
  calcDistance(x, distMat, n);
  
  if(corrFn == "Matern"){
    covMatern(n, distMat, sigma2, phi, nu, Sigma);
  } else if(corrFn == "Exp"){
    covPowExp(n, distMat, sigma2, phi, 1.0, Sigma);
  } else if(corrFn == "PowExp"){
    covPowExp(n, distMat, sigma2, phi, gamma, Sigma);
  }
  GetRNGstate();
  
  if(model=="Brown-Resnick"){
    mspBrownResnick(keepPsi, n, count, curr_storage, cstar, Z, zero_mean, Sigma, psiStore, zetaStore);
  } else if(model=="extremeGauss"){
    mspExtremeT(keepPsi, n, count, curr_storage, cstar, Z, zero_mean, Sigma, 1.0, psiStore, zetaStore);
  } else if(model=="extremeT"){
    mspExtremeT(keepPsi, n, count, curr_storage, cstar, Z, zero_mean, Sigma, alpha, psiStore, zetaStore);
  }
  
  if(keepPsi){
    mat psiOut(psiStore.rows(0, count-1));
    vec zetaOut(zetaStore.subvec(0, count-1));
    
    return Rcpp::List::create(
      Rcpp::Named("dim") = d,
      Rcpp::Named("num_functions") = count-1,
      Rcpp::Named("X") = x,
      Rcpp::Named("Z") = Z,
      Rcpp::Named("psi") = psiOut,
      Rcpp::Named("zeta") = zetaOut);
  }else{
    return Rcpp::List::create(
      Rcpp::Named("dim") = d,
      Rcpp::Named("num_functions") = count-1,
      Rcpp::Named("X") = x,
      Rcpp::Named("Z") = Z);
  }
}

// [[Rcpp::export]]
SEXP Sim_Reich(
    SEXP x_, 
    SEXP knots_,
    SEXP Sigma_, 
    SEXP keepPsi_,
    SEXP bw_,
    SEXP alpha_){
  /***** Initialization *****/
  // Variable declarations
  bool keepPsi;
  int n;                // Total # of grid points (n in 1D, n^2 in 2D)
  int nknots;           // Total # of knot points
  int d;                // # of dimensions
  int i, j;
  double bw;
  //double alpha;         // 
  double temp;
  
  mat x;                // Location matrix
  mat knots;            // Knot locations
  mat Sigma;            // Covariance matrix of underlying process
  mat kernelVal;        // Matrix to store kernel values at each grid point
  mat alpha;            // Spatial nugget effect
  vec A;                // Vector to store A values
  vec U;                // Vector to store U values
  vec theta;            // Vector to store theta values
  vec Z;                // Output
  
  // Convert SEXP objects to Rcpp objects
  NumericMatrix X = as<NumericMatrix>(x_);
  NumericMatrix knotsR = as<NumericMatrix>(knots_);
  NumericMatrix SigmaR = as<NumericMatrix>(Sigma_);
  NumericMatrix alp = as<NumericMatrix>(alpha_);
  
  keepPsi = as<bool>(keepPsi_);
  bw      = as<double>(bw_);
  //alpha   = as<double>(alpha_);
  
  // Initialize other variables
  x     = mat(X.begin(), X.nrow(), X.ncol(), TRUE);
  knots = mat(knotsR.begin(), knotsR.nrow(), knotsR.ncol(), TRUE);
  Sigma = mat(SigmaR.begin(), SigmaR.nrow(), SigmaR.ncol(), TRUE);
  alpha = mat(alp.begin(), alp.nrow(), alp.ncol(), TRUE);
  
  n   = X.nrow();
  d   = X.ncol();
  nknots = knotsR.nrow();
  kernelVal = mat(n, nknots);
  
  // Allocate output storage
  Z = vec(n, fill::zeros);
  A = vec(nknots);
  U = vec(n);
  theta = vec(n);
  
  /**** Calculations ****/
  // Calculate kernel values at each grid point
  
  const double normConst1 = 1/(2*PI*bw);
  const double normConst2 = 1/(2*bw);
  vec kernelSum = vec(n);
  
  for(i=0; i<n; i++){
    for(j=0; j<nknots; j++){
      kernelVal.at(i,j) = normConst1 * exp(-normConst2 * as_scalar((x.row(i) - knots.row(j)) * Sigma * (x.row(i) - knots.row(j)).t()));
    }
    kernelSum.at(i) = sum(kernelVal.row(i));
  }
  
  for(i=0; i<n; i++){
    for(j=0; j<nknots; j++){
      kernelVal.at(i,j) = kernelVal.at(i,j) / kernelSum.at(i);
    }
  }
  
  // Simulations 
  GetRNGstate();
  
  for(i=0; i < n; i++){
    rgev(alpha.at(i,3), n, U);          // Simulate values for U (nugget effect)    
  }
  for(i = n; i < (n+nknots); i++){
    rpstable(alpha.at(i,3), nknots, A); // Simulate values for A/theta (basis function coefficients)
  }

  for(i = 0; i<n; i++){
    temp = 0;
    for(j = 0; j<nknots; j++){
      temp += A.at(j) * pow(kernelVal.at(i,j), 1/alpha);
    }
    theta.at(i) = pow(temp, alpha);
  }
  
  Z = U % theta;
  
  if(keepPsi){
    return Rcpp::List::create(
      Rcpp::Named("dim") = d,
      Rcpp::Named("X") = x,
      Rcpp::Named("Z") = Z,
      Rcpp::Named("A") = A,
      Rcpp::Named("U") = U,
      Rcpp::Named("theta") = theta,
      Rcpp::Named("kernel") = kernelVal);
  }else{
    return Rcpp::List::create(
      Rcpp::Named("dim") = d,
      Rcpp::Named("X") = x,
      Rcpp::Named("Z") = Z);
  }
}

// [[Rcpp::export]]
SEXP Sim_Smith(
    SEXP x_, 
    SEXP cstar_,
    SEXP keepPsi_,
    SEXP Sigma_, 
    SEXP radius_){
  /***** Initialization *****/
  // Variable declarations
  bool keepPsi;
  int n;                // Total # of grid points (n in 1D, n^2 in 2D)
  int d;                // # of dimensions
  int count = 0;
  int curr_storage = 100;
  double cstar;         // Approx max value of stochastic process
  double radius;
  double norm_range;
  mat x;                // Location matrix
  mat Sigma;            // Covariance matrix of underlying process
  mat SStore;           // Matrix to store storm locations
  mat psiStore;         // Matrix to store psi functions
  vec S;                // Storm locations
  vec unitVec;          // Unit vector in direction of storm location from center of X
  vec psi;              // 
  vec Z;                // Output
  vec center;           // Center of X space
  vec zetaStore;        // Matrix to store zeta values
  
  // Convert SEXP objects to Rcpp objects
  NumericMatrix X = as<NumericMatrix>(x_);
  NumericMatrix Sig = as<NumericMatrix>(Sigma_);
  
  cstar   = as<double>(cstar_);
  radius  = as<double>(radius_);
  keepPsi = as<bool>(keepPsi_);
  
  // Initialize other variables
  x   = mat(X.begin(), X.nrow(), X.ncol(), TRUE);
  Sigma = mat(Sig.begin(), Sig.nrow(), Sig.ncol(), TRUE);
  n   = X.nrow();
  d   = X.ncol();
  
  center = vec(d);
  unitVec = vec(d);
  center.fill(0.5);
  S = vec(d);
  
  // Allocate output storage
  psi = vec(n, fill::zeros);
  Z = vec(n, fill::zeros);
  
  zetaStore = vec(curr_storage);
  psiStore  = mat(curr_storage, n);
  SStore = mat(curr_storage, d);
  
  /**** Calculations ****/

  double minZ = 0;
  double zetaInv = R::rexp(1);
  
  if(d==2){
    double theta;
    
    norm_range = PI * radius * radius;
    
    while(cstar/zetaInv > minZ/norm_range){
      theta = R::runif(0.0, 2*PI);
      unitVec[0] = cos(theta);
      unitVec[1] = sin(theta);
      
      S = center + R::runif(0.0, radius)*unitVec;
      dmvnorm(x, S, Sigma, psi, n, d);
      
      Z = max(Z, norm_range * psi / zetaInv);
      
      zetaInv += R::rexp(1);
      minZ = min(Z);
      
      if(keepPsi){
        storeCalcSmith(psiStore, SStore, zetaStore, psi, S, zetaInv, n, d, count, curr_storage);
      }
    }
  }else{
    norm_range =  2 * radius;
    while(cstar/zetaInv > minZ/norm_range){
      
      S = center + R::runif(-radius, radius);
      dmvnorm(x, S, Sigma, psi, n, d);
      Z = max(Z, norm_range * psi / zetaInv);
      
      zetaInv += R::rexp(1);
      minZ = min(Z);
      
      if(keepPsi){
        storeCalcSmith(psiStore, SStore, zetaStore, psi, S, zetaInv, n, d, count, curr_storage);
      }
    }
  }
  
  GetRNGstate();
  
  if(keepPsi){
    mat psiOut(psiStore.rows(0, count-1));
    mat SOut(SStore.rows(0, count-1));
    vec zetaOut(zetaStore.subvec(0, count-1));
    
    return Rcpp::List::create(
      Rcpp::Named("dim") = d,
      Rcpp::Named("num_functions") = count-1,
      Rcpp::Named("X") = x,
      Rcpp::Named("Z") = Z,
      Rcpp::Named("psi") = psiOut,
      Rcpp::Named("S") = SOut,
      Rcpp::Named("zeta") = zetaOut);
  }else{
    return Rcpp::List::create(
      Rcpp::Named("dim") = d,
      Rcpp::Named("num_functions") = count-1,
      Rcpp::Named("X") = x,
      Rcpp::Named("Z") = Z);
  }
}

void dmvnorm(const mat &x, const vec &mu, const mat &Sigma, vec &result, int n, int d){

  mat SigmaInv = inv_sympd(Sigma);
  const double norm_fac =  pow(det(Sigma), -0.5) * pow(2.0*PI, (double)-d/2);
  
  for(int i=0; i < n; i++){
    result.at(i) = as_scalar(-0.5*(x.row(i) - mu.t()) * SigmaInv * (x.row(i).t() - mu));
  }
  
  result = norm_fac * exp(result);
  
  return;
}

void rmvnorm(const vec &mu, const mat &Sigma, vec &result, int d){
  mat L = arma::chol(Sigma, "lower");
  
  vec draw = vec(d);
  for(int i=0; i < d; i++){
    draw.at(i) = R::rnorm(0, 1);
  }
  
  result = mu + L * draw;
  
  return;
}

void rgev(const double &alpha, const int &n, vec & result){
  // Simulate from a GEV(1, alpha, alpha) distribution
  double u;
  
  for(int i=0; i < n; i++){
    u = R::runif(0, 1);
    result.at(i) = pow(-log(u),-alpha);
  }
  
  return;
  
}

void rpstable(const double &alpha, const int &nknots, vec &result){
  // Chambers/Mallow/Stuck algorithm for simulating from stable(alpha, beta=1, mu=0, c=1)
  // Parameterization with beta=1, mu=0, alpha < 1 restricts support to [0, Inf)
  double u, w;
  const double zeta = -tan(PI*alpha/2);
  const double xi   = PI/2;
  const double z1   = pow(1+zeta*zeta, 1/(2*alpha));
  double z2, z3, z4;
  
  for(int i=0; i < nknots; i++){
    u = R::runif(-PI/2, PI/2);
    w = R::rexp(1.0);
    
    z2 = sin(alpha*(u+xi));
    z3 = pow(cos(u), 1/alpha);
    z4 = cos(u-alpha*(u+xi))/w;
    
    //result.at(i) = z1 * (z2 / z3) * pow(z4, (1-alpha)/alpha);
    result.at(i) = (z2 / z3) * pow(z4, (1-alpha)/alpha);
  }
  
  return;
}

void calcDistance(mat &x, mat &dist, int n){
  for(int i = 0; i < n; i++){
    for(int j = i; j < n; j++){
      dist(i,j) = norm(x.row(i) - x.row(j));
    }
  }
  dist = symmatu(dist);
  
  return;
}

void covMatern(int n, mat &distMat, double sigma2, double phi, double nu, mat &Sigma){
  double temp = sigma2 / (pow(2, nu-1) * R::gammafn(nu));
  
  for(int i = 0; i < n; i++){
    Sigma.at(i,i) = sigma2;
    for(int j = i+1; j < n; j++){
      Sigma.at(i,j) = temp * pow(phi*distMat.at(i,j), nu) * R::bessel_k(phi*distMat.at(i,j), nu, 1);
    }
  }
  
  Sigma = symmatu(Sigma);
  
  return;
}

void covPowExp(int n, mat &distMat, double sigma2, double phi, double gamma, mat &Sigma){
  for(int i = 0; i < n; i++){
    Sigma.at(i,i) = sigma2;
    for(int j = i+1; j < n; j++){
      Sigma.at(i,j) = sigma2 * exp(-1 * pow(phi*distMat.at(i,j), gamma));
    }
  }
  
  Sigma = symmatu(Sigma);
  
  return;
}

void mspBrownResnick(bool keepPsi, int n, int &count, int &curr_storage, double cstar, 
                     vec &Z, const vec &zero, const mat &Sigma,
                     mat &psiStore, vec &zetaStore){
  double minZ = 0;
  double zetaInv = R::rexp(1);
  vec psi = vec(n);
  vec tempVec = vec(n);
  
  while(cstar/zetaInv > minZ){
    rmvnorm(zero, Sigma, tempVec, n);
    psi = exp(tempVec - Sigma.diag(0)/2.0);
    Z = max(Z, psi / zetaInv);
    
    zetaInv += R::rexp(1);
    minZ = min(Z);
    
    if(keepPsi){
      storeCalc(psiStore, zetaStore, psi, zetaInv, n, count, curr_storage);
    }
  }
  return; 
}

void mspExtremeT(bool keepPsi, int n, int &count, int &curr_storage, double cstar, 
                 vec &Z, const vec &zero, const mat &Sigma, const double alpha, 
                 mat &psiStore, vec &zetaStore){
  double minZ = 0;
  double zetaInv = R::rexp(1);
  double norm_constant = sqrt(PI) * pow(2, -(alpha-2)/2) / tgamma((alpha+1)/2);
  
  vec psi = vec(n);
  vec tempVec = vec(n);
  
  while(cstar/zetaInv > minZ){
    rmvnorm(zero, Sigma, tempVec, n);
    psi = pow(max(tempVec, zero), alpha);
    Z = max(Z, norm_constant * psi / zetaInv);
    
    zetaInv += R::rexp(1);
    minZ = min(Z);
    
    if(keepPsi){
      storeCalc(psiStore, zetaStore, psi, zetaInv, n, count, curr_storage);
    }
  }
  return; 
}

void storeCalc(mat &psiStore, vec &zetaStore, vec &psi, double &zetaInv, int &n, int &count, int &curr_storage){
  psiStore.row(count) = psi.t();
  zetaStore.at(count) = 1/zetaInv;
  count++;
  
  if(count == curr_storage){
    psiStore.resize(2*curr_storage, n);
    zetaStore.resize(2*curr_storage);
    curr_storage = 2*curr_storage;
  }
}

void storeCalcSmith(mat &psiStore, mat &SStore, vec &zetaStore, vec &psi, vec &S, double &zetaInv, int &n, int &d, int &count, int &curr_storage){
  psiStore.row(count) = psi.t();
  SStore.row(count) = S.t();
  zetaStore.at(count) = 1/zetaInv;
  count++;
  
  if(count == curr_storage){
    psiStore.resize(2*curr_storage, n);
    SStore.resize(2*curr_storage, d);
    zetaStore.resize(2*curr_storage);
    curr_storage = 2*curr_storage;
  }
}