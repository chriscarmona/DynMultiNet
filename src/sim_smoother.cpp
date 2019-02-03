#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "DynMultiNet_shared.h"

// [[Rcpp::interfaces(r, cpp)]]

//' @description
//'    Performs Kalman filtering and smoothing or simulation smoothing.
//' 
//' @param y matrix. Matrix with observations. dim(y)=c(nobs,T_steps).
//' 
//' @param dd matrix. Trend parameter of the observation equation. dim(dd)=c(nobs,T_steps).
//' @param ZZ cube (3D array). Mapping matrix of the states into observations. dim(ZZ)=c(nobs,nstates,T_steps).
//' @param HHchol cube. Cholesky of the covariance of the shocks for the observations. dim(HHchol)=c(nobs,nobsshocks,T_steps).
//' 
//' @param cc matrix. Trend parameter of the state equation. dim(cc)=c(nstates,T_steps).
//' @param TT cube. Mapping matrix of the states. dim(TT)=c(nstates,nstates,T_steps).
//' @param RR cube. Cholesky of the covariance of the shocks for the states. dim(RR)=c(nstates,nstatesshocks,T_steps).
//' 
//' @param a1 matrix. Mean of the initial state. dim(a1)=c(nstates,1).
//' @param P1chol matrix. Cholesky of the covariance of the initial states. dim(P1chol)=c(nstates,nstates).
//' 
//' @param ind_output integer. Indicates the output that for the function, see details.
//' @param verbose boolean. Verbose execution.
//' 
//' @details
//'    The model assumes a latent variable approach.
//'    
//'    y(t) = dd(t) + ZZ(t) * a(t) + eps(t),  eps(t) ~ N(0,HH(t))
//'    a(t+1) = cc(t) + TT(t) * a(t) + RR(t) * eta(t), eta(t) ~ N(0,QQ(t))
//'    a(1) ~ N(a1,P1)
//'    
//'    The simulation smoother is a version of Algorithm 2 of Durbin and Koopman (2002)
//'    *WITH* the efficient modification that saves one pass of the filter/smoother.
//'    
//'    The initialization of the simulation smoother follows:
//'    Jarocinski (2015) A note on implementing the Durbin and Koopman simulation smoother,
//'    
//'    - Missing data not allowed (y does not contain any NaN)
//' 
//' @return
//'    A list with the following components:
//' \describe{
//'     \item{\code{aaa}}{if \code{ind_output=1}, the mean of the states conditional on y, nstates x T.}
//'     \item{\code{loglik}}{log likelihood of each observation of y, 1 x T.}
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List kfsim( const arma::mat y,
                  
                  const arma::mat dd,
                  const arma::cube ZZ,
                  const arma::cube HHchol,
                  
                  const arma::mat cc,
                  const arma::cube TT,
                  const arma::cube RR,
                  const arma::cube QQchol,
                  
                  const arma::colvec a1,
                  const arma::mat P1chol,
                  
                  const unsigned int ind_output,
                  const bool verbose=false ) {
  // PURPOSE: This function performs Kalman filtering and smoothing or simulation smoothing.
  // The simulation smoother is a version of Algorithm 2 of Durbin and Koopman (2002)
  // *WITH* the efficient modification that saves one pass of the filter/smoother.
  // The initialization of the simulation smoother follows Jarocinski (2015) A
  // note on implementing the Durbin and Koopman simulation smoother,
  // 
  // - Missing data not allowed (y does not contain any NaN)
  // We assume the model is:
  // y(t) = dd(t) + ZZ a(t) + HHchol * eps(t),  eps(t) ~ N(0,I)
  // a(t+1) = cc(t) + TT a(t) + RR * QQchol * eta(t), eta(t) ~ N(0,I)
  // a(1) ~ a1 + P1chol * N(0,I)
  // 
  // INPUTS:
  // y - data, nobs x T (where nobs is the number of observables and T is the
  //     number of time periods)
  // dd, ZZ, HH - parameters of the observation equation
  // cc, TT, RR - parameters of the state equation
  // a1, P1chol - parameters of the distribution of the initial state
  // ind_output - 0: return the log likelihood [NaN,loglik]
  //              1: return the smoothed state and the log likelihood [aaa,loglik]
  //              2: return a draw of the states from the simulation smoother [aaa,NaN]
  // 
  // OUTPUTS:
  // if ind_output = 0
  //   aaa = []
  //   loglik = log likelihood of each observation of y, 1 x T
  // if ind_output = 1
  //   aaa = the mean of the states conditional on y, nstates x T
  //   loglik = log likelihood of each observation of y, 1 x T
  // if ind_output = 2
  //   aaa = a draw of the states conditional on y, nstates x T
  //   loglik = []
  //
  // Log:
  // 2013-10-06 created in Matlab by Marek Jarocinski
  // 2014-01-16 check pos-def of Ft, added the option of reporting only loglik
  // 2019-01-30 modified by Chris Carmona: translated from Matlab to C++, added time-varying ZZ, HH, TT, RR, QQ
  
  if(verbose){ Rcpp::Rcout << "Simulation smoother for state space models" << std::endl;}
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  arma::cube aux_cube_1;
  arma::mat aux_mat_1;
  arma::vec aux_vec_1;
  
  // measure dimensions
  unsigned int nobs= y.n_rows;
  unsigned int T_steps = y.n_cols;
  unsigned int nobsshocks= HHchol.n_rows;
  unsigned int nstates= RR.n_rows;
  unsigned int nstatesshocks= QQchol.n_rows;
  
  // Output objects
  arma::mat aaa = arma::zeros<arma::mat>(nstates,T_steps);
  arma::mat loglik = arma::zeros<arma::mat>(1,T_steps);
  
  // Checking dimensions
  if(nobs!=nobsshocks){ throw std::range_error("You are using a different number of observations and shocks, maybe a model that multiplies S(t)*eps(t)?");}
  if((dd.n_rows!=nobs)|(dd.n_cols!=T_steps)){ throw std::range_error("Input size mismatch, check dd");}
  if((ZZ.n_rows!=nobs)|(ZZ.n_cols!=nstates)|(ZZ.n_slices!=T_steps)){ throw std::range_error("Input size mismatch, check ZZ");}
  if((HHchol.n_rows!=nobsshocks)|(HHchol.n_cols!=nobsshocks)|(HHchol.n_slices!=T_steps)){ throw std::range_error("Input size mismatch, check HHchol");}
  if((cc.n_rows!=nstates)|(cc.n_cols!=T_steps)){ throw std::range_error("Input size mismatch, check cc");}
  if((TT.n_rows!=nstates)|(TT.n_cols!=nstates)|(TT.n_slices!=T_steps)){ throw std::range_error("Input size mismatch, check TT");}
  if((RR.n_rows!=nstates)|(RR.n_cols!=nstatesshocks)|(RR.n_slices!=T_steps)){ throw std::range_error("Input size mismatch, check RR");}
  if((QQchol.n_rows!=nstatesshocks)|(QQchol.n_cols!=nstatesshocks)|(QQchol.n_slices!=T_steps)){ throw std::range_error("Input size mismatch, check QQchol");}
  if(a1.n_rows!=nstates){ throw std::range_error("Input size mismatch, check a1");}
  if((P1chol.n_rows!=nstates)|(P1chol.n_cols!=nstates)){ throw std::range_error("Input size mismatch, check P1chol");}
  
  arma::mat yplus = arma::zeros<arma::mat>(nobs,T_steps);
  arma::mat aplus = arma::zeros<arma::mat>(nstates,T_steps+1);
  arma::mat yy = arma::zeros<arma::mat>(nobs,T_steps);
  
  if(ind_output==3){
    // Simulate y and alpha from their unconditional distribution
    aplus.col(0) = a1+P1chol*arma::randn(nstates,1); // draw the first state with a1
    if(verbose){ Rcpp::Rcout << "Generating y, alpha" << std::endl << " t = ";}
    for( t=0; t<T_steps; t++ ) {
      if(verbose){ Rcpp::Rcout << t+1 << ",";}
      yplus.col(t) = dd.col(t) + ZZ.slice(t)*aplus.col(t) + HHchol.slice(t)*arma::randn(nobsshocks,1);
      aplus.col(t+1) = cc.col(t) + TT.slice(t)*aplus.col(t) + RR.slice(t)*(QQchol.slice(t)*arma::randn(nstatesshocks,1));
    }
    aplus.shed_col(aplus.n_cols-1);
    
    return Rcpp::List::create( Rcpp::Named("y") = yplus,
                               Rcpp::Named("alpha") = aplus  );
  }
  
  if(ind_output==2){
    // Durbin, Koopman 2002, Algorithm 2.
    // Generate yplus and aplus - y and a drawn from their unconditional
    // distribution *using the 'demeaned' model*,
    // i.e. zero initial state and zero constant terms!
    // For explanation see Jarocinski (2015), A note on implementing 
    // the Durbin and Koopman simulation smoother.
    
    aplus.col(0) = P1chol*arma::randn(nstates,1); // draw the first state with a1=0
    
    if(verbose){ Rcpp::Rcout << "Generating yplus, aplus" << std::endl << " t = ";}
    for( t=0; t<T_steps; t++ ) {
      if(verbose){ Rcpp::Rcout << t+1 << ",";}
      yplus.col(t) = ZZ.slice(t)*aplus.col(t) + HHchol.slice(t)*arma::randn(nobsshocks,1);
      aplus.col(t+1) = TT.slice(t)*aplus.col(t) + RR.slice(t)*(QQchol.slice(t)*arma::randn(nstatesshocks,1));
    }
    aplus.shed_col(aplus.n_cols-1);
    yy = y - yplus;
    if(verbose){ Rcpp::Rcout << std::endl << "done." << std::endl;}
  } else {
    yy = y;
  }
  
  // allocate space
  arma::mat vvv = arma::zeros<arma::mat>(nobs,T_steps);
  arma::cube FFFinv = arma::zeros<arma::cube>(nobs,nobs,T_steps);
  arma::cube KKK = arma::zeros<arma::cube>(nstates,nobs,T_steps);
  arma::cube LLL = arma::zeros<arma::cube>(nstates,nstates,T_steps);
  
  // Kalman filter on yy
  // compute frequently used matrix RQR'
  arma::cube RQR = arma::zeros<arma::cube>(nstates,nstates,T_steps);
  for( t=0; t<T_steps; t++ ) {
    RQR.slice(t) = RR.slice(t)*QQchol.slice(t)*trans(QQchol.slice(t))*trans(RR.slice(t));
  }
  // initialize Kalman filter
  arma::mat at = a1; // at|I(t-1)
  arma::mat Pt = P1chol*trans(P1chol); // Pt|I(t-1)
  arma::mat vt;
  arma::mat Ftinv;
  arma::mat Kt;
  arma::mat Lt;
  if(verbose){ Rcpp::Rcout << "Kalman filter, on yy" << std::endl << " t = ";}
  for( t=0; t<T_steps; t++ ) {
    if(verbose){ Rcpp::Rcout << t+1 << ",";}
    vt = yy.col(t) - dd.col(t) - ZZ.slice(t) * at;
    Ftinv = arma::inv_sympd( ZZ.slice(t)*Pt*trans(ZZ.slice(t)) + HHchol.slice(t)*trans(HHchol.slice(t)) );
    Kt = TT.slice(t)*Pt*trans(ZZ.slice(t))*Ftinv;
    Lt = TT.slice(t) - Kt * ZZ.slice(t);
    // update at,Pt; from now on their interpretation is "a(t+1),P(t+1)"
    at = cc.col(t) + TT.slice(t)*at + Kt*vt;
    Pt = TT.slice(t)*Pt*trans(Lt) + RQR.slice(t);
    // store the quantities needed later for smoothing
    vvv.col(t) = vt;
    FFFinv.slice(t) = Ftinv;
    KKK.slice(t) = Kt;
    LLL.slice(t) = Lt;
    if (ind_output<2){
      // store the log likelihood
      arma::mat cholFtinv = chol(Ftinv);
      loglik.col(t) = -0.5*nobs*log(2*arma::datum::pi) + sum(log(cholFtinv.diag())) - 0.5*(trans(vt)*Ftinv*vt);
      // Note that: 0.5*log(det(Ftinv)) = sum(log(diag(chol(Ftinv)))
    }
  }
  if(verbose){ Rcpp::Rcout << std::endl << "done." << std::endl;}
  
  if(ind_output==0) {
    aaa.fill(arma::datum::nan);
    return Rcpp::List::create( Rcpp::Named("aaa") = aaa,
                               Rcpp::Named("loglik") = loglik );
  }
  
  // Kalman smoother
  // backwards recursion on r
  arma::mat rrr = arma::zeros<arma::mat>(nstates,T_steps);
  if(verbose){ Rcpp::Rcout << "Kalman smoother, backwards recursion on r" << std::endl << " t = ";}
  for( t=T_steps-1; t>0; t-- ) {
    if(verbose){ Rcpp::Rcout << t+1 << ",";}
    // rrr.col(t-1) = trans(ZZ.slice(t))*FFFinv.slice(t)*vvv.col(t) + trans(TT.slice(t) - KKK.slice(t)*ZZ.slice(t))*rrr.col(t);
    rrr.col(t-1) = trans(ZZ.slice(t))*FFFinv.slice(t)*vvv.col(t) + trans(LLL.slice(t))*rrr.col(t);
  }
  // one more iteration to get r0
  arma::mat r0 = trans(ZZ.slice(0))*FFFinv.slice(0)*vvv.col(0) + trans(TT.slice(0) - KKK.slice(0)*ZZ.slice(0))*rrr.col(0);
  if(verbose){ Rcpp::Rcout << std::endl << "done." << std::endl;}
  
  // forwards recursion to compute smoothed states from r - DK(2002),eq.8
  aaa = arma::zeros<arma::mat>(nstates,T_steps);
  aaa.col(0) = a1 + P1chol*trans(P1chol)*r0; // initialize the forward recursion
  if(verbose){ Rcpp::Rcout << "Kalman smoother, computing smoothed states" << std::endl << " t = ";}
  for( t=1; t<T_steps; t++ ) {
    if(verbose){ Rcpp::Rcout << t+1 << ",";}
    aaa.col(t) = cc.col(t) + TT.slice(t)*aaa.col(t-1) + RQR.slice(t)*(rrr.col(t-1));
  }
  if(verbose){ Rcpp::Rcout << std::endl << "done." << std::endl;}
  
  if (ind_output==2){
    aaa = aaa + aplus;
    loglik.fill(arma::datum::nan);
  }
  
  if(verbose){ Rcpp::Rcout << "Finished!!!" << std::endl;}
  
  return Rcpp::List::create( Rcpp::Named("aaa") = aaa,
                             Rcpp::Named("loglik") = loglik );
}

// [[Rcpp::export]]
arma::mat kfsim_cpp( const arma::mat y,
                     
                     const arma::mat dd,
                     const arma::cube ZZ,
                     const arma::cube HHchol,
                     
                     const arma::mat cc,
                     const arma::cube TT,
                     const arma::cube RR,
                     const arma::cube QQchol,
                     
                     const arma::colvec a1,
                     const arma::mat P1chol ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  arma::cube aux_cube_1;
  arma::mat aux_mat_1;
  arma::vec aux_vec_1;
  
  // measure dimensions
  unsigned int nobs= y.n_rows;
  unsigned int T_steps = y.n_cols;
  unsigned int nobsshocks= HHchol.n_rows;
  unsigned int nstates= RR.n_rows;
  unsigned int nstatesshocks= QQchol.n_rows;
  
  // Output objects
  arma::mat aaa = arma::zeros<arma::mat>(nstates,T_steps);
  
  // Checking dimensions
  if(nobs!=nobsshocks){ throw std::range_error("You are using a different number of observations and shocks, maybe a model that multiplies S(t)*eps(t)?");}
  if((dd.n_rows!=nobs)|(dd.n_cols!=T_steps)){ throw std::range_error("Input size mismatch, check dd");}
  if((ZZ.n_rows!=nobs)|(ZZ.n_cols!=nstates)|(ZZ.n_slices!=T_steps)){ throw std::range_error("Input size mismatch, check ZZ");}
  if((HHchol.n_rows!=nobsshocks)|(HHchol.n_cols!=nobsshocks)|(HHchol.n_slices!=T_steps)){ throw std::range_error("Input size mismatch, check HHchol");}
  if((cc.n_rows!=nstates)|(cc.n_cols!=T_steps)){ throw std::range_error("Input size mismatch, check cc");}
  if((TT.n_rows!=nstates)|(TT.n_cols!=nstates)|(TT.n_slices!=T_steps)){ throw std::range_error("Input size mismatch, check TT");}
  if((RR.n_rows!=nstates)|(RR.n_cols!=nstatesshocks)|(RR.n_slices!=T_steps)){ throw std::range_error("Input size mismatch, check RR");}
  if((QQchol.n_rows!=nstatesshocks)|(QQchol.n_cols!=nstatesshocks)|(QQchol.n_slices!=T_steps)){ throw std::range_error("Input size mismatch, check QQchol");}
  if(a1.n_rows!=nstates){ throw std::range_error("Input size mismatch, check a1");}
  if((P1chol.n_rows!=nstates)|(P1chol.n_cols!=nstates)){ throw std::range_error("Input size mismatch, check P1chol");}
  
  arma::mat yplus = arma::zeros<arma::mat>(nobs,T_steps);
  arma::mat aplus = arma::zeros<arma::mat>(nstates,T_steps+1);
  arma::mat yy = arma::zeros<arma::mat>(nobs,T_steps);
  
  // Durbin, Koopman 2002, Algorithm 2.
  // Generate yplus and aplus - y and a drawn from their unconditional
  // distribution *using the 'demeaned' model*,
  // i.e. zero initial state and zero constant terms!
  // For explanation see Jarocinski (2015), A note on implementing 
  // the Durbin and Koopman simulation smoother.
  
  aplus.col(0) = P1chol*arma::randn(nstates,1); // draw the first state with a1=0
  
  for( t=0; t<T_steps; t++ ) {
    yplus.col(t) = ZZ.slice(t)*aplus.col(t) + HHchol.slice(t)*arma::randn(nobsshocks,1);
    aplus.col(t+1) = TT.slice(t)*aplus.col(t) + RR.slice(t)*(QQchol.slice(t)*arma::randn(nstatesshocks,1));
  }
  aplus.shed_col(aplus.n_cols-1);
  yy = y - yplus;
  
  // allocate space
  arma::mat vvv = arma::zeros<arma::mat>(nobs,T_steps);
  arma::cube FFFinv = arma::zeros<arma::cube>(nobs,nobs,T_steps);
  arma::cube KKK = arma::zeros<arma::cube>(nstates,nobs,T_steps);
  arma::cube LLL = arma::zeros<arma::cube>(nstates,nstates,T_steps);
  
  // Kalman filter on yy
  // compute frequently used matrix RQR'
  arma::cube RQR = arma::zeros<arma::cube>(nstates,nstates,T_steps);
  for( t=0; t<T_steps; t++ ) {
    RQR.slice(t) = RR.slice(t)*QQchol.slice(t)*trans(QQchol.slice(t))*trans(RR.slice(t));
  }
  // initialize Kalman filter
  arma::mat at = a1; // at|I(t-1)
  arma::mat Pt = P1chol*trans(P1chol); // Pt|I(t-1)
  arma::mat vt;
  arma::mat Ftinv;
  arma::mat Kt;
  arma::mat Lt;
  for( t=0; t<T_steps; t++ ) {
    vt = yy.col(t) - dd.col(t) - ZZ.slice(t) * at;
    Ftinv = arma::inv_sympd( ZZ.slice(t)*Pt*trans(ZZ.slice(t)) + HHchol.slice(t)*trans(HHchol.slice(t)) );
    Kt = TT.slice(t)*Pt*trans(ZZ.slice(t))*Ftinv;
    Lt = TT.slice(t) - Kt * ZZ.slice(t);
    // update at,Pt; from now on their interpretation is "a(t+1),P(t+1)"
    at = cc.col(t) + TT.slice(t)*at + Kt*vt;
    Pt = TT.slice(t)*Pt*trans(Lt) + RQR.slice(t);
    // store the quantities needed later for smoothing
    vvv.col(t) = vt;
    FFFinv.slice(t) = Ftinv;
    KKK.slice(t) = Kt;
    LLL.slice(t) = Lt;
  }
  
  // Kalman smoother
  // backwards recursion on r
  arma::mat rrr = arma::zeros<arma::mat>(nstates,T_steps);
  for( t=T_steps-1; t>0; t-- ) {
    // rrr.col(t-1) = trans(ZZ.slice(t))*FFFinv.slice(t)*vvv.col(t) + trans(TT.slice(t) - KKK.slice(t)*ZZ.slice(t))*rrr.col(t);
    rrr.col(t-1) = trans(ZZ.slice(t))*FFFinv.slice(t)*vvv.col(t) + trans(LLL.slice(t))*rrr.col(t);
  }
  // one more iteration to get r0
  arma::mat r0 = trans(ZZ.slice(0))*FFFinv.slice(0)*vvv.col(0) + trans(TT.slice(0) - KKK.slice(0)*ZZ.slice(0))*rrr.col(0);
  
  // forwards recursion to compute smoothed states from r - DK(2002),eq.8
  aaa = arma::zeros<arma::mat>(nstates,T_steps);
  aaa.col(0) = a1 + P1chol*trans(P1chol)*r0; // initialize the forward recursion
  for( t=1; t<T_steps; t++ ) {
    aaa.col(t) = cc.col(t) + TT.slice(t)*aaa.col(t-1) + RQR.slice(t)*(rrr.col(t-1));
  }
  
  aaa = aaa + aplus;
  
  return aaa;
}
