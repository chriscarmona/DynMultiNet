#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "DynMultiNet_shared.h"

// [[Rcpp::interfaces(r, cpp)]]



// [[Rcpp::export]]
Rcpp::List sample_baseline_t_link_nGP_cpp( arma::colvec eta_t,
                                           const arma::cube y_ijt,
                                           const arma::cube w_ijt,
                                           arma::cube gamma_ijt,
                                           const arma::cube nGP_G_t,
                                           const arma::cube nGP_H_t,
                                           const arma::cube nGP_Wchol_t,
                                           const bool directed=false ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  arma::uvec aux_uvec_1;
  arma::colvec aux_vec;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::cube aux_cube_1;
  
  // Rcpp::Rcout << "(0)" << std::endl;
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  
  // initialising objects
  arma::mat alpha_t = arma::zeros<arma::mat>(3,T_net);
  alpha_t.row(0) = eta_t.t(); // state space have T_net columns
  arma::mat Y_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat PG_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat linpred_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat C_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat Z_t = arma::zeros<arma::mat>(1,T_net);
  arma::cube X_t = arma::ones<arma::cube>(1,1,T_net); 
  
  if( directed ){
    // Model matrix
    X_t = arma::zeros<arma::cube>(V_net*(V_net-1),3,T_net); 
    X_t.col(0).ones();
    
    // Network data
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = y_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
      aux_mat_2.shed_row(i);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    Y_t = aux_mat_1;
    
    // Polya-Gamma data
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=0; i<V_net; i++ ) {
      // Rcpp::Rcout << "i=" << i << std::endl;
      aux_mat_2 = w_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
      aux_mat_2.shed_row(i);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    PG_t = aux_mat_1;
    
    // linear predictor
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
      aux_mat_2.shed_row(i);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    linpred_t = aux_mat_1;
  } else {
    // Model matrix
    // Has three columns for: U, U', A, only the first one goes is 1 to have effect on Y
    X_t = arma::zeros<arma::cube>(V_net*(V_net-1)/2,3,T_net); 
    X_t.col(0).ones();
    
    // Network data
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = y_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    Y_t = aux_mat_1;
    
    // Polya-Gamma data
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = w_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    PG_t = aux_mat_1;
    
    // linear predictor
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = gamma_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    linpred_t = aux_mat_1;
  }
  // Constant term for theta in the linear predictor
  C_t = linpred_t; C_t.zeros();
  
  for( t=0; t<T_net; t++ ) {
    // Rcpp::Rcout << "t=" << t << std::endl;
    C_t.col(t) = linpred_t.col(t) - ( X_t.slice(t) * alpha_t.col(t) );
  }
  
  // Adjusted response
  Z_t = (Y_t-0.5)/PG_t - C_t;
  arma::cube Z_t_cov_chol = arma::zeros<arma::cube>(Z_t.n_rows,Z_t.n_rows,T_net);
  for( t=0; t<T_net; t++ ) {
    Z_t_cov_chol.slice(t).diag() = pow( PG_t.col(t) , -1);
  }
  
  // Sampling coefficient using State-Space model and simulation smoother
  arma::mat dd = Z_t; dd.zeros();
  arma::mat cc = arma::zeros<arma::mat>(3,T_net);
  arma::colvec a1 = arma::zeros<arma::colvec>(3);
  arma::mat P1chol = arma::eye<arma::mat>(3,3); P1chol.diag().fill(1/100);
  
  alpha_t = kfsim_cpp( Z_t,
                       
                       dd,
                       X_t,
                       Z_t_cov_chol,
                       
                       cc,
                       nGP_G_t,
                       nGP_H_t,
                       nGP_Wchol_t,
                       
                       a1,
                       P1chol );
  
  // Recalculate linpred with the new values of eta_t
  for( t=0; t<T_net; t++ ) {
    linpred_t.col(t) = C_t.col(t) + ( X_t.slice(t) * alpha_t.col(t) );
  }
  
  // updating value of eta
  eta_t = trans(alpha_t.row(0));
  
  // Update gamma_ijt with the new values of linpred
  if( directed ){
    aux_mat_1 = linpred_t;
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = aux_mat_1.rows(i*(V_net-1),(i+1)*(V_net-1)-1);
      aux_mat_2.insert_rows(i,arma::zeros<arma::mat>(1,T_net));
      gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_2;
    }
  } else {
    aux_mat_1 = linpred_t;
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = aux_mat_1.rows((i-1)*V_net-((i-1)*i)/2,i*V_net-(i*(i+1))/2-1);
      gamma_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1) = aux_mat_2;
    }
  }
  
  return Rcpp::List::create( Rcpp::Named("eta_t") = eta_t,
                             Rcpp::Named("gamma_ijt") = gamma_ijt,
                             Rcpp::Named("alpha_eta_t") = alpha_t );
}
