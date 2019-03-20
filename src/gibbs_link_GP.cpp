#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "DynMultiNet_shared.h"

// [[Rcpp::interfaces(r, cpp)]]



// [[Rcpp::export]]
Rcpp::List sample_baseline_t_link_GP_cpp( arma::colvec eta_t,
                                          
                                          const arma::cube y_ijt,
                                          const arma::cube w_ijt,
                                          arma::cube gamma_ijt,
                                          
                                          const arma::mat eta_t_cov_prior_inv,
                                          
                                          const bool lat_mean,
                                          double eta_t_bar,
                                          const double sigma_eta_bar,
                                          
                                          const bool directed=false ) {
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  arma::uvec aux_uvec_1;
  arma::colvec aux_vec;
  arma::colvec aux_colvec;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::cube aux_cube_1;
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  
  arma::colvec Y = arma::zeros<arma::colvec>(1);
  arma::colvec W = arma::zeros<arma::colvec>(1);
  arma::colvec linpred = arma::zeros<arma::colvec>(1);
  
  arma::colvec C = arma::zeros<arma::colvec>(1);
  arma::colvec C_all = C;
  arma::colvec Z = arma::zeros<arma::colvec>(1);
  
  aux_mat_1 = arma::zeros<arma::mat>(T_net,T_net); aux_mat_1.eye();
  if( directed ){
    aux_mat_2 = arma::ones<arma::mat>(V_net*(V_net-1),1);
  } else {
    aux_mat_2 = arma::ones<arma::mat>(V_net*(V_net-1)/2,1);
  }
  arma::mat X = arma::kron( aux_mat_1, aux_mat_2 );
  arma::sp_mat X_sp = arma::sp_mat(X);
  arma::sp_mat X_sp_valid = X_sp;
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  
  // Marginal posterior covariance matrix for eta_t
  arma::mat eta_t_cov_inv;
  arma::mat eta_t_cov;
  
  // GP prior mean at zero
  if(!lat_mean){
    eta_t_bar=0;
  }
  // prior mean vector for theta_t
  arma::colvec eta_t_mean_prior = arma::zeros<arma::colvec>(T_net);
  // Marginal posterior mean vector for theta_t
  arma::colvec eta_t_mean = arma::zeros<arma::colvec>(T_net);
  
  double eta_t_bar_mean;
  double eta_t_bar_var;
  
  if( directed ){
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = y_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
      aux_mat_2.shed_row(i);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
    Y = aux_mat_1;
    
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = w_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
      aux_mat_2.shed_row(i);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
    W = aux_mat_1;
    
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
      aux_mat_2.shed_row(i);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
    linpred = aux_mat_1;
  } else {
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = y_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
    Y = aux_mat_1;
    
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = w_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
    W = aux_mat_1;
    
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = gamma_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
    linpred = aux_mat_1;
  }
  
  // Constant term for theta in the linear predictor
  C = linpred - (X_sp * eta_t);
  // Adjusted Z for the covariance
  Z = (Y-0.5)/W - C;
  
  
  // identifies valid obs
  valid_obs = find_finite(Y); // find valid elements
  
  // Keep only valid cases
  C_all = C;
  Y = Y.rows(valid_obs);
  W = W.rows(valid_obs);
  linpred = linpred.rows(valid_obs);
  C = C.rows(valid_obs);
  Z = Z.rows(valid_obs);
  X_sp_valid = arma::sp_mat(X.rows(valid_obs));
  
  // marginal posterior covariance
  // Way 1:
  // Omega_sp.diag() = W;
  // eta_t_cov_inv = X_sp_valid.t() * Omega_sp * X_sp_valid;
  // Way 2 (more efficient)
  // eta_t_cov_inv is a diagonal matrix with the element (t,t) equal to thesum w_ijtk for all i,j
  eta_t_cov_inv = arma::zeros<arma::mat>(T_net,T_net);
  if( directed ){
    aux_vec = arma::zeros<arma::vec>(T_net*V_net*(V_net-1));
    for( t=0; t<T_net; t++ ) { aux_vec.rows( t*V_net*(V_net-1) , (t+1)*V_net*(V_net-1)-1 ).fill(t); }
  } else {
    aux_vec = arma::zeros<arma::vec>(T_net*V_net*(V_net-1)/2);
    for( t=0; t<T_net; t++ ) { aux_vec( arma::span(t*V_net*(V_net-1)/2,(t+1)*V_net*(V_net-1)/2-1) ).fill(t); }
  }
  aux_vec = aux_vec.rows(valid_obs);
  for( t=0; t<T_net; t++ ) { eta_t_cov_inv(t,t) = sum( W.elem( find(aux_vec==t)  ) ); }
  eta_t_cov_inv += eta_t_cov_prior_inv;
  eta_t_cov = arma::inv_sympd(eta_t_cov_inv);
  
  // marginal posterior mean
  aux_colvec = arma::zeros<arma::colvec>(T_net);
  eta_t_mean_prior = eta_t_bar * aux_colvec;
  eta_t_mean = eta_t_cov*(X_sp_valid.t() * (W % Z) + eta_t_cov_prior_inv * eta_t_mean_prior);
  
  // Sampling eta_t
  eta_t = arma::mvnrnd( eta_t_mean , eta_t_cov );
  
  // return eta_t;
  
  // Recalculate linpred with the new values of mu
  linpred = X_sp * eta_t + C_all;
  // Redefine gamma_ijt with the new values of linpred
  if( directed ){
    aux_mat_1 = linpred;
    aux_mat_1.reshape(V_net*(V_net-1),T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = aux_mat_1.rows(i*(V_net-1),(i+1)*(V_net-1)-1);
      aux_mat_2.insert_rows(i,arma::zeros<arma::mat>(1,T_net));
      gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_2;
    }
  } else {
    aux_mat_1 = linpred;
    aux_mat_1.reshape(V_net*(V_net-1)/2,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = aux_mat_1.rows((i-1)*V_net-((i-1)*i)/2,i*V_net-(i*(i+1))/2-1);
      gamma_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1) = aux_mat_2;
    }
  }
  
  // Sample GP prior mean
  if(lat_mean){
    aux_colvec = arma::randn(1);
    eta_t_bar = aux_colvec(0);
    
    eta_t_bar_var = 1/( accu(eta_t_cov_prior_inv)+pow(sigma_eta_bar,-2) );
    
    aux_mat_1 = arma::ones<arma::mat>(1,T_net);
    aux_colvec = eta_t_bar_var * (aux_mat_1*eta_t_cov_prior_inv*eta_t);
    eta_t_bar_mean = aux_colvec(0);
    
    eta_t_bar = eta_t_bar_mean + eta_t_bar * sqrt(eta_t_bar_var);
  }
  
  return Rcpp::List::create( Rcpp::Named("eta_t") = eta_t,
                             Rcpp::Named("gamma_ijt") = gamma_ijt,
                             Rcpp::Named("eta_t_bar") = eta_t_bar );
}


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_link_GP_cpp( arma::cube ab_ith,
                                         const arma::cube y_ijt,
                                         const arma::cube w_ijt,
                                         arma::cube gamma_ijt,
                                         const arma::mat ab_t_sigma_prior_inv,
                                         const arma::colvec tau_h ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = ab_ith.n_slices;
  
  ab_ith = reshape(ab_ith,V_net,T_net*H_dim,1);
  arma::mat ab_ith_mat = ab_ith.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec W = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec linpred = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec C_all = C;
  arma::colvec Z = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::sp_mat Omega_sp=arma::speye<arma::sp_mat>((V_net-1)*T_net,(V_net-1)*T_net);
  
  arma::mat tau_h_diag(tau_h.n_rows,tau_h.n_rows); tau_h_diag.eye();
  tau_h_diag.diag() = tau_h;
  
  arma::mat x_i_cov_prior = kron( tau_h_diag, ab_t_sigma_prior_inv );
  arma::mat x_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat x_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  
  arma::colvec x_i_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  // Performing polya-gamma sampling in ab_ith...
  // The rows in ab_ith_mat will act as the targeted "coefficients"
  // The columns in X_all will act as the "covariates"
  // The order of the response gamma_ijth and the latent w_ijtk will depend on the order of rows in X_all
  
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  arma::mat X = X_all;
  arma::sp_mat X_sp = arma::sp_mat(X);
  arma::sp_mat X_sp_valid = X_sp;
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  
  for( i=0; i<V_net; i++ ) {
    aux_mat_1 = y_ijt.subcube(i,0,0, i,i,T_net-1);
    aux_mat_2 = y_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
    aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
    aux_mat_1.shed_rows(i,i+1);
    aux_mat_1 = aux_mat_1.t();
    aux_mat_1.reshape((V_net-1)*T_net,1);
    Y = aux_mat_1;
    
    aux_mat_1 = w_ijt.subcube(i,0,0, i,i,T_net-1);
    aux_mat_2 = w_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
    aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
    aux_mat_1.shed_rows(i,i+1);
    aux_mat_1 = aux_mat_1.t();
    aux_mat_1.reshape((V_net-1)*T_net,1);
    W = aux_mat_1;
    
    aux_mat_1 = gamma_ijt.subcube(i,0,0, i,i,T_net-1);
    aux_mat_2 = gamma_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
    aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
    aux_mat_1.shed_rows(i,i+1);
    aux_mat_1 = aux_mat_1.t();
    aux_mat_1.reshape((V_net-1)*T_net,1);
    linpred = aux_mat_1;
    
    // X
    X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
    for( t=0; t<T_net; t++ ) {
      X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = ab_ith_mat.submat( aux_uvec_3, aux_uvec_2+t );
    }
    X = X_all;
    X.shed_rows(T_net*i,T_net*(i+1)-1);
    X_sp = arma::sp_mat(X);
    
    // Constant term in linear predictor
    C = linpred - X_sp * trans(ab_ith_mat.row(i));
    
    // identifies valid obs
    valid_obs = find_finite(Y); // find valid elements
    
    // Keep only valid cases
    C_all = C;
    Y = Y.rows(valid_obs);
    W = W.rows(valid_obs);
    linpred = linpred.rows(valid_obs);
    C = C.rows(valid_obs);
    X_sp_valid = arma::sp_mat(X.rows(valid_obs));
    
    // Marginal posterior covariance
    Omega_sp = arma::speye<arma::sp_mat>(W.n_rows,W.n_rows);
    Omega_sp.diag() = W;
    x_i_cov_inv = X_sp_valid.t() * Omega_sp * X_sp_valid;
    x_i_cov_inv = x_i_cov_inv + x_i_cov_prior ;
    x_i_cov = arma::inv_sympd(x_i_cov_inv);
    
    // Marginal posterior mean
    Z = (Y-0.5)/W - C;
    x_i_mean = x_i_cov*(X_sp_valid.t() * (W % Z));
    
    // Sampling ab_ith_mat
    ab_ith_mat.row(i) = trans(arma::mvnrnd( x_i_mean , x_i_cov ));
    
    // Recalculate linpred with the new values of ab_ith_mat
    linpred = X_sp * trans(ab_ith_mat.row(i)) + C_all;
    // Redefine gamma_ijt with the new values of linpred
    aux_mat_1 = linpred;
    aux_mat_1.reshape(T_net,(V_net-1));
    aux_mat_1 = aux_mat_1.t();
    aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
    gamma_ijt.subcube(i,0,0, i,i,T_net-1) = aux_mat_1.rows(0,i);
    gamma_ijt.subcube(i,i,0, V_net-1,i,T_net-1) = aux_mat_1.rows(i,V_net-1);
  }
  
  // get ab_ith from ab_ith_mat
  ab_ith.slice(0)=ab_ith_mat;
  ab_ith = reshape(ab_ith,V_net,T_net,H_dim);
  
  return Rcpp::List::create( Rcpp::Named("ab_ith") = ab_ith,
                             Rcpp::Named("gamma_ijt") = gamma_ijt );
}


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_shared_link_GP_cpp( arma::cube ab_ith,
                                                const arma::field<arma::cube> y_ijtk,
                                                const arma::field<arma::cube> w_ijtk,
                                                arma::field<arma::cube> gamma_ijtk,
                                                const arma::mat ab_t_sigma_prior_inv,
                                                const arma::colvec tau_h ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int k=0;
  unsigned int t=0;
  
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  arma::cube y_ijt = y_ijtk(0);
  arma::cube w_ijt = w_ijtk(0);
  arma::cube gamma_ijt = gamma_ijtk(0);
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int K_net = y_ijtk.n_rows;
  unsigned int H_dim = ab_ith.n_slices;
  
  ab_ith = reshape(ab_ith,V_net,T_net*H_dim,1);
  arma::mat ab_ith_mat = ab_ith.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec W = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec linpred = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec C_all = C;
  arma::colvec Z = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  arma::sp_mat Omega_sp=arma::speye<arma::sp_mat>((V_net-1)*T_net*K_net,(V_net-1)*T_net*K_net);
  
  arma::mat tau_h_diag(tau_h.n_rows,tau_h.n_rows); tau_h_diag.eye();
  tau_h_diag.diag() = tau_h;
  
  arma::mat x_i_cov_prior = kron( tau_h_diag, ab_t_sigma_prior_inv );
  arma::mat x_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat x_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  
  arma::colvec x_i_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  arma::mat X = X_all;
  arma::sp_mat X_sp = arma::sp_mat(X);
  arma::sp_mat X_sp_valid = X_sp;
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  
  for( i=0; i<V_net; i++ ) {
    // Rcpp::Rcout << "i=" << i << std::endl;
    
    arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
    arma::colvec W = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
    arma::colvec linpred = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
    
    for( k=0; k<K_net; k++ ) {
      // Rcpp::Rcout << "k=" << k << std::endl;
      y_ijt = y_ijtk(k);
      aux_mat_1 = y_ijt.subcube(i,0,0, i,i,T_net-1);
      aux_mat_2 = y_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
      aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
      aux_mat_1.shed_rows(i,i+1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape((V_net-1)*T_net,1);
      Y.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
      
      w_ijt = w_ijtk(k);
      aux_mat_1 = w_ijt.subcube(i,0,0, i,i,T_net-1);
      aux_mat_2 = w_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
      aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
      aux_mat_1.shed_rows(i,i+1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape((V_net-1)*T_net,1);
      W.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
      
      gamma_ijt = gamma_ijtk(k);
      aux_mat_1 = gamma_ijt.subcube(i,0,0, i,i,T_net-1);
      aux_mat_2 = gamma_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
      aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
      aux_mat_1.shed_rows(i,i+1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape((V_net-1)*T_net,1);
      linpred.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
      
    }
    
    // X
    // Update matrix with covariate X
    X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
    for( t=0; t<T_net; t++ ) {
      X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = ab_ith_mat.submat( aux_uvec_3, aux_uvec_2+t );
    }
    X = X_all;
    X.shed_rows(T_net*i,T_net*(i+1)-1);
    X = repmat(X,K_net,1);
    X_sp = arma::sp_mat(X);
    
    // Constant term in linear predictor
    C = linpred - X_sp * trans(ab_ith_mat.row(i));
    
    // identifies valid obs
    valid_obs = find_finite(Y); // find valid elements
    
    // Keep only valid cases
    C_all = C;
    Y = Y.rows(valid_obs);
    W = W.rows(valid_obs);
    linpred = linpred.rows(valid_obs);
    C = C.rows(valid_obs);
    X_sp_valid = arma::sp_mat(X.rows(valid_obs));
    
    // marginal posterior covariance
    Omega_sp = arma::speye<arma::sp_mat>(W.n_rows,W.n_rows);
    Omega_sp.diag() = W;
    x_i_cov_inv = X_sp_valid.t() * Omega_sp * X_sp_valid;
    x_i_cov_inv = x_i_cov_inv + x_i_cov_prior ;
    x_i_cov = arma::inv_sympd(x_i_cov_inv);
    
    // marginal posterior mean
    Z = (Y-0.5)/W - C;
    x_i_mean = x_i_cov*(X_sp_valid.t() * (W % Z));
    
    // Sampling ab_ith_mat
    ab_ith_mat.row(i) = trans(arma::mvnrnd( x_i_mean , x_i_cov ));
    
    // Recalculate linpred with the new values of ab_ith_mat
    linpred = X_sp * trans(ab_ith_mat.row(i)) + C_all;
    // Redefine gamma_ijtk with the new values of linpred
    for( k=0; k<K_net; k++ ) {
      aux_mat_1 = linpred.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1);
      aux_mat_1.reshape(T_net,V_net-1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
      
      gamma_ijt = gamma_ijtk(k);
      gamma_ijt.subcube(i,0,0, i,i,T_net-1) = aux_mat_1.rows(0,i);
      gamma_ijt.subcube(i,i,0, V_net-1,i,T_net-1) = aux_mat_1.rows(i,V_net-1);
      gamma_ijtk(k)=gamma_ijt;
    }
    
  }
  
  // get ab_ith from ab_ith_mat
  ab_ith.slice(0)=ab_ith_mat;
  ab_ith = reshape(ab_ith,V_net,T_net,H_dim);
  
  return Rcpp::List::create( Rcpp::Named("ab_ith") = ab_ith,
                             Rcpp::Named("gamma_ijtk") = gamma_ijtk );
}


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_link_dir_GP_cpp( arma::cube a_ith,
                                             arma::cube b_ith,
                                             
                                             const arma::cube y_ijt,
                                             const arma::cube w_ijt,
                                             arma::cube gamma_ijt,
                                             
                                             const arma::mat ab_t_sigma_prior_inv,
                                             const bool lat_mean,
                                             arma::mat a_ith_bar,
                                             arma::mat b_ith_bar,
                                             const double sigma_ab_bar,
                                             
                                             
                                             const arma::colvec tau_h_send,
                                             const arma::colvec tau_h_receive ) {
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = a_ith.n_slices;
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  unsigned int h=0;
  unsigned int dir=0;
  arma::colvec aux_colvec = arma::ones<arma::colvec>(T_net);;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  a_ith = reshape(a_ith,V_net,T_net*H_dim,1);
  arma::mat a_ith_mat = a_ith.slice(0);
  
  b_ith = reshape(b_ith,V_net,T_net*H_dim,1);
  arma::mat b_ith_mat = b_ith.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec W = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec linpred = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec C_all = C;
  arma::colvec Z = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::sp_mat Omega_sp=arma::speye<arma::sp_mat>((V_net-1)*T_net,(V_net-1)*T_net);
  arma::mat tau_h_diag(tau_h_send.n_rows,tau_h_send.n_rows); tau_h_diag.eye();
  arma::mat x_i_cov_prior = kron( tau_h_diag, ab_t_sigma_prior_inv );
  arma::mat x_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat x_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  
  // prior mean vector for uv
  arma::colvec x_i_mean_prior = arma::zeros<arma::colvec>(T_net*H_dim);
  
  // Posterior mean vector for uv
  arma::colvec x_i_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  // Model matrix for all latent coordinates
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  
  // Model matrix that will be changed (V_net*2 times) for each node in the cycle, in their two positions: sender/receiver
  arma::mat X = X_all;
  arma::sp_mat X_sp = arma::sp_mat(X);
  arma::sp_mat X_sp_valid = X_sp;
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  
  // GP prior mean equal to zero
  if(!lat_mean){
    a_ith_bar.fill(0);
    b_ith_bar.fill(0);
  }
  arma::vec a_ith_bar_mean = arma::zeros<arma::vec>(1);
  arma::vec b_ith_bar_mean = arma::zeros<arma::vec>(1);
  double ab_ith_bar_var=1;
  
  
  for( dir=0; dir<2; dir++ ) {
    
    for( i=0; i<V_net; i++ ) {
      
      if( dir==0 ) {
        aux_mat_1 = y_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
        aux_mat_1.shed_row(i);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        Y = aux_mat_1;
        
        aux_mat_1 = w_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
        aux_mat_1.shed_row(i);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        W = aux_mat_1;
        
        aux_mat_1 = gamma_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
        aux_mat_1.shed_row(i);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        linpred = aux_mat_1;
        
        // X
        // Model matrix for all senders, made out of receiver coordinates
        if(i==0){ // The receiver coordinates doesn't change while sampling the senders
          X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
          for( t=0; t<T_net; t++ ) {
            X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = b_ith_mat.submat( aux_uvec_3, aux_uvec_2+t );
          }
          tau_h_diag.diag() = tau_h_send;
        }
        X = X_all;
        X.shed_rows(T_net*i,T_net*(i+1)-1);
        X_sp = arma::sp_mat(X);
        
        // Constant term in linear predictor
        C = linpred - X_sp * trans(a_ith_mat.row(i));
        
      } else if( dir==1 ) {
        aux_mat_1 = y_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_1.shed_row(i);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        Y = aux_mat_1;
        
        aux_mat_1 = w_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_1.shed_row(i);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        W = aux_mat_1;
        
        aux_mat_1 = gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_1.shed_row(i);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        linpred = aux_mat_1;
        
        // X
        // Model matrix for all receivers, made out of sender coordinates
        if(i==0){ // The sender coordinates doesn't change while sampling the receivers
          X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
          for( t=0; t<T_net; t++ ) {
            X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = a_ith_mat.submat( aux_uvec_3, aux_uvec_2+t );
          }
          tau_h_diag.diag() = tau_h_receive;
        }
        X = X_all;
        X.shed_rows(T_net*i,T_net*(i+1)-1);
        X_sp = arma::sp_mat(X);
        
        // Constant term in linear predictor
        C = linpred - X_sp * trans(b_ith_mat.row(i));
        
      } else {
        throw std::range_error("direction not supported");
      }
      
      // identifies valid obs
      valid_obs = find_finite(Y); // find valid elements
      
      // Keep only valid cases
      C_all = C;
      Y = Y.rows(valid_obs);
      W = W.rows(valid_obs);
      linpred = linpred.rows(valid_obs);
      C = C.rows(valid_obs);
      X_sp_valid = arma::sp_mat(X.rows(valid_obs));
      
      // Marginal posterior covariance
      Omega_sp = arma::speye<arma::sp_mat>(W.n_rows,W.n_rows);
      Omega_sp.diag() = W;
      x_i_cov_inv = X_sp_valid.t() * Omega_sp * X_sp_valid;
      x_i_cov_inv = x_i_cov_inv + x_i_cov_prior ;
      x_i_cov = arma::inv_sympd(x_i_cov_inv);
      
      // Z for Polson posterior mean
      Z = (Y-0.5)/W - C;
      
      if( dir==0 ) {
        
        // Marginal posterior mean
        x_i_mean_prior = kron( a_ith_bar.row(i).t(), aux_colvec );
        x_i_mean = x_i_cov * (X_sp_valid.t() * (W % Z) + x_i_cov_prior * x_i_mean_prior);
        
        // Sampling a_ith_mat
        a_ith_mat.row(i) = trans(arma::mvnrnd( x_i_mean , x_i_cov ));
        
        // Recalculate linpred with the new values of a_ith_mat
        linpred = X_sp * trans(a_ith_mat.row(i)) + C_all;
        // Redefine gamma_ijt with the new values of linpred
        aux_mat_1 = linpred;
        aux_mat_1.reshape(T_net,V_net-1);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
        gamma_ijt.subcube(i,0,0, i,V_net-1,T_net-1) = aux_mat_1;
        
      } else if( dir==1 ) {
        
        // Marginal posterior mean
        x_i_mean_prior = kron( b_ith_bar.row(i).t(), aux_colvec );
        x_i_mean = x_i_cov * (X_sp_valid.t() * (W % Z) + x_i_cov_prior * x_i_mean_prior);
        
        // Sampling b_ith_mat
        b_ith_mat.row(i) = trans(arma::mvnrnd( x_i_mean , x_i_cov ));
        
        // Recalculate linpred with the new values of a_ith_mat
        linpred = X_sp * trans(b_ith_mat.row(i)) + C_all;
        // Redefine gamma_ijt with the new values of linpred
        aux_mat_1 = linpred;
        aux_mat_1.reshape(T_net,V_net-1);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
        gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_1;
      } else {
        throw std::range_error("direction not supported");
      }
    }
  }
  
  // get ab_ith from ab_ith_mat
  a_ith.slice(0)=a_ith_mat;
  a_ith = reshape(a_ith,V_net,T_net,H_dim);
  
  b_ith.slice(0)=b_ith_mat;
  b_ith = reshape(b_ith,V_net,T_net,H_dim);
  
  // Sample GP prior mean 
  if(lat_mean){
    a_ith_bar.randn();
    
    b_ith_bar.randn();
    
    ab_ith_bar_var = 1/( accu(ab_t_sigma_prior_inv)+pow(sigma_ab_bar,-2) );
    
    aux_mat_1 = arma::ones<arma::mat>(1,T_net);
    aux_mat_1 = aux_mat_1*ab_t_sigma_prior_inv;
    for( i=0; i<V_net; i++ ) {
      // Rcpp::Rcout << " i=" << i << std::endl;
      for( h=0; h<H_dim; h++ ) {
        // Rcpp::Rcout << " h=" << h << std::endl;
        
        // Sender
        aux_mat_2 = a_ith.subcube(i,0,h, i,T_net-1,h);
        a_ith_bar_mean = ab_ith_bar_var * (aux_mat_1*aux_mat_2.t());
        a_ith_bar(i,h) = a_ith_bar_mean(0) + a_ith_bar(i,h) * sqrt(ab_ith_bar_var);
        
        // Receiver
        aux_mat_2 = b_ith.subcube(i,0,h, i,T_net-1,h);
        b_ith_bar_mean = ab_ith_bar_var * (aux_mat_1*aux_mat_2.t());
        b_ith_bar(i,h) = b_ith_bar_mean(0) + b_ith_bar(i,h) * sqrt(ab_ith_bar_var);
      }
    }
    
  }
  
  return Rcpp::List::create( Rcpp::Named("a_ith") = a_ith,
                             Rcpp::Named("b_ith") = b_ith,
                             Rcpp::Named("gamma_ijt") = gamma_ijt,
                             Rcpp::Named("a_ith_bar") = a_ith_bar,
                             Rcpp::Named("b_ith_bar") = b_ith_bar );
}


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_shared_link_dir_GP_cpp( arma::cube a_ith_shared,
                                                    arma::cube b_ith_shared,
                                                    
                                                    const arma::field<arma::cube> y_ijtk,
                                                    const arma::field<arma::cube> w_ijtk,
                                                    arma::field<arma::cube> gamma_ijtk,
                                                    
                                                    const arma::mat ab_t_sigma_prior_inv,
                                                    const bool lat_mean,
                                                    arma::mat a_ith_shared_bar,
                                                    arma::mat b_ith_shared_bar,
                                                    const double sigma_ab_bar,
                                                    
                                                    const arma::colvec tau_h_shared_send,
                                                    const arma::colvec tau_h_shared_receive ) {
  
  // Single layer objects
  arma::cube y_ijt = y_ijtk(0);
  arma::cube w_ijt = w_ijtk(0);
  arma::cube gamma_ijt = gamma_ijtk(0);
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = a_ith_shared.n_slices;
  unsigned int K_net = y_ijtk.n_rows;
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int k=0;
  unsigned int t=0;
  unsigned int h=0;
  unsigned int dir=0;
  arma::colvec aux_colvec = arma::ones<arma::colvec>(T_net);
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  a_ith_shared = reshape(a_ith_shared,V_net,T_net*H_dim,1);
  arma::mat a_ith_shared_mat = a_ith_shared.slice(0);
  
  b_ith_shared = reshape(b_ith_shared,V_net,T_net*H_dim,1);
  arma::mat b_ith_shared_mat = b_ith_shared.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec W = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec linpred = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec C_all = C;
  arma::colvec Z = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  arma::sp_mat Omega_sp=arma::speye<arma::sp_mat>((V_net-1)*T_net*K_net,(V_net-1)*T_net*K_net);
  arma::mat tau_h_diag(tau_h_shared_send.n_rows,tau_h_shared_send.n_rows); tau_h_diag.eye();
  arma::mat x_i_cov_prior = kron( tau_h_diag, ab_t_sigma_prior_inv );
  arma::mat x_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat x_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  
  // prior mean vector for ab
  arma::colvec x_i_mean_prior = arma::zeros<arma::colvec>(T_net*H_dim);
  
  // Posterior mean vector for ab
  arma::colvec x_i_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  // Model matrix for all senders, made out of receiver coordinates
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  
  // Model matrix that will be changed (V_net*2 times) for each node in the cycle, in their two positions: sender/receiver
  arma::mat X = X_all;
  arma::sp_mat X_sp = arma::sp_mat(X);
  arma::sp_mat X_sp_valid = X_sp;
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  
  // GP prior mean equal to zero
  if(!lat_mean){
    a_ith_shared_bar.fill(0);
    b_ith_shared_bar.fill(0);
  }
  arma::vec a_ith_shared_bar_mean = arma::zeros<arma::vec>(1);
  arma::vec b_ith_shared_bar_mean = arma::zeros<arma::vec>(1);
  double ab_ith_shared_bar_var=1;
  
  for( dir=0; dir<2; dir++ ) {
    // Rcpp::Rcout << "dir=" << dir << std::endl;
    
    for( i=0; i<V_net; i++ ) {
      // Rcpp::Rcout << " i=" << i << std::endl;
      
      Y = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
      W = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
      linpred = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
      C = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
      Z = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
      
      if( dir==0 ) {
        // Updating Y, W, linpred as vectors, lower triangular adjacency
        for( k=0; k<K_net; k++ ) {
          
          y_ijt = y_ijtk(k);
          aux_mat_1 = y_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
          aux_mat_1.shed_row(i);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          Y.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
          
          w_ijt = w_ijtk(k);
          aux_mat_1 = w_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
          aux_mat_1.shed_row(i);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          W.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
          
          gamma_ijt = gamma_ijtk(k);
          aux_mat_1 = gamma_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
          aux_mat_1.shed_row(i);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          linpred.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
        }
        // X
        if(i==0){
          // The matrix of covariates X will only be updated at the beginning of the cycle
          // as it does not change with the newly sampled values
          X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
          for( t=0; t<T_net; t++ ) {
            X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = b_ith_shared_mat.submat( aux_uvec_3, aux_uvec_2+t );
          }
          tau_h_diag.diag() = tau_h_shared_send;
        }
        X = X_all;
        X.shed_rows(T_net*i,T_net*(i+1)-1);
        X = repmat(X,K_net,1);
        X_sp = arma::sp_mat(X);
        
        // Constant term in linear predictor
        C = linpred - X_sp * trans(a_ith_shared_mat.row(i));
        
      } else if( dir==1 ) {
        // Updating Y, W, linpred as vectors, upper triangular adjacency
        for( k=0; k<K_net; k++ ) {
          y_ijt = y_ijtk(k);
          aux_mat_1 = y_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
          aux_mat_1.shed_row(i);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          Y.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
          
          w_ijt = w_ijtk(k);
          aux_mat_1 = w_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
          aux_mat_1.shed_row(i);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          W.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
          
          gamma_ijt = gamma_ijtk(k);
          aux_mat_1 = gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
          aux_mat_1.shed_row(i);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          linpred.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
        }
        // X
        if(i==0){
          // The full matrix of covariates X will only be updated at the beginning of the cycle
          // as it does not change with the newly sampled values
          X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
          for( t=0; t<T_net; t++ ) {
            X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = a_ith_shared_mat.submat( aux_uvec_3, aux_uvec_2+t );
          }
          tau_h_diag.diag() = tau_h_shared_receive;
        }
        X = X_all;
        X.shed_rows(T_net*i,T_net*(i+1)-1);
        X = repmat(X,K_net,1);
        X_sp = arma::sp_mat(X);
        
        // Constant term in linear predictor
        C = linpred - X_sp * trans(b_ith_shared_mat.row(i));
        
      } else {
        throw std::range_error("direction not supported");
      }
      
      // identifies valid obs
      valid_obs = find_finite(Y); // find valid elements
      
      // Keep only valid cases
      C_all = C;
      Y = Y.rows(valid_obs);
      W = W.rows(valid_obs);
      linpred = linpred.rows(valid_obs);
      C = C.rows(valid_obs);
      X_sp_valid = arma::sp_mat(X.rows(valid_obs));
      
      // Marginal posterior covariance
      Omega_sp = arma::speye<arma::sp_mat>(W.n_rows,W.n_rows);
      Omega_sp.diag() = W;
      x_i_cov_inv = X_sp_valid.t() * Omega_sp * X_sp_valid;
      x_i_cov_inv = x_i_cov_inv + x_i_cov_prior ;
      x_i_cov = arma::inv_sympd(x_i_cov_inv);
      
      // Z for Polson posterior mean
      Z = (Y-0.5)/W - C;
      
      if( dir==0 ) {
        
        // Marginal posterior mean
        x_i_mean_prior = kron( a_ith_shared_bar.row(i).t(), aux_colvec );
        x_i_mean = x_i_cov * (X_sp_valid.t() * (W % Z) + x_i_cov_prior * x_i_mean_prior);
        
        // Sampling a_ith_shared_mat
        a_ith_shared_mat.row(i) = trans(arma::mvnrnd( x_i_mean , x_i_cov ));
        
        // Recalculate linpred with the new values of a_ith_shared_mat
        linpred = X_sp * trans(a_ith_shared_mat.row(i)) + C_all;
        // Redefine gamma_ijt with the new values of linpred
        for( k=0; k<K_net; k++ ) {
          aux_mat_1 = linpred.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1);
          aux_mat_1.reshape(T_net,(V_net-1));
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
          
          gamma_ijt = gamma_ijtk(k);
          gamma_ijt.subcube(i,0,0, i,V_net-1,T_net-1) = aux_mat_1;
          gamma_ijtk(k)= gamma_ijt;
        }
        
      } else if( dir==1 ) {
        
        // Marginal posterior mean
        x_i_mean_prior = kron( b_ith_shared_bar.row(i).t(), aux_colvec );
        x_i_mean = x_i_cov * (X_sp_valid.t() * (W % Z) + x_i_cov_prior * x_i_mean_prior);
        
        // Sampling b_ith_shared_mat
        b_ith_shared_mat.row(i) = trans(arma::mvnrnd( x_i_mean , x_i_cov ));
        
        // Recalculate linpred with the new values of b_ith_shared_mat
        linpred = X_sp * trans(b_ith_shared_mat.row(i)) + C_all;
        // Redefine gamma_ijt with the new values of linpred
        for( k=0; k<K_net; k++ ) {
          aux_mat_1 = linpred.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1);
          aux_mat_1.reshape(T_net,(V_net-1));
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
          
          gamma_ijt = gamma_ijtk(k);
          gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_1;
          gamma_ijtk(k)= gamma_ijt;
        }
      } else {
        throw std::range_error("direction not supported");
      }
    }
  }
  
  // get ab_ith from ab_ith_mat
  a_ith_shared.slice(0)=a_ith_shared_mat;
  a_ith_shared = reshape(a_ith_shared,V_net,T_net,H_dim);
  
  b_ith_shared.slice(0)=b_ith_shared_mat;
  b_ith_shared = reshape(b_ith_shared,V_net,T_net,H_dim);
  
  // Sample GP prior mean 
  if(lat_mean){
    a_ith_shared_bar.randn();
    
    b_ith_shared_bar.randn();
    
    ab_ith_shared_bar_var = 1/( accu(ab_t_sigma_prior_inv)+pow(sigma_ab_bar,-2) );
    
    aux_mat_1 = arma::ones<arma::mat>(1,T_net);
    aux_mat_1 = aux_mat_1*ab_t_sigma_prior_inv;
    for( i=0; i<V_net; i++ ) {
      // Rcpp::Rcout << " i=" << i << std::endl;
      for( h=0; h<H_dim; h++ ) {
        // Rcpp::Rcout << " h=" << h << std::endl;
        
        // Sender
        aux_mat_2 = a_ith_shared.subcube(i,0,h, i,T_net-1,h);
        a_ith_shared_bar_mean = ab_ith_shared_bar_var * (aux_mat_1*aux_mat_2.t());
        a_ith_shared_bar(i,h) = a_ith_shared_bar_mean(0) + a_ith_shared_bar(i,h) * sqrt(ab_ith_shared_bar_var);
        
        // Receiver
        aux_mat_2 = b_ith_shared.subcube(i,0,h, i,T_net-1,h);
        b_ith_shared_bar_mean = ab_ith_shared_bar_var * (aux_mat_1*aux_mat_2.t());
        b_ith_shared_bar(i,h) = b_ith_shared_bar_mean(0) + b_ith_shared_bar(i,h) * sqrt(ab_ith_shared_bar_var);
      }
    }
    
  }
  
  return Rcpp::List::create( Rcpp::Named("a_ith_shared") = a_ith_shared,
                             Rcpp::Named("b_ith_shared") = b_ith_shared,
                             Rcpp::Named("gamma_ijtk") = gamma_ijtk,
                             Rcpp::Named("a_ith_shared_bar") = a_ith_shared_bar,
                             Rcpp::Named("b_ith_shared_bar") = b_ith_shared_bar );
}


// [[Rcpp::export]]
Rcpp::List sample_add_eff_it_link_cpp( arma::colvec sp_it,
                                       const arma::mat sp_t_cov_prior_inv,
                                       const arma::cube y_ijt,
                                       arma::cube w_ijt,
                                       arma::cube gamma_ijt,
                                       const bool directed=false ) {
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  unsigned int dir=0;
  arma::uvec aux_uvec_1;
  arma::colvec aux_vec;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::mat aux_mat_3;
  arma::cube aux_cube_1;
  
  // Set irrelevants elements of w_ijt equal to zero //
  if( !directed ){
    for( t=0; t<T_net; t++ ) {
      w_ijt.slice(t) = trimatl(w_ijt.slice(t),-1);
      w_ijt.slice(t).elem( find_nonfinite(y_ijt.slice(t)) ).zeros(); // set W=0 where Y is na
    }
  } else {
    for( t=0; t<T_net; t++ ) {
      w_ijt.slice(t).diag().zeros();
      w_ijt.slice(t).elem( find_nonfinite(y_ijt.slice(t)) ).zeros(); // set W=0 where Y is na
    }
  }
  
  // Objects for the calculation
  arma::colvec Y = arma::zeros<arma::colvec>(1);
  arma::colvec W = arma::zeros<arma::colvec>(1);
  arma::colvec linpred = arma::zeros<arma::colvec>(1);
  
  arma::colvec C = arma::zeros<arma::colvec>(1);
  arma::colvec C_valid = C;
  arma::colvec Z = arma::zeros<arma::colvec>(1);
  
  arma::sp_mat Omega_sp;
  
  // Prior covariance
  arma::mat sp_it_cov_prior_inv = kron( sp_t_cov_prior_inv, arma::eye<arma::mat>(V_net,V_net) );
  
  // Design matrix //
  arma::mat X = arma::zeros<arma::mat>(1,1);
  // Setting design matrix //
  aux_mat_1 = arma::eye<arma::mat>(V_net,V_net); // create diagonal matrix
  if( !directed ){
    aux_mat_1.insert_cols(aux_mat_1.n_cols,(T_net-1)*V_net); // insert zero columns to end up having TV columns
  } else {
    aux_mat_1.insert_cols(aux_mat_1.n_cols,(2*T_net-1)*V_net); // insert zero columns to end up having 2TV columns
  }
  // aux_mat_1 will be the mold duplicated
  aux_mat_3 = arma::zeros<arma::mat>(1,aux_mat_1.n_cols);
  for( i=0; i<V_net; i++ ) {
    aux_mat_2 = aux_mat_1;
    if( !directed ){
      aux_mat_2.col(i) = arma::ones<arma::mat>(V_net,1); // set i-th column to 1
      aux_mat_2.shed_rows(0,i);
    } else {
      aux_mat_2.col(T_net*V_net+i) = arma::ones<arma::mat>(V_net,1); // set TV+i-th column to 1
      aux_mat_2.shed_row(i);
    }
    aux_mat_3.insert_rows(aux_mat_3.n_rows,aux_mat_2);
  }
  aux_mat_3.shed_row(0);
  
  X = aux_mat_3;
  for( t=1; t<T_net; t++ ) {
    aux_mat_3.shed_cols(aux_mat_3.n_cols-V_net,aux_mat_3.n_cols-1);
    aux_mat_3.insert_cols(0,V_net);
    X.insert_rows(X.n_rows,aux_mat_3);
  }
  arma::sp_mat X_sp = arma::sp_mat(X);
  arma::sp_mat X_sp_valid = X_sp;
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  
  // Marginal posterior covariance matrix for sp_it
  arma::mat sp_it_cov_inv=arma::eye<arma::mat>(T_net*V_net,T_net*V_net);
  arma::mat sp_it_cov;
  arma::colvec sp_it_mean = arma::zeros<arma::colvec>(2*V_net*T_net);
  
  if( !directed ){
    
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = y_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
    Y = aux_mat_1;
    
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = w_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
    W = aux_mat_1;
    
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = gamma_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
    linpred = aux_mat_1;
    
  } else {
    
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = y_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
      aux_mat_2.shed_row(i);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
    Y = aux_mat_1;
    
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = w_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
      aux_mat_2.shed_row(i);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
    W = aux_mat_1;
    
    aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
      aux_mat_2.shed_row(i);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
    linpred = aux_mat_1;
    
  }
  
  // identifies valid obs
  valid_obs = find_finite(Y); // find valid elements
  
  // Keep only valid cases
  Y = Y.rows(valid_obs);
  W = W.rows(valid_obs);
  X_sp_valid = arma::sp_mat(X.rows(valid_obs));
  
  for( dir=0; dir<2; dir++ ) { // dir=0 samples p ; dir=1 samples s
    if( (dir==1) | (directed&(dir==0)) ) { // if directed, sample s and p. if undirected sample only s
      // Rcpp::Rcout << "dir=" << dir << std::endl;
      
      // marginal posterior covariance
      if(false){
        // // Way 1:
        // Omega_sp=arma::speye<arma::sp_mat>(W.n_rows/2,W.n_rows/2);
        // Omega_sp.diag() = W;
        // if( dir==0 ){
        //   sp_it_cov_inv = X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1).t() * Omega_sp * X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1);
        // } else if( dir==1 ){
        //   sp_it_cov_inv = X_sp_valid.cols(0,T_net*V_net-1).t() * Omega_sp * X_sp_valid.cols(0,T_net*V_net-1);
        // }
      } else {
        // Way 2:
        sp_it_cov_inv=arma::zeros<arma::mat>(T_net*V_net,T_net*V_net);
        if(!directed){
          for( t=0; t<T_net; t++ ) {
            aux_mat_1 = symmatl(w_ijt.slice(t)); // symmetric matrix from w_ijt, reflecting lower triangular matrix
            aux_mat_1.diag() = sum(aux_mat_1,1); // diagonal equal to the sum by row or column
            sp_it_cov_inv.submat(t*V_net,t*V_net,(t+1)*V_net-1,(t+1)*V_net-1) = aux_mat_1;
          }
        } else {
          aux_mat_1 = sum(w_ijt,dir);
          aux_mat_1.reshape(T_net*V_net,1);
          sp_it_cov_inv.diag()=aux_mat_1;
        }
      }
      // adding prior covariance
      sp_it_cov_inv = sp_it_cov_inv + sp_it_cov_prior_inv;
      sp_it_cov = arma::inv_sympd(sp_it_cov_inv);
      
      
      if( dir==0 ){
        // marginal posterior mean
        C = linpred - (X_sp.cols(T_net*V_net,2*T_net*V_net-1) * sp_it.rows(T_net*V_net,2*T_net*V_net-1));
        C_valid = C.rows(valid_obs);
        Z = (Y-0.5)/W - C_valid;
        sp_it_mean = sp_it_cov * (X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1).t() * (W % Z));
        
        // Sampling sp_it
        sp_it.rows(T_net*V_net,2*T_net*V_net-1) = arma::mvnrnd( sp_it_mean , sp_it_cov );
        
        // Recalculate linpred with the new values of sp_it
        linpred = X_sp.cols(T_net*V_net,2*T_net*V_net-1) * sp_it.rows(T_net*V_net,2*T_net*V_net-1) + C;
        
      } else if( dir==1 ){
        // marginal posterior mean
        C = linpred - (X_sp.cols(0,T_net*V_net-1) * sp_it.rows(0,T_net*V_net-1));
        C_valid = C.rows(valid_obs);
        Z = (Y-0.5)/W - C_valid;
        sp_it_mean = sp_it_cov * (X_sp_valid.cols(0,T_net*V_net-1).t() * (W % Z));
        
        // Sampling sp_it
        sp_it.rows(0,T_net*V_net-1) = arma::mvnrnd( sp_it_mean , sp_it_cov );
        
        // Recalculate linpred with the new values of sp_it
        linpred = X_sp.cols(0,T_net*V_net-1) * sp_it.rows(0,T_net*V_net-1) + C;
        
      }
    }
  }
  
  // Redefine gamma_ijt with the new values of linpred
  if( !directed ){
    aux_mat_1 = linpred;
    aux_mat_1.reshape(V_net*(V_net-1)/2,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = aux_mat_1.rows((i-1)*V_net-((i-1)*i)/2,i*V_net-(i*(i+1))/2-1);
      gamma_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1) = aux_mat_2;
    }
  } else {
    aux_mat_1 = linpred;
    aux_mat_1.reshape(V_net*(V_net-1),T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = aux_mat_1.rows(i*(V_net-1),(i+1)*(V_net-1)-1);
      aux_mat_2.insert_rows(i,arma::zeros<arma::mat>(1,T_net));
      gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_2;
    }
  }
  
  return Rcpp::List::create( Rcpp::Named("sp_it") = sp_it,
                             Rcpp::Named("gamma_ijt") = gamma_ijt );
}


// [[Rcpp::export]]
Rcpp::List sample_add_eff_it_shared_link_cpp( arma::colvec sp_it,
                                              const arma::mat sp_t_cov_prior_inv,
                                              const arma::field<arma::cube> y_ijtk,
                                              arma::field<arma::cube> w_ijtk,
                                              arma::field<arma::cube> gamma_ijtk,
                                              const bool directed=false ) {
  
  arma::cube y_ijt = y_ijtk(0);
  arma::cube w_ijt = w_ijtk(0);
  arma::cube gamma_ijt = gamma_ijtk(0);
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int K_net = y_ijtk.n_rows;
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  unsigned int k=0;
  unsigned int dir=0;
  arma::uvec aux_uvec_1;
  arma::colvec aux_vec;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::mat aux_mat_3;
  arma::cube aux_cube_1;
  
  // Set irrelevants elements of w_ijt equal to zero //
  for( k=0; k<K_net; k++ ) {
    w_ijt = w_ijtk(k);
    if( !directed ){
      for( t=0; t<T_net; t++ ) {
        w_ijt.slice(t) = trimatl(w_ijt.slice(t),-1);
        w_ijt.slice(t).elem( find_nonfinite(y_ijt.slice(t)) ).zeros(); // set W=0 where Y is na
      }
    } else {
      for( t=0; t<T_net; t++ ) {
        w_ijt.slice(t).diag().zeros();
        w_ijt.slice(t).elem( find_nonfinite(y_ijt.slice(t)) ).zeros(); // set W=0 where Y is na
      }
    }
    w_ijtk(k)=w_ijt;
  }
  
  // Objects for the calculation
  arma::colvec Y = arma::zeros<arma::colvec>(1);
  arma::colvec W = arma::zeros<arma::colvec>(1);
  arma::colvec linpred = arma::zeros<arma::colvec>(1);
  arma::colvec linpred_all = linpred;
  
  arma::colvec C = arma::zeros<arma::colvec>(1);
  arma::colvec C_valid = C;
  arma::colvec Z = arma::zeros<arma::colvec>(1);
  
  // arma::sp_mat Omega_sp;
  
  // Prior covariance
  arma::mat sp_it_cov_prior_inv = kron( sp_t_cov_prior_inv, arma::eye<arma::mat>(V_net,V_net) );
  
  // Design matrix //
  arma::mat X = arma::zeros<arma::mat>(1,1);
  // Setting design matrix //
  aux_mat_1 = arma::eye<arma::mat>(V_net,V_net); // create diagonal matrix
  if( !directed ){
    aux_mat_1.insert_cols(aux_mat_1.n_cols,(T_net-1)*V_net); // insert zero columns to end up having TV columns
  } else {
    aux_mat_1.insert_cols(aux_mat_1.n_cols,(2*T_net-1)*V_net); // insert zero columns to end up having 2TV columns
  }
  // aux_mat_1 will be the mold duplicated
  aux_mat_3 = arma::zeros<arma::mat>(1,aux_mat_1.n_cols);
  for( i=0; i<V_net; i++ ) {
    aux_mat_2 = aux_mat_1;
    if( !directed ){
      aux_mat_2.col(i) = arma::ones<arma::mat>(V_net,1); // set i-th column to 1
      aux_mat_2.shed_rows(0,i);
    } else {
      aux_mat_2.col(T_net*V_net+i) = arma::ones<arma::mat>(V_net,1); // set TV+i-th column to 1
      aux_mat_2.shed_row(i);
    }
    aux_mat_3.insert_rows(aux_mat_3.n_rows,aux_mat_2);
  }
  aux_mat_3.shed_row(0);
  
  X = aux_mat_3;
  for( t=1; t<T_net; t++ ) {
    aux_mat_3.shed_cols(aux_mat_3.n_cols-V_net,aux_mat_3.n_cols-1);
    aux_mat_3.insert_cols(0,V_net);
    X.insert_rows(X.n_rows,aux_mat_3);
  }
  X = repmat(X,K_net,1);
  arma::sp_mat X_sp = arma::sp_mat(X);
  arma::sp_mat X_sp_valid = X_sp;
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  
  // Marginal posterior covariance matrix for sp_it
  arma::mat sp_it_cov_inv=arma::zeros<arma::mat>(T_net*V_net,T_net*V_net);
  arma::mat sp_it_cov;
  arma::colvec sp_it_mean = arma::zeros<arma::colvec>(2*V_net*T_net);
  
  
  // Create vectors Y, W and linpred tht serve as response and  //
  if( !directed ){
    Y = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
    W = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
    linpred = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
  } else {
    Y = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
    W = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
    linpred = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
  }
  
  
  for( k=0; k<K_net; k++ ) {
    // Rcpp::Rcout << "k=" << k << std::endl;
    y_ijt = y_ijtk(k);
    w_ijt = w_ijtk(k);
    gamma_ijt = gamma_ijtk(k);
    
    if( !directed ){
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=1; i<V_net; i++ ) {
        aux_mat_2 = y_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
      Y.rows(k*T_net*V_net*(V_net-1)/2, (k+1)*T_net*V_net*(V_net-1)/2-1) = aux_mat_1;
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=1; i<V_net; i++ ) {
        aux_mat_2 = w_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
      W.rows(k*T_net*V_net*(V_net-1)/2, (k+1)*T_net*V_net*(V_net-1)/2-1) = aux_mat_1;
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=1; i<V_net; i++ ) {
        aux_mat_2 = gamma_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
      linpred.rows(k*T_net*V_net*(V_net-1)/2, (k+1)*T_net*V_net*(V_net-1)/2-1) = aux_mat_1;
      
    } else {
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=0; i<V_net; i++ ) {
        aux_mat_2 = y_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_2.shed_row(i);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
      Y.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1) = aux_mat_1;
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=0; i<V_net; i++ ) {
        aux_mat_2 = w_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_2.shed_row(i);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
      W.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1) = aux_mat_1;
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=0; i<V_net; i++ ) {
        aux_mat_2 = gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_2.shed_row(i);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
      linpred.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1) = aux_mat_1;
      
    }
  }
  
  // identifies valid obs
  valid_obs = find_finite(Y); // find valid elements
  
  // Keep only valid cases
  Y = Y.rows(valid_obs);
  W = W.rows(valid_obs);
  X_sp_valid = arma::sp_mat(X.rows(valid_obs));
  
  for( dir=0; dir<2; dir++ ) { // dir=0 samples p ; dir=1 samples s
    if( (dir==1) | (directed&(dir==0)) ) { // if directed, sample s and p. if undirected sample only s
      // Rcpp::Rcout << "dir=" << dir << std::endl;
      
      // marginal posterior covariance
      if(false){
        // // Way 1:
        // Omega_sp=arma::speye<arma::sp_mat>(W.n_rows,W.n_rows);
        // Omega_sp.diag() = W;
        // if( dir==0 ){
        //   sp_it_cov_inv = X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1).t() * Omega_sp * X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1);
        // } else if( dir==1 ){
        //   sp_it_cov_inv = X_sp_valid.cols(0,T_net*V_net-1).t() * Omega_sp * X_sp_valid.cols(0,T_net*V_net-1);
        // }
      } else {
        // Way 2:
        sp_it_cov_inv=arma::zeros<arma::mat>(T_net*V_net,T_net*V_net);
        if(!directed){
          for( k=0; k<K_net; k++ ) {
            for( t=0; t<T_net; t++ ) {
              aux_mat_1 = symmatl(w_ijtk(k).slice(t)); // symmetric matrix from w_ijt, reflecting lower triangular matrix
              aux_mat_1.diag() = sum(aux_mat_1,1); // diagonal equal to the sum by row or column
              sp_it_cov_inv.submat(t*V_net,t*V_net,(t+1)*V_net-1,(t+1)*V_net-1) += aux_mat_1;
            }
          }
        } else {
          for( k=0; k<K_net; k++ ) {
            aux_mat_1 = sum(w_ijtk(k),dir);
            aux_mat_1.reshape(T_net*V_net,1);
            sp_it_cov_inv.diag()+=aux_mat_1;
          }
        }
      }
      // adding prior covariance
      sp_it_cov_inv = sp_it_cov_inv + sp_it_cov_prior_inv;
      sp_it_cov = arma::inv_sympd(sp_it_cov_inv);
      
      if( dir==0 ){
        // marginal posterior mean
        C = linpred - (X_sp.cols(T_net*V_net,2*T_net*V_net-1) * sp_it.rows(T_net*V_net,2*T_net*V_net-1));
        C_valid = C.rows(valid_obs);
        Z = (Y-0.5)/W - C_valid;
        sp_it_mean = sp_it_cov * (X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1).t() * (W % Z));
        
        // Sampling sp_it
        sp_it.rows(T_net*V_net,2*T_net*V_net-1) = arma::mvnrnd( sp_it_mean , sp_it_cov );
        
        // Recalculate linpred with the new values of sp_it
        linpred = X_sp.cols(T_net*V_net,2*T_net*V_net-1) * sp_it.rows(T_net*V_net,2*T_net*V_net-1) + C;
      } else if( dir==1 ){
        // marginal posterior mean
        C = linpred - (X_sp.cols(0,T_net*V_net-1) * sp_it.rows(0,T_net*V_net-1));
        C_valid = C.rows(valid_obs);
        Z = (Y-0.5)/W - C_valid;
        sp_it_mean = sp_it_cov * (X_sp_valid.cols(0,T_net*V_net-1).t() * (W % Z));
        
        // Sampling sp_it
        sp_it.rows(0,T_net*V_net-1) = arma::mvnrnd( sp_it_mean , sp_it_cov );
        
        // Recalculate linpred with the new values of sp_it
        linpred = X_sp.cols(0,T_net*V_net-1) * sp_it.rows(0,T_net*V_net-1) + C;
      }
    }
  }
  
  // Redefine gamma_ijt with the new values of linpred
  for( k=0; k<K_net; k++ ) {
    // Rcpp::Rcout << "k=" << k << std::endl;
    gamma_ijt = gamma_ijtk(k);
    
    if( !directed ){
      aux_mat_1 = linpred.rows(k*T_net*V_net*(V_net-1)/2, (k+1)*T_net*V_net*(V_net-1)/2-1);
      aux_mat_1.reshape(V_net*(V_net-1)/2,T_net);
      for( i=1; i<V_net; i++ ) {
        aux_mat_2 = aux_mat_1.rows((i-1)*V_net-((i-1)*i)/2,i*V_net-(i*(i+1))/2-1);
        gamma_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1) = aux_mat_2;
      }
    } else {
      aux_mat_1 = linpred.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1);
      aux_mat_1.reshape(V_net*(V_net-1),T_net);
      for( i=0; i<V_net; i++ ) {
        aux_mat_2 = aux_mat_1.rows(i*(V_net-1),(i+1)*(V_net-1)-1);
        aux_mat_2.insert_rows(i,arma::zeros<arma::mat>(1,T_net));
        gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_2;
      }
    }
    
    gamma_ijtk(k) = gamma_ijt;
  }
  
  return Rcpp::List::create( Rcpp::Named("sp_it") = sp_it,
                             Rcpp::Named("gamma_ijtk") = gamma_ijtk );
}


// [[Rcpp::export]]
Rcpp::List sample_coeff_tp_link_cpp( arma::mat beta_tp,
                                     const arma::mat beta_t_cov_prior_inv,
                                     const arma::field<arma::cube> y_ijtk,
                                     arma::field<arma::cube> w_ijtk,
                                     arma::field<arma::cube> gamma_ijtk,
                                     const arma::mat x_ijtkp_mat,
                                     const bool directed=false ) {
  
  // This function block sample
  // beta_tp = (beta_tp[1,1],...,beta_tp[T_net,1],...,beta_tp[1,P_pred],...,beta_tp[T_net,P_pred])
  
  arma::cube y_ijt = y_ijtk(0);
  arma::cube w_ijt = w_ijtk(0);
  arma::cube gamma_ijt = gamma_ijtk(0);
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int K_net = y_ijtk.n_rows;
  unsigned int P_pred = beta_tp.n_cols;
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  unsigned int k=0;
  arma::uvec aux_uvec_1;
  arma::colvec aux_vec;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::mat aux_mat_3;
  arma::cube aux_cube_1;
  
  // reshape beta
  beta_tp.reshape(T_net*P_pred,1);
  
  // Set irrelevants elements of w_ijt equal to zero //
  for( k=0; k<K_net; k++ ) {
    w_ijt = w_ijtk(k);
    if( !directed ){
      for( t=0; t<T_net; t++ ) {
        w_ijt.slice(t) = trimatl(w_ijt.slice(t),-1);
        w_ijt.slice(t).elem( find_nonfinite(y_ijt.slice(t)) ).zeros(); // set W=0 where Y is na
      }
    } else {
      for( t=0; t<T_net; t++ ) {
        w_ijt.slice(t).diag().zeros();
        w_ijt.slice(t).elem( find_nonfinite(y_ijt.slice(t)) ).zeros(); // set W=0 where Y is na
      }
    }
    w_ijtk(k)=w_ijt;
  }
  
  // Objects for the calculation
  arma::colvec Y = arma::zeros<arma::colvec>(1);
  arma::colvec W = arma::zeros<arma::colvec>(1);
  arma::colvec linpred = arma::zeros<arma::colvec>(1);
  arma::colvec linpred_all = linpred;
  
  arma::colvec C = arma::zeros<arma::colvec>(1);
  arma::colvec C_valid = C;
  arma::colvec Z = arma::zeros<arma::colvec>(1);
  
  arma::sp_mat Omega_sp;
  
  // Prior covariance
  arma::mat beta_tp_cov_prior_inv = kron( beta_t_cov_prior_inv, arma::eye<arma::mat>(P_pred,P_pred) );
  
  arma::sp_mat X_sp = arma::sp_mat(x_ijtkp_mat);
  arma::sp_mat X_sp_valid = X_sp;
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  
  // Marginal posterior covariance matrix for sp_it
  arma::mat beta_tp_cov_inv=arma::zeros<arma::mat>(T_net*V_net,T_net*V_net);
  arma::mat beta_tp_cov;
  arma::colvec beta_tp_mean = arma::zeros<arma::colvec>(P_pred*T_net);
  
  
  // Create vectors Y, W and linpred tht serve as response and  //
  if( !directed ){
    Y = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
    W = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
    linpred = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
  } else {
    Y = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
    W = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
    linpred = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
  }
  
  
  for( k=0; k<K_net; k++ ) {
    // Rcpp::Rcout << "k=" << k << std::endl;
    y_ijt = y_ijtk(k);
    w_ijt = w_ijtk(k);
    gamma_ijt = gamma_ijtk(k);
    
    if( !directed ){
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=1; i<V_net; i++ ) {
        aux_mat_2 = y_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
      Y.rows(k*T_net*V_net*(V_net-1)/2, (k+1)*T_net*V_net*(V_net-1)/2-1) = aux_mat_1;
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=1; i<V_net; i++ ) {
        aux_mat_2 = w_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
      W.rows(k*T_net*V_net*(V_net-1)/2, (k+1)*T_net*V_net*(V_net-1)/2-1) = aux_mat_1;
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=1; i<V_net; i++ ) {
        aux_mat_2 = gamma_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
      linpred.rows(k*T_net*V_net*(V_net-1)/2, (k+1)*T_net*V_net*(V_net-1)/2-1) = aux_mat_1;
      
    } else {
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=0; i<V_net; i++ ) {
        aux_mat_2 = y_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_2.shed_row(i);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
      Y.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1) = aux_mat_1;
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=0; i<V_net; i++ ) {
        aux_mat_2 = w_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_2.shed_row(i);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
      W.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1) = aux_mat_1;
      
      aux_mat_1 = arma::zeros<arma::mat>(1,T_net);
      for( i=0; i<V_net; i++ ) {
        aux_mat_2 = gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_2.shed_row(i);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
      linpred.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1) = aux_mat_1;
      
    }
  }
  
  // identifies valid obs
  valid_obs = find_finite(Y); // find valid elements
  
  // Keep only valid cases
  Y = Y.rows(valid_obs);
  W = W.rows(valid_obs);
  X_sp_valid = arma::sp_mat(x_ijtkp_mat.rows(valid_obs));
  
  // marginal posterior covariance
  Omega_sp=arma::speye<arma::sp_mat>(W.n_rows,W.n_rows);
  Omega_sp.diag() = W;
  beta_tp_cov_inv = X_sp_valid.t() * (Omega_sp * X_sp_valid);
  // adding prior covariance
  beta_tp_cov_inv = beta_tp_cov_inv + beta_tp_cov_prior_inv;
  beta_tp_cov = arma::inv_sympd(beta_tp_cov_inv);
  
  // marginal posterior mean
  C = linpred - (X_sp * beta_tp);
  C_valid = C.rows(valid_obs);
  Z = (Y-0.5)/W - C_valid;
  beta_tp_mean = beta_tp_cov * (X_sp_valid.t() * (W % Z));
  
  // Sampling sp_it
  beta_tp = arma::mvnrnd( beta_tp_mean , beta_tp_cov );
  
  // Recalculate linpred with the new values of sp_it
  linpred = X_sp * beta_tp + C;
  
  // Redefine gamma_ijt with the new values of linpred
  for( k=0; k<K_net; k++ ) {
    // Rcpp::Rcout << "k=" << k << std::endl;
    gamma_ijt = gamma_ijtk(k);
    
    if( !directed ){
      aux_mat_1 = linpred.rows(k*T_net*V_net*(V_net-1)/2, (k+1)*T_net*V_net*(V_net-1)/2-1);
      aux_mat_1.reshape(V_net*(V_net-1)/2,T_net);
      for( i=1; i<V_net; i++ ) {
        aux_mat_2 = aux_mat_1.rows((i-1)*V_net-((i-1)*i)/2,i*V_net-(i*(i+1))/2-1);
        gamma_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1) = aux_mat_2;
      }
    } else {
      aux_mat_1 = linpred.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1);
      aux_mat_1.reshape(V_net*(V_net-1),T_net);
      for( i=0; i<V_net; i++ ) {
        aux_mat_2 = aux_mat_1.rows(i*(V_net-1),(i+1)*(V_net-1)-1);
        aux_mat_2.insert_rows(i,arma::zeros<arma::mat>(1,T_net));
        gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_2;
      }
    }
    
    gamma_ijtk(k) = gamma_ijt;
  }
  
  // reshape beta
  beta_tp.reshape(T_net,P_pred);
  
  return Rcpp::List::create( Rcpp::Named("beta_tp") = beta_tp,
                             Rcpp::Named("gamma_ijtk") = gamma_ijtk );
}
