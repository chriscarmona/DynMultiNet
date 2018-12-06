#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "DynMultiNet_shared.h"

// [[Rcpp::interfaces(r, cpp)]]



// [[Rcpp::export]]
arma::mat sample_mu_t_DynMultiNet_bin_cpp( arma::colvec mu_t,
                                           const arma::mat mu_t_cov_prior_inv,
                                           const arma::cube y_ijt,
                                           const arma::cube w_ijt,
                                           const arma::cube s_ijt,
                                           const bool directed=false ) {
  // Auxiliar objects
  unsigned int i=0;
  
  arma::cube aux_cube_1;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::uvec aux_uvec_1;
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  
  arma::colvec Y = arma::zeros<arma::colvec>(1);
  arma::colvec W = arma::zeros<arma::colvec>(1);
  arma::colvec S = arma::zeros<arma::colvec>(1);
  
  arma::colvec C = arma::zeros<arma::colvec>(1);
  arma::colvec Z = arma::zeros<arma::colvec>(1);
  
  aux_mat_1 = arma::ones<arma::mat>(V_net*(V_net-1)/2,1);
  aux_mat_2 = arma::zeros<arma::mat>(T_net,T_net); aux_mat_2.eye();
  arma::mat X = arma::kron( aux_mat_1, aux_mat_2 );
  arma::sp_mat X_sp = arma::sp_mat(X);
  
  arma::sp_mat Omega_sp=arma::speye<arma::sp_mat>(T_net*V_net*(V_net-1)/2,T_net*V_net*(V_net-1)/2);
  
  if( directed ){
    X_sp = repmat(X_sp,2,1);
    Omega_sp=arma::speye<arma::sp_mat>(T_net*V_net*(V_net-1),T_net*V_net*(V_net-1)); // too large
  }
  
  arma::mat mu_t_cov_inv;
  arma::mat mu_t_cov;
  arma::colvec aux_vec_mean = arma::zeros<arma::colvec>(T_net);
  
  if( directed ){
    for( i=0; i<V_net; i++ ) {
      aux_mat_1 = y_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
      aux_mat_1.shed_row(i);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape((V_net-1)*T_net,1);
      Y.insert_rows( Y.n_rows, aux_mat_1);
      
      aux_mat_1 = w_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
      aux_mat_1.shed_row(i);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape((V_net-1)*T_net,1);
      W.insert_rows( W.n_rows, aux_mat_1);
      
      aux_mat_1 = s_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
      aux_mat_1.shed_row(i);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape((V_net-1)*T_net,1);
      S.insert_rows( S.n_rows, aux_mat_1);
    }
  } else {
    for( i=1; i<V_net; i++ ) {
      aux_mat_1 = y_ijt.subcube(i,0,0, i,i-1,T_net-1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape(i*T_net,1);
      Y.insert_rows( Y.n_rows, aux_mat_1);
      
      aux_mat_1 = w_ijt.subcube(i,0,0, i,i-1,T_net-1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape(i*T_net,1);
      W.insert_rows( W.n_rows, aux_mat_1);
      
      aux_mat_1 = s_ijt.subcube(i,0,0, i,i-1,T_net-1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape(i*T_net,1);
      S.insert_rows( S.n_rows, aux_mat_1);
    }
  }
  Y.shed_row(0); W.shed_row(0); S.shed_row(0);
  Omega_sp.diag() = W;
  
  mu_t_cov_inv = X_sp.t() * Omega_sp * X_sp;
  mu_t_cov_inv = mu_t_cov_inv + mu_t_cov_prior_inv;
  
  mu_t_cov = arma::inv_sympd(mu_t_cov_inv);
  
  C = S - (X_sp * mu_t);
  Z = (Y-0.5)/W - C;
  
  aux_vec_mean = X_sp.t() * (Omega_sp * Z);
  
  mu_t = arma::mvnrnd( mu_t_cov*aux_vec_mean , mu_t_cov );
  
  return mu_t;
}

// [[Rcpp::export]]
arma::mat sample_mu_t_DynMultiNet_bin_v2_cpp( arma::colvec mu_t,
                                              const arma::mat mu_t_cov_prior_inv,
                                              const arma::cube y_ijt,
                                              const arma::cube w_ijt,
                                              const arma::cube s_ijt,
                                              const bool directed=false ) {
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  unsigned int aux_int=0;
  
  arma::cube aux_cube_1;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::uvec aux_uvec_1;
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  
  arma::colvec Y;
  arma::colvec W;
  arma::colvec S;
  
  arma::colvec C = arma::zeros<arma::colvec>(1);
  arma::colvec Z = arma::zeros<arma::colvec>(1);
  
  aux_mat_1 = arma::zeros<arma::mat>(T_net,T_net); aux_mat_1.eye();
  if( directed ){
    aux_mat_2 = arma::ones<arma::mat>(V_net*(V_net-1),1);
  } else {
    aux_mat_2 = arma::ones<arma::mat>(V_net*(V_net-1)/2,1);
  }
  arma::mat X = arma::kron( aux_mat_1, aux_mat_2 );
  arma::sp_mat X_sp = arma::sp_mat(X);
  
  arma::mat mu_t_cov_inv;
  arma::mat mu_t_cov;
  arma::colvec aux_vec_mean = arma::zeros<arma::colvec>(T_net);
  
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
      aux_mat_2 = s_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
      aux_mat_2.shed_row(i);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
    S = aux_mat_1;
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
      aux_mat_2 = s_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
    S = aux_mat_1;
  }
  
  // Omega_sp.diag() = W;
  // mu_t_cov_inv = X_sp.t() * Omega_sp * X_sp;
  // more efficient:
  mu_t_cov_inv = arma::zeros<arma::mat>(T_net,T_net);
  if( directed ){
    aux_int = V_net*(V_net-1);
  } else {
    aux_int = V_net*(V_net-1)/2;
  }
  for( t=0; t<T_net; t++ ) { mu_t_cov_inv(t,t) = sum( W.subvec( t*aux_int, (t+1)*aux_int-1  ) ); }
  mu_t_cov_inv = mu_t_cov_inv + mu_t_cov_prior_inv;
  mu_t_cov = arma::inv_sympd(mu_t_cov_inv);
  
  C = S - (X_sp * mu_t);
  
  Z = (Y-0.5)/W - C;
  
  aux_vec_mean = X_sp.t() * (W % Z);
  
  mu_t = arma::mvnrnd( mu_t_cov*aux_vec_mean , mu_t_cov );
  
  return mu_t;
}

// [[Rcpp::export]]
arma::colvec sample_beta_z_layer_DynMultiNet_bin_cpp( arma::colvec beta_t,
                                                   arma::colvec z_t,
                                                   const arma::mat beta_t_cov_prior_inv,
                                                   const arma::cube y_ijt,
                                                   const arma::cube w_ijt,
                                                   const arma::cube s_ijt,
                                                   const bool directed=false ) {
  // Auxiliar objects
  unsigned int i=0;
  
  arma::cube aux_cube_1;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::uvec aux_uvec_1;
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  
  arma::colvec Y = arma::zeros<arma::colvec>(1);
  arma::colvec W = arma::zeros<arma::colvec>(1);
  arma::colvec S = arma::zeros<arma::colvec>(1);
  
  arma::colvec C = arma::zeros<arma::colvec>(1);
  arma::colvec Z = arma::zeros<arma::colvec>(1);
  
  aux_mat_1 = arma::ones<arma::mat>(V_net*(V_net-1)/2,1);
  aux_mat_2 = arma::zeros<arma::mat>(T_net,T_net); aux_mat_2.diag()=z_t;
  arma::mat X = arma::kron( aux_mat_1, aux_mat_2 );
  arma::sp_mat X_sp = arma::sp_mat(X);
  
  arma::sp_mat Omega_sp = arma::speye<arma::sp_mat>(T_net*V_net*(V_net-1)/2,T_net*V_net*(V_net-1)/2);
  
  if( directed ){
    X_sp = repmat(X_sp,2,1);
    Omega_sp=arma::speye<arma::sp_mat>(T_net*V_net*(V_net-1),T_net*V_net*(V_net-1));
  }
  
  
  arma::mat beta_t_cov_inv;
  arma::mat beta_t_cov;
  arma::colvec aux_vec_mean = arma::zeros<arma::colvec>(T_net);
  
  if( directed ){
    for( i=0; i<V_net; i++ ) {
      aux_mat_1 = y_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
      aux_mat_1.shed_row(i);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape((V_net-1)*T_net,1);
      Y.insert_rows( Y.n_rows, aux_mat_1);
      
      aux_mat_1 = w_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
      aux_mat_1.shed_row(i);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape((V_net-1)*T_net,1);
      W.insert_rows( W.n_rows, aux_mat_1);
      
      aux_mat_1 = s_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
      aux_mat_1.shed_row(i);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape((V_net-1)*T_net,1);
      S.insert_rows( S.n_rows, aux_mat_1);
    }
  } else {
    for( i=1; i<V_net; i++ ) {
      aux_mat_1 = y_ijt.subcube(i,0,0, i,i-1,T_net-1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape(i*T_net,1);
      Y.insert_rows( Y.n_rows, aux_mat_1);
      
      aux_mat_1 = w_ijt.subcube(i,0,0, i,i-1,T_net-1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape(i*T_net,1);
      W.insert_rows( W.n_rows, aux_mat_1);
      
      aux_mat_1 = s_ijt.subcube(i,0,0, i,i-1,T_net-1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape(i*T_net,1);
      S.insert_rows( S.n_rows, aux_mat_1);
    }
  }
  Y.shed_row(0); W.shed_row(0); S.shed_row(0);
  Omega_sp.diag() = W;
  
  beta_t_cov_inv = X_sp.t() * Omega_sp * X_sp;
  beta_t_cov_inv = beta_t_cov_inv + beta_t_cov_prior_inv;
  
  beta_t_cov = arma::inv_sympd(beta_t_cov_inv);
  
  C = S - (X_sp * beta_t);
  Z = (Y-0.5)/W - C;
  
  aux_vec_mean = X_sp.t() * (Omega_sp * Z);
  
  beta_t = arma::mvnrnd( beta_t_cov*aux_vec_mean , beta_t_cov );
  
  return beta_t;
}



// [[Rcpp::export]]
arma::cube sample_x_ith_DynMultiNet_bin_cpp( arma::cube x_ith,
                                             const arma::mat x_t_sigma_prior_inv,
                                             const arma::colvec tau_h,
                                             const arma::cube y_ijt,
                                             const arma::cube w_ijt,
                                             const arma::cube s_ijt ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = x_ith.n_slices;
  
  x_ith = reshape(x_ith,V_net,T_net*H_dim,1);
  arma::mat x_ith_mat = x_ith.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec W = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec S = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec Z = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::sp_mat Omega_sp=arma::speye<arma::sp_mat>((V_net-1)*T_net,(V_net-1)*T_net);
  
  arma::mat tau_h_diag(tau_h.n_rows,tau_h.n_rows); tau_h_diag.eye();
  tau_h_diag.diag() = tau_h;
  
  arma::mat x_i_cov_prior = kron( tau_h_diag, x_t_sigma_prior_inv );
  arma::mat x_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat x_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  
  arma::colvec aux_vec_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  for( t=0; t<T_net; t++ ) {
    X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = x_ith_mat.submat( aux_uvec_3, aux_uvec_2+t );
  }
  
  arma::sp_mat X_all_sp = arma::sp_mat(X_all);
  arma::sp_mat X_sp = X_all_sp;
  
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
    
    aux_mat_1 = s_ijt.subcube(i,0,0, i,i,T_net-1);
    aux_mat_2 = s_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
    aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
    aux_mat_1.shed_rows(i,i+1);
    aux_mat_1 = aux_mat_1.t();
    aux_mat_1.reshape((V_net-1)*T_net,1);
    S = aux_mat_1;
    
    Omega_sp.diag() = W;
    
    // X
    X_sp = X_all_sp;
    X_sp.shed_rows(T_net*i,T_net*(i+1)-1);
    
    x_i_cov_inv = X_sp.t() * Omega_sp * X_sp;
    x_i_cov_inv = x_i_cov_inv + x_i_cov_prior ;
    x_i_cov = arma::inv_sympd(x_i_cov_inv);
    
    C = S - X_sp * trans(x_ith_mat.row(i));
    Z = (Y-0.5)/W - C;
    
    aux_vec_mean = X_sp.t() * (W % Z);
    
    x_ith_mat.row(i) = trans(arma::mvnrnd( x_i_cov*aux_vec_mean , x_i_cov ));
  }
  
  // get x_ith from x_ith_mat
  x_ith.slice(0)=x_ith_mat;
  x_ith = reshape(x_ith,V_net,T_net,H_dim);
  
  return x_ith;
}

// [[Rcpp::export]]
arma::cube sample_x_ith_shared_DynMultiNet_bin_cpp( arma::cube x_ith_shared,
                                                    const arma::mat x_t_sigma_prior_inv,
                                                    const arma::colvec tau_h,
                                                    const arma::field<arma::cube> y_ijtk,
                                                    const arma::field<arma::cube> w_ijtk,
                                                    const arma::field<arma::cube> s_ijtk ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int k=0;
  unsigned int t=0;
  
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  arma::cube y_ijt = y_ijtk(0);
  arma::cube w_ijt = w_ijtk(0);
  arma::cube s_ijt = s_ijtk(0);
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = x_ith_shared.n_slices;
  unsigned int K_net = y_ijtk.n_rows;
  
  x_ith_shared = reshape(x_ith_shared,V_net,T_net*H_dim,1);
  arma::mat x_ith_shared_mat = x_ith_shared.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec W = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec S = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec Z = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  arma::sp_mat Omega_sp=arma::speye<arma::sp_mat>((V_net-1)*T_net*K_net,(V_net-1)*T_net*K_net);
  
  arma::mat tau_h_diag(tau_h.n_rows,tau_h.n_rows); tau_h_diag.eye();
  tau_h_diag.diag() = tau_h;
  
  arma::mat x_i_cov_prior = kron( tau_h_diag, x_t_sigma_prior_inv );
  arma::mat x_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat x_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  
  arma::colvec aux_vec_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  for( t=0; t<T_net; t++ ) {
    X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = x_ith_shared_mat.submat( aux_uvec_3, aux_uvec_2+t );
  }
  
  arma::sp_mat X_all_sp = arma::sp_mat(X_all);
  arma::sp_mat X_sp = X_all_sp;
  
  for( i=0; i<V_net; i++ ) {
    for( k=0; k<K_net; k++ ) {
      
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
      
      s_ijt = s_ijtk(k);
      aux_mat_1 = s_ijt.subcube(i,0,0, i,i,T_net-1);
      aux_mat_2 = s_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
      aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
      aux_mat_1.shed_rows(i,i+1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.reshape((V_net-1)*T_net,1);
      S.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
    }
    Omega_sp.diag() = W;
    
    // X
    X_sp = X_all_sp;
    X_sp.shed_rows(T_net*i,T_net*(i+1)-1);
    X_sp = repmat(X_sp,K_net,1);
    
    x_i_cov_inv = X_sp.t() * Omega_sp * X_sp;
    x_i_cov_inv = x_i_cov_inv + x_i_cov_prior ;
    x_i_cov = arma::inv_sympd(x_i_cov_inv);
    
    C = S - X_sp * trans(x_ith_shared_mat.row(i));
    Z = (Y-0.5)/W - C;
    
    aux_vec_mean = X_sp.t() * (W % Z);
    
    x_ith_shared_mat.row(i) = trans(arma::mvnrnd( x_i_cov*aux_vec_mean , x_i_cov ));
  }
  
  // get x_ith_shared from x_ith_shared_mat
  x_ith_shared.slice(0)=x_ith_shared_mat;
  x_ith_shared = reshape(x_ith_shared,V_net,T_net,H_dim);
  
  return x_ith_shared;
}


// [[Rcpp::export]]
Rcpp::List sample_x_ith_DynMultiNet_bin_dir_cpp( arma::cube x_ith_send,
                                                 arma::cube x_ith_receive,
                                                 const arma::mat x_t_sigma_prior_inv,
                                                 const arma::colvec tau_h_send,
                                                 const arma::colvec tau_h_receive,
                                                 const arma::cube y_ijt,
                                                 const arma::cube w_ijt,
                                                 const arma::cube s_ijt ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  
  unsigned int dir=0;
  
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = x_ith_send.n_slices;
  
  x_ith_send = reshape(x_ith_send,V_net,T_net*H_dim,1);
  arma::mat x_ith_send_mat = x_ith_send.slice(0);
  
  x_ith_receive = reshape(x_ith_receive,V_net,T_net*H_dim,1);
  arma::mat x_ith_receive_mat = x_ith_receive.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec W = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec S = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec Z = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::sp_mat Omega_sp=arma::speye<arma::sp_mat>((V_net-1)*T_net,(V_net-1)*T_net);
  arma::mat tau_h_diag(tau_h_send.n_rows,tau_h_send.n_rows); tau_h_diag.eye();
  arma::mat x_i_cov_prior = kron( tau_h_diag, x_t_sigma_prior_inv );
  arma::mat x_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat x_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::colvec aux_vec_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  // Model matrix for all senders, made out of receiver coordinates
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  for( t=0; t<T_net; t++ ) {
    X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = x_ith_receive_mat.submat( aux_uvec_3, aux_uvec_2+t );
  }
  arma::sp_mat X_all_sp_send = arma::sp_mat(X_all);
  
  // Model matrix for all receivers, made out of sender coordinates
  X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  for( t=0; t<T_net; t++ ) {
    X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = x_ith_send_mat.submat( aux_uvec_3, aux_uvec_2+t );
  }
  arma::sp_mat X_all_sp_receive = arma::sp_mat(X_all);
  
  // Model matrix that will be changed (V_net*2 times) for each node in the cycle, in their two positions: sender/receiver
  arma::sp_mat X_sp = X_all_sp_send;
  
  for( dir=0; dir<2; dir++ ) {
    
    Y = arma::zeros<arma::colvec>((V_net-1)*T_net);
    W = arma::zeros<arma::colvec>((V_net-1)*T_net);
    S = arma::zeros<arma::colvec>((V_net-1)*T_net);
    C = arma::zeros<arma::colvec>((V_net-1)*T_net);
    Z = arma::zeros<arma::colvec>((V_net-1)*T_net);
    
    for( i=0; i<V_net; i++ ) {
      if( dir==0 ) {
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
        
        aux_mat_1 = s_ijt.subcube(i,0,0, i,i,T_net-1);
        aux_mat_2 = s_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
        aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
        aux_mat_1.shed_rows(i,i+1);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        S = aux_mat_1;
        
        // X
        X_sp = X_all_sp_send;
        X_sp.shed_rows(T_net*i,T_net*(i+1)-1);
        
        tau_h_diag.diag() = tau_h_send;
        
      } else if( dir==1 ) {
        aux_mat_1 = y_ijt.subcube(0,i,0, i,i,T_net-1);
        aux_mat_2 = y_ijt.subcube(i,i,0, i,V_net-1,T_net-1);
        aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
        aux_mat_1.shed_rows(i,i+1);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        Y = aux_mat_1;
        
        aux_mat_1 = w_ijt.subcube(0,i,0, i,i,T_net-1);
        aux_mat_2 = w_ijt.subcube(i,i,0, i,V_net-1,T_net-1);
        aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
        aux_mat_1.shed_rows(i,i+1);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        W = aux_mat_1;
        
        aux_mat_1 = s_ijt.subcube(0,i,0, i,i,T_net-1);
        aux_mat_2 = s_ijt.subcube(i,i,0, i,V_net-1,T_net-1);
        aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
        aux_mat_1.shed_rows(i,i+1);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        S = aux_mat_1;
        
        // X
        X_sp = X_all_sp_receive;
        X_sp.shed_rows(T_net*i,T_net*(i+1)-1);
        
        tau_h_diag.diag() = tau_h_receive;
        
      } else {
        throw std::range_error("direction not supported");
      }
      
      Omega_sp.diag() = W;
      
      x_i_cov_prior = kron( tau_h_diag, x_t_sigma_prior_inv );
      
      x_i_cov_inv = X_sp.t() * Omega_sp * X_sp;
      x_i_cov_inv = x_i_cov_inv + x_i_cov_prior ;
      x_i_cov = arma::inv_sympd(x_i_cov_inv);
      
      if( dir==0 ) {
        C = S - X_sp * trans(x_ith_send_mat.row(i));
        Z = (Y-0.5)/W - C;
        aux_vec_mean = X_sp.t() * (W % Z);
        x_ith_send_mat.row(i) = trans(arma::mvnrnd( x_i_cov*aux_vec_mean , x_i_cov ));
      } else if( dir==1 ) {
        C = S - X_sp * trans(x_ith_receive_mat.row(i));
        Z = (Y-0.5)/W - C;
        aux_vec_mean = X_sp.t() * (W % Z);
        x_ith_receive_mat.row(i) = trans(arma::mvnrnd( x_i_cov*aux_vec_mean , x_i_cov ));
      } else {
        throw std::range_error("direction not supported");
      }
    }
  }
  
  // get x_ith_shared from x_ith_shared_mat
  x_ith_send.slice(0)=x_ith_send_mat;
  x_ith_send = reshape(x_ith_send,V_net,T_net,H_dim);
  
  x_ith_receive.slice(0)=x_ith_receive_mat;
  x_ith_receive = reshape(x_ith_receive,V_net,T_net,H_dim);
  
  return Rcpp::List::create( Rcpp::Named("send") = x_ith_send,
                             Rcpp::Named("receive") = x_ith_receive );
}


// [[Rcpp::export]]
Rcpp::List sample_x_ith_shared_DynMultiNet_bin_dir_cpp( arma::cube x_ith_shared_send,
                                                        arma::cube x_ith_shared_receive,
                                                        const arma::mat x_t_sigma_prior_inv,
                                                        const arma::colvec tau_h_shared_send,
                                                        const arma::colvec tau_h_shared_receive,
                                                        const arma::field<arma::cube> y_ijtk,
                                                        const arma::field<arma::cube> w_ijtk,
                                                        const arma::field<arma::cube> s_ijtk ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int k=0;
  unsigned int t=0;
  
  unsigned int dir=0;
  
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  arma::cube y_ijt = y_ijtk(0);
  arma::cube w_ijt = w_ijtk(0);
  arma::cube s_ijt = s_ijtk(0);
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = x_ith_shared_send.n_slices;
  unsigned int K_net = y_ijtk.n_rows;
  
  x_ith_shared_send = reshape(x_ith_shared_send,V_net,T_net*H_dim,1);
  arma::mat x_ith_shared_send_mat = x_ith_shared_send.slice(0);
  
  x_ith_shared_receive = reshape(x_ith_shared_receive,V_net,T_net*H_dim,1);
  arma::mat x_ith_shared_receive_mat = x_ith_shared_receive.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec W = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec S = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec Z = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  arma::sp_mat Omega_sp=arma::speye<arma::sp_mat>((V_net-1)*T_net*K_net,(V_net-1)*T_net*K_net);
  arma::mat tau_h_diag(tau_h_shared_send.n_rows,tau_h_shared_send.n_rows); tau_h_diag.eye();
  arma::mat x_i_cov_prior = kron( tau_h_diag, x_t_sigma_prior_inv );
  arma::mat x_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat x_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::colvec aux_vec_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  // Model matrix for all senders, made out of receiver coordinates
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  for( t=0; t<T_net; t++ ) {
    X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = x_ith_shared_receive_mat.submat( aux_uvec_3, aux_uvec_2+t );
  }
  arma::sp_mat X_all_sp_send = arma::sp_mat(X_all);
  
  // Model matrix for all receivers, made out of sender coordinates
  X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  for( t=0; t<T_net; t++ ) {
    X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = x_ith_shared_send_mat.submat( aux_uvec_3, aux_uvec_2+t );
  }
  arma::sp_mat X_all_sp_receive = arma::sp_mat(X_all);
  
  // Model matrix that will be changed (V_net*2 times) for each node in the cycle, in their two positions: sender/receiver
  arma::sp_mat X_sp = X_all_sp_send;
  
  for( dir=0; dir<2; dir++ ) {
    
    Y = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
    W = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
    S = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
    C = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
    Z = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
    
    for( i=0; i<V_net; i++ ) {
      if( dir==0 ) {
        for( k=0; k<K_net; k++ ) {
          
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
          
          s_ijt = s_ijtk(k);
          aux_mat_1 = s_ijt.subcube(i,0,0, i,i,T_net-1);
          aux_mat_2 = s_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
          aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
          aux_mat_1.shed_rows(i,i+1);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          S.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
        }
        // X
        X_sp = X_all_sp_send;
        X_sp.shed_rows(T_net*i,T_net*(i+1)-1);
        X_sp = repmat(X_sp,K_net,1);
        
        tau_h_diag.diag() = tau_h_shared_send;
        
      } else if( dir==1 ) {
        for( k=0; k<K_net; k++ ) {
          y_ijt = y_ijtk(k);
          aux_mat_1 = y_ijt.subcube(0,i,0, i,i,T_net-1);
          aux_mat_2 = y_ijt.subcube(i,i,0, i,V_net-1,T_net-1);
          aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
          aux_mat_1.shed_rows(i,i+1);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          Y.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
          
          w_ijt = w_ijtk(k);
          aux_mat_1 = w_ijt.subcube(0,i,0, i,i,T_net-1);
          aux_mat_2 = w_ijt.subcube(i,i,0, i,V_net-1,T_net-1);
          aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
          aux_mat_1.shed_rows(i,i+1);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          W.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
          
          s_ijt = s_ijtk(k);
          aux_mat_1 = s_ijt.subcube(0,i,0, i,i,T_net-1);
          aux_mat_2 = s_ijt.subcube(i,i,0, i,V_net-1,T_net-1);
          aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
          aux_mat_1.shed_rows(i,i+1);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          S.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
        }
        
        // X
        X_sp = X_all_sp_receive;
        X_sp.shed_rows(T_net*i,T_net*(i+1)-1);
        X_sp = repmat(X_sp,K_net,1);
        
        tau_h_diag.diag() = tau_h_shared_receive;
        
      } else {
        throw std::range_error("direction not supported");
      }
      
      Omega_sp.diag() = W;
      
      x_i_cov_prior = kron( tau_h_diag, x_t_sigma_prior_inv );
      
      x_i_cov_inv = X_sp.t() * Omega_sp * X_sp;
      x_i_cov_inv = x_i_cov_inv + x_i_cov_prior ;
      x_i_cov = arma::inv_sympd(x_i_cov_inv);
      
      if( dir==0 ) {
        C = S - X_sp * trans(x_ith_shared_send_mat.row(i));
        Z = (Y-0.5)/W - C;
        aux_vec_mean = X_sp.t() * (W % Z);
        x_ith_shared_send_mat.row(i) = trans(arma::mvnrnd( x_i_cov*aux_vec_mean , x_i_cov ));
      } else if( dir==1 ) {
        C = S - X_sp * trans(x_ith_shared_receive_mat.row(i));
        Z = (Y-0.5)/W - C;
        aux_vec_mean = X_sp.t() * (W % Z);
        x_ith_shared_receive_mat.row(i) = trans(arma::mvnrnd( x_i_cov*aux_vec_mean , x_i_cov ));
      } else {
        throw std::range_error("direction not supported");
      }
    }
  }
  
  // get x_ith_shared from x_ith_shared_mat
  x_ith_shared_send.slice(0)=x_ith_shared_send_mat;
  x_ith_shared_send = reshape(x_ith_shared_send,V_net,T_net,H_dim);
  
  x_ith_shared_receive.slice(0)=x_ith_shared_receive_mat;
  x_ith_shared_receive = reshape(x_ith_shared_receive,V_net,T_net,H_dim);
  
  return Rcpp::List::create( Rcpp::Named("send") = x_ith_shared_send,
                             Rcpp::Named("receive") = x_ith_shared_receive );
}
