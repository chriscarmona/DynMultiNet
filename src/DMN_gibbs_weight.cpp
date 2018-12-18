#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "DynMultiNet_shared.h"

// [[Rcpp::interfaces(r, cpp)]]


// [[Rcpp::export]]
Rcpp::List sample_baseline_tk_weight_cpp( arma::colvec theta_t,
                                          const arma::mat theta_t_cov_prior_inv,
                                          const arma::cube y_ijt,
                                          arma::cube mu_ijt,
                                          const double sigma_k,
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
  
  // column matrix for the continuous response
  arma::colvec Y = arma::zeros<arma::colvec>(1);
  
  // column matrix with the linear predictor
  arma::colvec linpred = arma::zeros<arma::colvec>(1);
  
  // matrix with the part of the linear predictor that does not depend on theta
  arma::colvec C = arma::zeros<arma::colvec>(1);
  
  // Design matrix X that makes
  // linpred = X * theta_t + C
  aux_mat_1 = arma::zeros<arma::mat>(T_net,T_net); aux_mat_1.eye();
  if( directed ){
    aux_mat_2 = arma::ones<arma::mat>(V_net*(V_net-1),1);
  } else {
    aux_mat_2 = arma::ones<arma::mat>(V_net*(V_net-1)/2,1);
  }
  arma::mat X = arma::kron( aux_mat_1, aux_mat_2 );
  arma::sp_mat X_sp = arma::sp_mat(X);
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  arma::colvec Y_valid = Y;
  arma::sp_mat X_sp_valid = X_sp;
  arma::colvec C_valid = C;
  
  // Marginal posterior mean vector for theta_t
  arma::colvec theta_t_mean = arma::zeros<arma::colvec>(T_net);
  
  // Marginal posterior covariance matrix for theta_t
  arma::mat theta_t_cov = arma::zeros<arma::mat>(T_net,T_net);
  arma::mat theta_t_cov_inv = arma::zeros<arma::mat>(T_net,T_net);
  
  
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
      aux_mat_2 = mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
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
      aux_mat_2 = mu_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1)/2,1);
    linpred = aux_mat_1;
  }
  
  // Constant term for theta in the linear predictor
  C = linpred - ( X_sp * theta_t );
  
  // identifies valid obs
  Y.elem( find_nonfinite(Y) ).zeros(); // change NAs to zero
  valid_obs = find(Y); // find elements different than zero
  Y_valid = Y.rows(valid_obs);
  X_sp_valid = arma::sp_mat(X.rows(valid_obs));
  C_valid = C.rows(valid_obs);
  
  // Marginal Posterior
  // Covariance
  theta_t_cov_inv = (1/pow(sigma_k,2)) * X_sp_valid.t() * X_sp_valid + theta_t_cov_prior_inv ;
  theta_t_cov = arma::inv_sympd(theta_t_cov_inv);
  // Mean
  theta_t_mean = theta_t_cov * ( (1/pow(sigma_k,2)) * X_sp_valid.t() * (Y_valid-C_valid) );
  
  // sampling theta_t
  theta_t = arma::mvnrnd( theta_t_mean , theta_t_cov );
  
  // return theta_t;
  
  // Recalculate linpred with the new values of mu
  linpred = X_sp * theta_t + C;
  // Redefine mu_ijt with the new values of linpred
  if( directed ){
    aux_mat_1 = linpred;
    aux_mat_1.reshape(V_net*(V_net-1),T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = aux_mat_1.rows(i*(V_net-1),(i+1)*(V_net-1)-1);
      aux_mat_2.insert_rows(i,arma::zeros<arma::mat>(1,T_net));
      mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_2;
    }
  } else {
    aux_mat_1 = linpred;
    aux_mat_1.reshape(V_net*(V_net-1)/2,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = aux_mat_1.rows((i-1)*V_net-((i-1)*i)/2,i*V_net-(i*(i+1))/2-1);
      mu_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1) = aux_mat_2;
    }
  }
  
  return Rcpp::List::create( Rcpp::Named("theta_t") = theta_t,
                             Rcpp::Named("mu_ijt") = mu_ijt );
}


// To be implemented
// [[Rcpp::export]]
Rcpp::List sample_coord_ith_weight_cpp( arma::cube uv_ith_shared,
                                        const arma::mat uv_t_sigma_prior_inv,
                                        const arma::colvec tau_h,
                                        const arma::cube y_ijt,
                                        arma::cube mu_ijt,
                                        const double sigma_k ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = uv_ith_shared.n_slices;
  
  uv_ith_shared = reshape(uv_ith_shared,V_net,T_net*H_dim,1);
  arma::mat uv_ith_shared_mat = uv_ith_shared.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec linpred = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::mat tau_h_diag(tau_h.n_rows,tau_h.n_rows); tau_h_diag.eye();
  tau_h_diag.diag() = tau_h;
  
  arma::mat uv_i_cov_prior = kron( tau_h_diag, uv_t_sigma_prior_inv );
  arma::mat uv_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat uv_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  
  arma::colvec uv_i_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  // Performing sampling in uv_ith_shared...
  // The rows in uv_ith_shared_mat will act as the targeted "coefficients"
  // The columns in X_all will act as the "covariates"
  // The order of the response mu_ijth will depend on the order of rows in X_all
  
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  arma::mat X = X_all;
  arma::sp_mat X_sp = arma::sp_mat(X);
  
  for( i=0; i<V_net; i++ ) {
    aux_mat_1 = y_ijt.subcube(i,0,0, i,i,T_net-1);
    aux_mat_2 = y_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
    aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
    aux_mat_1.shed_rows(i,i+1);
    aux_mat_1 = aux_mat_1.t();
    aux_mat_1.reshape((V_net-1)*T_net,1);
    Y = aux_mat_1;
    
    aux_mat_1 = mu_ijt.subcube(i,0,0, i,i,T_net-1);
    aux_mat_2 = mu_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
    aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
    aux_mat_1.shed_rows(i,i+1);
    aux_mat_1 = aux_mat_1.t();
    aux_mat_1.reshape((V_net-1)*T_net,1);
    linpred = aux_mat_1;
    
    // X
    X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
    for( t=0; t<T_net; t++ ) {
      X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = uv_ith_shared_mat.submat( aux_uvec_3, aux_uvec_2+t );
    }
    X = X_all;
    X.shed_rows(T_net*i,T_net*(i+1)-1);
    X_sp = arma::sp_mat(X);
    
    uv_i_cov_inv = X_sp.t() * X_sp;
    uv_i_cov_inv = uv_i_cov_inv + uv_i_cov_prior ;
    uv_i_cov = arma::inv_sympd(uv_i_cov_inv);
    
    C = linpred - X_sp * trans(uv_ith_shared_mat.row(i));
    
    // Sampling uv_ith_shared_mat
    uv_ith_shared_mat.row(i) = trans(arma::mvnrnd( uv_i_mean , uv_i_cov ));
    
    // Recalculate linpred with the new values of uv_ith_shared_mat
    linpred = X_sp * trans(uv_ith_shared_mat.row(i)) + C;
    // Redefine mu_ijt with the new values of linpred
    aux_mat_1 = linpred;
    aux_mat_1.reshape(T_net,(V_net-1));
    aux_mat_1 = aux_mat_1.t();
    aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
    mu_ijt.subcube(i,0,0, i,i,T_net-1) = aux_mat_1.rows(0,i);
    mu_ijt.subcube(i,i,0, V_net-1,i,T_net-1) = aux_mat_1.rows(i,V_net-1);
  }
  
  // get uv_ith_shared from uv_ith_shared_mat
  uv_ith_shared.slice(0)=uv_ith_shared_mat;
  uv_ith_shared = reshape(uv_ith_shared,V_net,T_net,H_dim);
  
  return Rcpp::List::create( Rcpp::Named("uv_ith_shared") = uv_ith_shared,
                             Rcpp::Named("mu_ijt") = mu_ijt );
}


// To be implemented
// [[Rcpp::export]]
Rcpp::List sample_coord_ith_shared_weight_cpp( arma::cube uv_ith_shared,
                                               const arma::mat uv_t_sigma_prior_inv,
                                               const arma::colvec tau_h,
                                               const arma::field<arma::cube> y_ijtk,
                                               arma::field<arma::cube> mu_ijtk,
                                               const double sigma_k ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int k=0;
  unsigned int t=0;
  
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  arma::cube y_ijt = y_ijtk(0);
  arma::cube mu_ijt = mu_ijtk(0);
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = uv_ith_shared.n_slices;
  unsigned int K_net = y_ijtk.n_rows;
  
  uv_ith_shared = reshape(uv_ith_shared,V_net,T_net*H_dim,1);
  arma::mat uv_ith_shared_mat = uv_ith_shared.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec linpred = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  arma::mat tau_h_diag(tau_h.n_rows,tau_h.n_rows); tau_h_diag.eye();
  tau_h_diag.diag() = tau_h;
  
  arma::mat uv_i_cov_prior = kron( tau_h_diag, uv_t_sigma_prior_inv );
  arma::mat uv_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat uv_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  
  arma::colvec uv_i_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  arma::sp_mat X_sp = arma::sp_mat(X_all);
  
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
      
      mu_ijt = mu_ijtk(k);
      aux_mat_1 = mu_ijt.subcube(i,0,0, i,i,T_net-1);
      aux_mat_2 = mu_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
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
      X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = uv_ith_shared_mat.submat( aux_uvec_3, aux_uvec_2+t );
    }
    X_sp = arma::sp_mat(X_all);
    X_sp.shed_rows(T_net*i,T_net*(i+1)-1);
    X_sp = repmat(X_sp,K_net,1);
    
    uv_i_cov_inv = X_sp.t() * X_sp;
    uv_i_cov_inv = uv_i_cov_inv + uv_i_cov_prior ;
    uv_i_cov = arma::inv_sympd(uv_i_cov_inv);
    
    C = linpred - X_sp * trans(uv_ith_shared_mat.row(i));
    
    // Sampling uv_ith_shared_mat
    uv_ith_shared_mat.row(i) = trans(arma::mvnrnd( uv_i_mean , uv_i_cov ));
    
    // Rcpp::Rcout << i << std::endl;
    
    // Recalculate linpred with the new values of uv_ith_shared_mat
    linpred = X_sp * trans(uv_ith_shared_mat.row(i)) + C;
    // Redefine mu_ijtk with the new values of linpred
    for( k=0; k<K_net; k++ ) {
      aux_mat_1 = linpred.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1);
      aux_mat_1.reshape(T_net,V_net-1);
      aux_mat_1 = aux_mat_1.t();
      aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
      
      mu_ijt = mu_ijtk(k);
      mu_ijt.subcube(i,0,0, i,i,T_net-1) = aux_mat_1.rows(0,i);
      mu_ijt.subcube(i,i,0, V_net-1,i,T_net-1) = aux_mat_1.rows(i,V_net-1);
      mu_ijtk(k)=mu_ijt;
    }
    
  }
  
  // get uv_ith_shared from uv_ith_shared_mat
  uv_ith_shared.slice(0)=uv_ith_shared_mat;
  uv_ith_shared = reshape(uv_ith_shared,V_net,T_net,H_dim);
  
  return Rcpp::List::create( Rcpp::Named("uv_ith_shared") = uv_ith_shared,
                             Rcpp::Named("mu_ijtk") = mu_ijtk );
}


// To be checked
// [[Rcpp::export]]
Rcpp::List sample_coord_ith_weight_dir_cpp( arma::cube u_ith,
                                            arma::cube v_ith,
                                            const arma::mat uv_t_sigma_prior_inv,
                                            const arma::colvec tau_h_send,
                                            const arma::colvec tau_h_receive,
                                            const arma::cube y_ijt,
                                            arma::cube mu_ijt,
                                            const double sigma_k ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  
  unsigned int dir=0;
  
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = u_ith.n_slices;
  
  u_ith = reshape(u_ith,V_net,T_net*H_dim,1);
  arma::mat u_ith_mat = u_ith.slice(0);
  
  v_ith = reshape(v_ith,V_net,T_net*H_dim,1);
  arma::mat v_ith_mat = v_ith.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net);
  arma::colvec linpred = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::mat tau_h_diag(tau_h_send.n_rows,tau_h_send.n_rows); tau_h_diag.eye();
  arma::mat uv_i_cov_prior = kron( tau_h_diag, uv_t_sigma_prior_inv );
  arma::mat uv_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat uv_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::colvec uv_i_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  // Model matrix for all latent coordinates
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  arma::mat X = X_all;
  
  // Model matrix that will be changed (V_net*2 times) for each node in the cycle, in their two positions: sender/receiver
  arma::sp_mat X_sp = arma::sp_mat(X_all);
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  arma::colvec Y_valid = Y;
  arma::sp_mat X_sp_valid = X_sp;
  arma::colvec C_valid = C;
  
  for( dir=0; dir<2; dir++ ) {
    
    for( i=0; i<V_net; i++ ) {
      if( dir==0 ) {
        aux_mat_1 = y_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
        aux_mat_1.shed_row(i);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        Y = aux_mat_1;
        
        aux_mat_1 = mu_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
        aux_mat_1.shed_row(i);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        linpred = aux_mat_1;
        
        // X
        // Model matrix for all senders, made out of receiver coordinates
        if(i==0){ // The receiver coordinates doesn't change while sampling the senders
          X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
          for( t=0; t<T_net; t++ ) {
            X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = v_ith_mat.submat( aux_uvec_3, aux_uvec_2+t );
          }
          tau_h_diag.diag() = tau_h_send;
        }
        X = X_all;
        X.shed_rows(T_net*i,T_net*(i+1)-1);
        X_sp = arma::sp_mat(X);
        
      } else if( dir==1 ) {
        aux_mat_1 = y_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_1.shed_row(i);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        Y = aux_mat_1;
        
        aux_mat_1 = mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_1.shed_row(i);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.reshape((V_net-1)*T_net,1);
        linpred = aux_mat_1;
        
        // X
        // Model matrix for all receivers, made out of sender coordinates
        if(i==0){ // The sender coordinates doesn't change while sampling the receivers
          X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
          for( t=0; t<T_net; t++ ) {
            X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = u_ith_mat.submat( aux_uvec_3, aux_uvec_2+t );
          }
          tau_h_diag.diag() = tau_h_receive;
        }
        X = X_all;
        X.shed_rows(T_net*i,T_net*(i+1)-1);
        X_sp = arma::sp_mat(X);
        
      } else {
        throw std::range_error("direction not supported");
      }
      
      // identifies valid obs
      Y.elem( find_nonfinite(Y) ).zeros(); // change NAs to zero
      valid_obs = find(Y); // find elements different than zero
      Y_valid = Y.rows(valid_obs);
      X_sp_valid = arma::sp_mat(X.rows(valid_obs));
      
      // The prior for the latent coordinates is uv_t_sigma_prior_inv repeated for each latent dimension
      uv_i_cov_prior = kron( tau_h_diag, uv_t_sigma_prior_inv );
      
      // Marginal Posterior
      // Covariance
      uv_i_cov_inv = (1/pow(sigma_k,2)) * X_sp_valid.t() * X_sp_valid + uv_i_cov_prior ;
      uv_i_cov = arma::inv_sympd(uv_i_cov_inv);
      
      if( dir==0 ) {
        C = linpred - X_sp * trans(u_ith_mat.row(i));
        C_valid = C.rows(valid_obs);
        
        // Marginal Posterior
        // Mean
        uv_i_mean = uv_i_cov * ( (1/pow(sigma_k,2)) * X_sp_valid.t() * (Y_valid-C_valid) );
        
        // Sampling u_ith_mat
        u_ith_mat.row(i) = trans(arma::mvnrnd( uv_i_mean , uv_i_cov ));
        
        // Recalculate linpred with the new values of u_ith_shared_mat
        linpred = X_sp * trans(u_ith_mat.row(i)) + C;
        // Redefine mu_ijt with the new values of linpred
        aux_mat_1 = linpred;
        aux_mat_1.reshape(T_net,V_net-1);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
        mu_ijt.subcube(i,0,0, i,V_net-1,T_net-1) = aux_mat_1;
        
      } else if( dir==1 ) {
        C = linpred - X_sp * trans(v_ith_mat.row(i));
        C_valid = C.rows(valid_obs);
        
        // Marginal Posterior
        // Mean
        uv_i_mean = uv_i_cov * ( (1/pow(sigma_k,2)) * X_sp_valid.t() * (Y_valid-C_valid) );
        
        // Sampling v_ith_mat
        v_ith_mat.row(i) = trans(arma::mvnrnd( uv_i_mean , uv_i_cov ));
        
        // Recalculate linpred with the new values of u_ith_shared_mat
        linpred = X_sp * trans(v_ith_mat.row(i)) + C;
        // Redefine mu_ijt with the new values of linpred
        aux_mat_1 = linpred;
        aux_mat_1.reshape(T_net,V_net-1);
        aux_mat_1 = aux_mat_1.t();
        aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
        mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_1;
      } else {
        throw std::range_error("direction not supported");
      }
    }
  }
  
  // get uv_ith_shared from uv_ith_shared_mat
  u_ith.slice(0)=u_ith_mat;
  u_ith = reshape(u_ith,V_net,T_net,H_dim);
  
  v_ith.slice(0)=v_ith_mat;
  v_ith = reshape(v_ith,V_net,T_net,H_dim);
  
  return Rcpp::List::create( Rcpp::Named("u_ith") = u_ith,
                             Rcpp::Named("v_ith") = v_ith,
                             Rcpp::Named("mu_ijt") = mu_ijt );
}


// To be checked
// [[Rcpp::export]]
Rcpp::List sample_coord_ith_shared_weight_dir_cpp( arma::cube u_ith_shared,
                                                   arma::cube v_ith_shared,
                                                   const arma::mat uv_t_sigma_prior_inv,
                                                   const arma::colvec tau_h_shared_send,
                                                   const arma::colvec tau_h_shared_receive,
                                                   const arma::field<arma::cube> y_ijtk,
                                                   arma::field<arma::cube> mu_ijtk,
                                                   const arma::colvec sigma_k ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int k=0;
  unsigned int t=0;
  
  unsigned int dir=0;
  
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  arma::cube y_ijt = y_ijtk(0);
  arma::cube mu_ijt = mu_ijtk(0);
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = u_ith_shared.n_slices;
  unsigned int K_net = y_ijtk.n_rows;
  
  u_ith_shared = reshape(u_ith_shared,V_net,T_net*H_dim,1);
  arma::mat u_ith_shared_mat = u_ith_shared.slice(0);
  
  v_ith_shared = reshape(v_ith_shared,V_net,T_net*H_dim,1);
  arma::mat v_ith_shared_mat = v_ith_shared.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  arma::colvec linpred = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  
  // Variance associated with each observation in Y
  arma::colvec sigma_Y_inv = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  for( k=0; k<K_net; k++ ) { sigma_Y_inv.subvec( k*(V_net-1)*T_net+1,(k+1)*(V_net-1)*T_net).fill(1/sigma_k(k)); }
  
  arma::mat tau_h_diag(tau_h_shared_send.n_rows,tau_h_shared_send.n_rows); tau_h_diag.eye();
  arma::mat uv_i_cov_prior = kron( tau_h_diag, uv_t_sigma_prior_inv );
  arma::mat uv_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat uv_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::colvec uv_i_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  // Model matrix for all senders, made out of receiver coordinates
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  arma::mat X = X_all;
  
  // Model matrix that will be changed (V_net*2 times) for each node in the cycle, in their two positions: sender/receiver
  arma::sp_mat X_sp = arma::sp_mat(X_all);
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  arma::colvec Y_valid = Y;
  arma::sp_mat X_sp_valid = X_sp;
  arma::colvec C_valid = C;
  arma::colvec sigma_Y_inv_valid = sigma_Y_inv;
  arma::sp_mat sigma_Y_inv_valid_mat = arma::speye<arma::sp_mat>(1,1);
  
  for( dir=0; dir<2; dir++ ) {
    // Rcpp::Rcout << "dir=" << dir << std::endl;
    
    Y = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
    linpred = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
    C = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
    
    X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
    X_sp = arma::sp_mat(X_all);
    
    for( i=0; i<V_net; i++ ) {
      // Rcpp::Rcout << " i=" << i << std::endl;
      
      if( dir==0 ) {
        // Updating Y, W, linpred as vectors, lower triangular adjacency
        for( k=0; k<K_net; k++ ) {
          
          y_ijt = y_ijtk(k);
          aux_mat_1 = y_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
          aux_mat_1.shed_row(i);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          Y.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
          
          mu_ijt = mu_ijtk(k);
          aux_mat_1 = mu_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
          aux_mat_1.shed_row(i);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          linpred.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
        }
        // X_sp
        if(i==0){
          // The matrix of covariates X will only be updated at the beginning of the cycle
          // as it does not change with the newly sampled values
          X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
          for( t=0; t<T_net; t++ ) {
            X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = v_ith_shared_mat.submat( aux_uvec_3, aux_uvec_2+t );
          }
          tau_h_diag.diag() = tau_h_shared_send;
        }
        X = X_all;
        X.shed_rows(T_net*i,T_net*(i+1)-1);
        X = repmat(X,K_net,1);
        X_sp = arma::sp_mat(X);
        
      } else if( dir==1 ) {
        // Updating Y, W, linpred as vectors, upper triangular adjacency
        for( k=0; k<K_net; k++ ) {
          y_ijt = y_ijtk(k);
          aux_mat_1 = y_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
          aux_mat_1.shed_row(i);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          Y.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
          
          mu_ijt = mu_ijtk(k);
          aux_mat_1 = mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
          aux_mat_1.shed_row(i);
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.reshape((V_net-1)*T_net,1);
          linpred.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1) = aux_mat_1;
        }
        // X_sp
        if(i==0){
          // The full matrix of covariates X will only be updated at the beginning of the cycle
          // as it does not change with the newly sampled values
          X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
          for( t=0; t<T_net; t++ ) {
            X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = u_ith_shared_mat.submat( aux_uvec_3, aux_uvec_2+t );
          }
          tau_h_diag.diag() = tau_h_shared_receive;
        }
        X = X_all;
        X.shed_rows(T_net*i,T_net*(i+1)-1);
        X = repmat(X,K_net,1);
        X_sp = arma::sp_mat(X);
        
      } else {
        throw std::range_error("direction not supported");
      }
      
      // identifies valid obs
      Y.elem( find_nonfinite(Y) ).zeros(); // change NAs to zero
      valid_obs = find(Y); // find elements different than zero
      Y_valid = Y.rows(valid_obs);
      X_sp_valid = arma::sp_mat(X.rows(valid_obs));
      sigma_Y_inv_valid = sigma_Y_inv.rows(valid_obs);
      sigma_Y_inv_valid_mat = arma::speye<arma::sp_mat>(sigma_Y_inv_valid.n_rows,sigma_Y_inv_valid.n_rows);
      sigma_Y_inv_valid_mat.diag() = sigma_Y_inv_valid;
      
      // The prior for the latent coordinates is uv_t_sigma_prior_inv repeated for each latent dimension
      uv_i_cov_prior = kron( tau_h_diag, uv_t_sigma_prior_inv );
      
      // Marginal Posterior
      // Covariance
      uv_i_cov_inv = X_sp_valid.t() * (sigma_Y_inv_valid_mat * X_sp_valid) + uv_i_cov_prior ;
      uv_i_cov = arma::inv_sympd(uv_i_cov_inv);
      
      if( dir==0 ) {
        C = linpred - X_sp * trans(u_ith_shared_mat.row(i));
        C_valid = C.rows(valid_obs);
        
        // Marginal Posterior
        // Mean
        uv_i_mean = uv_i_cov * ( X_sp_valid.t() * (sigma_Y_inv_valid % (Y_valid-C_valid)) );
        
        // Sampling u_ith_shared_mat
        u_ith_shared_mat.row(i) = trans(arma::mvnrnd( uv_i_mean , uv_i_cov ));
        
        // Recalculate linpred with the new values of u_ith_shared_mat
        linpred = X_sp * trans(u_ith_shared_mat.row(i)) + C;
        // Redefine mu_ijt with the new values of linpred
        for( k=0; k<K_net; k++ ) {
          aux_mat_1 = linpred.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1);
          aux_mat_1.reshape(T_net,(V_net-1));
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
          
          mu_ijt = mu_ijtk(k);
          mu_ijt.subcube(i,0,0, i,V_net-1,T_net-1) = aux_mat_1;
          mu_ijtk(k)= mu_ijt;
        }
        
      } else if( dir==1 ) {
        C = linpred - X_sp * trans(v_ith_shared_mat.row(i));
        C_valid = C.rows(valid_obs);
        
        // Marginal Posterior
        // Mean
        uv_i_mean = uv_i_cov * ( X_sp_valid.t() * (sigma_Y_inv_valid % (Y_valid-C_valid)) );
        
        // Sampling v_ith_shared_mat
        v_ith_shared_mat.row(i) = trans(arma::mvnrnd( uv_i_mean , uv_i_cov ));
        
        // Recalculate linpred with the new values of v_ith_shared_mat
        linpred = X_sp * trans(v_ith_shared_mat.row(i)) + C;
        // Redefine mu_ijt with the new values of linpred
        for( k=0; k<K_net; k++ ) {
          aux_mat_1 = linpred.rows((V_net-1)*T_net*k, (V_net-1)*T_net*(k+1)-1);
          aux_mat_1.reshape(T_net,(V_net-1));
          aux_mat_1 = aux_mat_1.t();
          aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
          
          mu_ijt = mu_ijtk(k);
          mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_1;
          mu_ijtk(k)= mu_ijt;
        }
      } else {
        throw std::range_error("direction not supported");
      }
    }
  }
  
  // get uv_ith_shared from uv_ith_shared_mat
  u_ith_shared.slice(0)=u_ith_shared_mat;
  u_ith_shared = reshape(u_ith_shared,V_net,T_net,H_dim);
  
  v_ith_shared.slice(0)=v_ith_shared_mat;
  v_ith_shared = reshape(v_ith_shared,V_net,T_net,H_dim);
  
  return Rcpp::List::create( Rcpp::Named("u_ith_shared") = u_ith_shared,
                             Rcpp::Named("v_ith_shared") = v_ith_shared,
                             Rcpp::Named("mu_ijtk") = mu_ijtk );
}

