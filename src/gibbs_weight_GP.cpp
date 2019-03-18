#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "DynMultiNet_shared.h"

// [[Rcpp::interfaces(r, cpp)]]


// [[Rcpp::export]]
Rcpp::List sample_baseline_tk_weight_cpp( arma::colvec theta_t,
                                          
                                          const arma::cube y_ijt,
                                          arma::cube mu_ijt,
                                          const double sigma_k,
                                          
                                          const arma::mat theta_t_cov_prior_inv,
                                          
                                          const bool lat_mean,
                                          double theta_t_bar,
                                          const double sigma_theta_bar,
                                          
                                          const bool directed=false ) {
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  
  // Auxiliar objects
  unsigned int i=0;
  arma::uvec aux_uvec_1;
  arma::colvec aux_colvec;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::cube aux_cube_1;
  
  double theta_t_bar_mean;
  double theta_t_bar_var;
  
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
  
  // prior mean vector for theta_t
  arma::colvec theta_t_mean_prior = arma::zeros<arma::colvec>(T_net);
  
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
  Y.replace(0, arma::datum::nan); // replace 0 with NA
  valid_obs = find_finite(Y); // find valid elements
  Y_valid = Y.rows(valid_obs);
  X_sp_valid = arma::sp_mat(X.rows(valid_obs));
  C_valid = C.rows(valid_obs);
  
  // GP prior mean
  if(!lat_mean){
    theta_t_bar=0;
  }
  theta_t_mean_prior.fill(theta_t_bar);
  
  // Marginal Posterior
  // Covariance
  theta_t_cov_inv = (1/pow(sigma_k,2)) * X_sp_valid.t() * X_sp_valid + theta_t_cov_prior_inv ;
  theta_t_cov = arma::inv_sympd(theta_t_cov_inv);
  // Mean
  theta_t_mean = theta_t_cov * ( (1/pow(sigma_k,2)) * X_sp_valid.t() * (Y_valid-C_valid) + theta_t_cov_prior_inv * theta_t_mean_prior );
  
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
  
  // Sample GP prior mean
  if(lat_mean){
    aux_colvec = arma::randn(1);
    theta_t_bar = aux_colvec(0);
    
    theta_t_bar_var = 1/( accu(theta_t_cov_prior_inv)+pow(sigma_theta_bar,-2) );
    
    aux_mat_1 = arma::ones<arma::mat>(1,T_net);
    aux_colvec = theta_t_bar_var * (aux_mat_1*theta_t_cov_prior_inv*theta_t);
    theta_t_bar_mean = aux_colvec(0);
    
    theta_t_bar = theta_t_bar_mean + theta_t_bar * sqrt(theta_t_bar_var);
  }
  
  return Rcpp::List::create( Rcpp::Named("theta_t") = theta_t,
                             Rcpp::Named("mu_ijt") = mu_ijt,
                             Rcpp::Named("theta_t_bar") = theta_t_bar );
}


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_weight_cpp( arma::cube uv_ith,
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
  unsigned int H_dim = uv_ith.n_slices;
  
  uv_ith = reshape(uv_ith,V_net,T_net*H_dim,1);
  arma::mat uv_ith_mat = uv_ith.slice(0);
  
  arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::colvec linpred = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::colvec C = arma::zeros<arma::colvec>((V_net-1)*T_net);
  
  arma::mat tau_h_diag(tau_h.n_rows,tau_h.n_rows); tau_h_diag.eye();
  tau_h_diag.diag() = tau_h;
  
  arma::mat uv_i_cov_prior = kron( tau_h_diag, uv_t_sigma_prior_inv );
  arma::mat uv_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat uv_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  
  arma::colvec uv_i_mean = arma::zeros<arma::colvec>(T_net*H_dim);
  
  // Performing sampling in uv_ith...
  // The rows in uv_ith_mat will act as the targeted "coefficients"
  // The columns in X_all will act as the "covariates"
  // The order of the response mu_ijth will depend on the order of rows in X_all
  
  arma::mat X_all = arma::zeros<arma::mat>(V_net*T_net,T_net*H_dim);
  arma::uvec aux_uvec_1 = T_net * arma::regspace<arma::uvec>( 0, V_net-1 );
  arma::uvec aux_uvec_2 = T_net * arma::regspace<arma::uvec>( 0, H_dim-1 );
  arma::uvec aux_uvec_3 = arma::regspace<arma::uvec>( 0, V_net-1 );
  
  arma::mat X = X_all;
  arma::sp_mat X_sp = arma::sp_mat(X);
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  arma::colvec Y_valid = Y;
  arma::sp_mat X_sp_valid = X_sp;
  arma::colvec C_valid = C;
  
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
      X_all.submat( aux_uvec_1+t , aux_uvec_2+t ) = uv_ith_mat.submat( aux_uvec_3, aux_uvec_2+t );
    }
    X = X_all;
    X.shed_rows(T_net*i,T_net*(i+1)-1);
    X_sp = arma::sp_mat(X);
    
    C = linpred - X_sp * trans(uv_ith_mat.row(i));
    
    // identifies valid obs
    Y.replace(0, arma::datum::nan); // replace 0 with NA
    valid_obs = find_finite(Y); // find valid elements
    Y_valid = Y.rows(valid_obs);
    X_sp_valid = arma::sp_mat(X.rows(valid_obs));
    C_valid = C.rows(valid_obs);
    
    // Marginal Posterior
    // Covariance
    uv_i_cov_inv = (1/pow(sigma_k,2)) * X_sp_valid.t() * X_sp_valid + uv_i_cov_prior ;
    uv_i_cov = arma::inv_sympd(uv_i_cov_inv);
    
    // Marginal Posterior
    // Mean
    uv_i_mean = uv_i_cov * ( (1/pow(sigma_k,2)) * X_sp_valid.t() * (Y_valid-C_valid) );
    
    // Sampling uv_ith_mat
    uv_ith_mat.row(i) = trans(arma::mvnrnd( uv_i_mean , uv_i_cov ));
    
    // Recalculate linpred with the new values of uv_ith_mat
    linpred = X_sp * trans(uv_ith_mat.row(i)) + C;
    // Redefine mu_ijt with the new values of linpred
    aux_mat_1 = linpred;
    aux_mat_1.reshape(T_net,(V_net-1));
    aux_mat_1 = aux_mat_1.t();
    aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
    mu_ijt.subcube(i,0,0, i,i,T_net-1) = aux_mat_1.rows(0,i);
    mu_ijt.subcube(i,i,0, V_net-1,i,T_net-1) = aux_mat_1.rows(i,V_net-1);
  }
  
  // get uv_ith from uv_ith_mat
  uv_ith.slice(0)=uv_ith_mat;
  uv_ith = reshape(uv_ith,V_net,T_net,H_dim);
  
  return Rcpp::List::create( Rcpp::Named("uv_ith") = uv_ith,
                             Rcpp::Named("mu_ijt") = mu_ijt );
}


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_shared_weight_cpp( arma::cube uv_ith_shared,
                                               const arma::mat uv_t_sigma_prior_inv,
                                               const arma::colvec tau_h,
                                               const arma::field<arma::cube> y_ijtk,
                                               arma::field<arma::cube> mu_ijtk,
                                               const arma::colvec sigma_k ) {
  
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
  
  // Variance associated with each observation in Y
  arma::colvec sigma_Y_inv = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  for( k=0; k<K_net; k++ ) {
    sigma_Y_inv.rows( k*(V_net-1)*T_net,(k+1)*(V_net-1)*T_net-1 ).fill( 1/sigma_k(k) );
  }
  
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
  
  arma::mat X = X_all;
  arma::sp_mat X_sp = arma::sp_mat(X_all);
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  arma::colvec Y_valid = Y;
  arma::sp_mat X_sp_valid = X_sp;
  arma::colvec C_valid = C;
  arma::colvec sigma_Y_inv_valid = sigma_Y_inv;
  arma::sp_mat sigma_Y_inv_valid_mat = arma::speye<arma::sp_mat>(1,1);
  
  for( i=0; i<V_net; i++ ) {
    // Rcpp::Rcout << "i=" << i << std::endl;
    arma::colvec Y = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
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
    X = X_all;
    X.shed_rows(T_net*i,T_net*(i+1)-1);
    X = repmat(X,K_net,1);
    X_sp = arma::sp_mat(X);
    
    C = linpred - X_sp * trans(uv_ith_shared_mat.row(i));
    
    // identifies valid obs
    Y.replace(0, arma::datum::nan); // replace 0 with NA
    valid_obs = find_finite(Y); // find valid elements
    Y_valid = Y.rows(valid_obs);
    X_sp_valid = arma::sp_mat(X.rows(valid_obs));
    C_valid = C.rows(valid_obs);
    
    X_sp_valid = arma::sp_mat(X.rows(valid_obs));
    sigma_Y_inv_valid = sigma_Y_inv.rows(valid_obs);
    sigma_Y_inv_valid_mat = arma::speye<arma::sp_mat>(sigma_Y_inv_valid.n_rows,sigma_Y_inv_valid.n_rows);
    sigma_Y_inv_valid_mat.diag() = sigma_Y_inv_valid;
    
    // Marginal Posterior
    // Covariance
    uv_i_cov_inv = X_sp_valid.t() * (sigma_Y_inv_valid_mat * X_sp_valid) + uv_i_cov_prior ;
    uv_i_cov = arma::inv_sympd(uv_i_cov_inv);
    
    // Marginal Posterior
    // Mean
    uv_i_mean = uv_i_cov * ( X_sp_valid.t() * (sigma_Y_inv_valid % (Y_valid-C_valid)) );
    
    // Sampling uv_ith_shared_mat
    uv_ith_shared_mat.row(i) = trans(arma::mvnrnd( uv_i_mean , uv_i_cov ));
    
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


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_weight_dir_cpp( arma::cube u_ith,
                                            arma::cube v_ith,
                                            
                                            const arma::cube y_ijt,
                                            arma::cube mu_ijt,
                                            const double sigma_k,
                                            
                                            const arma::mat uv_t_sigma_prior_inv,
                                            const bool lat_mean,
                                            arma::mat u_ith_bar,
                                            arma::mat v_ith_bar,
                                            const double sigma_uv_bar,
                                            
                                            const arma::colvec tau_h_send,
                                            const arma::colvec tau_h_receive ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  unsigned int h=0;
  
  unsigned int dir=0;
  
  arma::colvec aux_colvec;
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
  
  // prior mean vector for uv
  arma::colvec uv_i_mean_prior = arma::zeros<arma::colvec>(T_net*H_dim);
  
  // Posterior mean vector for uv
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
  
  // GP prior mean equal to zero
  if(!lat_mean){
    u_ith_bar.fill(0);
    v_ith_bar.fill(0);
  }
  arma::vec u_ith_mean = arma::zeros<arma::vec>(1);
  arma::vec v_ith_mean = arma::zeros<arma::vec>(1);
  double uv_ith_var=1;
  
  aux_colvec = arma::ones<arma::colvec>(T_net);
  
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
      Y.replace(0, arma::datum::nan); // replace 0 with NA
      valid_obs = find_finite(Y); // find valid elements
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
        uv_i_mean_prior = kron( u_ith_bar.row(i).t(), aux_colvec );
        uv_i_mean = uv_i_cov * ( (1/pow(sigma_k,2)) * X_sp_valid.t() * (Y_valid-C_valid) + uv_i_cov_prior * uv_i_mean_prior );
        
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
        uv_i_mean_prior = kron( v_ith_bar.row(i).t(), aux_colvec );
        uv_i_mean = uv_i_cov * ( (1/pow(sigma_k,2)) * X_sp_valid.t() * (Y_valid-C_valid) + uv_i_cov_prior * uv_i_mean_prior );
        
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
  
  // Sample GP prior mean 
  if(lat_mean){
    u_ith_bar.randn();
    
    v_ith_bar.randn();
    
    uv_ith_var = 1/( accu(uv_t_sigma_prior_inv)+pow(sigma_uv_bar,-2) );
    
    aux_mat_1 = arma::ones<arma::mat>(1,T_net);
    aux_mat_1 = aux_mat_1*uv_t_sigma_prior_inv;
    for( i=0; i<V_net; i++ ) {
      // Rcpp::Rcout << " i=" << i << std::endl;
      for( h=0; h<H_dim; h++ ) {
        // Rcpp::Rcout << " h=" << h << std::endl;
        
        // Sender
        aux_mat_2 = u_ith.subcube(i,0,h, i,T_net-1,h);
        u_ith_mean = uv_ith_var * (aux_mat_1*aux_mat_2.t());
        u_ith_bar(i,h) = u_ith_mean(0) + u_ith_bar(i,h) * sqrt(uv_ith_var);
        
        // Receiver
        aux_mat_2 = v_ith.subcube(i,0,h, i,T_net-1,h);
        v_ith_mean = uv_ith_var * (aux_mat_1*aux_mat_2.t());
        v_ith_bar(i,h) = v_ith_mean(0) + v_ith_bar(i,h) * sqrt(uv_ith_var);
      }
    }
    
  }
  
  return Rcpp::List::create( Rcpp::Named("u_ith") = u_ith,
                             Rcpp::Named("v_ith") = v_ith,
                             Rcpp::Named("mu_ijt") = mu_ijt,
                             Rcpp::Named("u_ith_bar") = u_ith_bar,
                             Rcpp::Named("v_ith_bar") = v_ith_bar );
}


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_shared_weight_dir_cpp( arma::cube u_ith_shared,
                                                   arma::cube v_ith_shared,
                                                   
                                                   const arma::field<arma::cube> y_ijtk,
                                                   arma::field<arma::cube> mu_ijtk,
                                                   const arma::colvec sigma_k,
                                                     
                                                   const arma::mat uv_t_sigma_prior_inv,
                                                   const bool lat_mean,
                                                   arma::mat u_ith_shared_bar,
                                                   arma::mat v_ith_shared_bar,
                                                   const double sigma_uv_bar,
                                                   
                                                   const arma::colvec tau_h_shared_send,
                                                   const arma::colvec tau_h_shared_receive ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int k=0;
  unsigned int t=0;
  unsigned int h=0;
  
  unsigned int dir=0;
  
  arma::colvec aux_colvec;
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
  for( k=0; k<K_net; k++ ) {
    sigma_Y_inv.rows( k*(V_net-1)*T_net,(k+1)*(V_net-1)*T_net-1 ).fill( 1/sigma_k(k) );
  }
  
  arma::mat tau_h_diag(tau_h_shared_send.n_rows,tau_h_shared_send.n_rows); tau_h_diag.eye();
  arma::mat uv_i_cov_prior = kron( tau_h_diag, uv_t_sigma_prior_inv );
  arma::mat uv_i_cov_inv = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  arma::mat uv_i_cov = arma::zeros<arma::mat>( T_net*H_dim, T_net*H_dim );
  
  // prior mean vector for uv
  arma::colvec uv_i_mean_prior = arma::zeros<arma::colvec>(T_net*H_dim);
  
  // Posterior mean vector for uv
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
  
  // GP prior mean equal to zero
  if(!lat_mean){
    u_ith_shared_bar.fill(0);
    v_ith_shared_bar.fill(0);
  }
  arma::vec u_ith_shared_mean = arma::zeros<arma::vec>(1);
  arma::vec v_ith_shared_mean = arma::zeros<arma::vec>(1);
  double uv_ith_shared_var=1;
  
  aux_colvec = arma::ones<arma::colvec>(T_net);
  
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
        // Updating Y, linpred as vectors, lower triangular adjacency
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
        // Updating Y, linpred as vectors, upper triangular adjacency
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
      Y.replace(0, arma::datum::nan); // replace 0 with NA
      valid_obs = find_finite(Y); // find valid elements
      Y_valid = Y.rows(valid_obs);
      X_sp_valid = arma::sp_mat(X.rows(valid_obs));
      sigma_Y_inv_valid = sigma_Y_inv.rows(valid_obs);
      sigma_Y_inv_valid_mat = arma::speye<arma::sp_mat>(sigma_Y_inv_valid.n_rows,sigma_Y_inv_valid.n_rows);
      sigma_Y_inv_valid_mat.diag() = sigma_Y_inv_valid;
      
      // Marginal Posterior
      // Covariance
      uv_i_cov_inv = X_sp_valid.t() * (sigma_Y_inv_valid_mat * X_sp_valid) + uv_i_cov_prior ;
      uv_i_cov = arma::inv_sympd(uv_i_cov_inv);
      
      if( dir==0 ) {
        C = linpred - X_sp * trans(u_ith_shared_mat.row(i));
        C_valid = C.rows(valid_obs);
        
        // Marginal Posterior
        // Mean
        uv_i_mean_prior = kron( u_ith_shared_bar.row(i).t(), aux_colvec );
        uv_i_mean = uv_i_cov * ( X_sp_valid.t() * (sigma_Y_inv_valid % (Y_valid-C_valid)) + uv_i_cov_prior * uv_i_mean_prior );
        
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
        uv_i_mean_prior = kron( v_ith_shared_bar.row(i).t(), aux_colvec );
        uv_i_mean = uv_i_cov * ( X_sp_valid.t() * (sigma_Y_inv_valid % (Y_valid-C_valid)) + uv_i_cov_prior * uv_i_mean_prior );
        
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
  u_ith_shared.slice(0) = u_ith_shared_mat;
  u_ith_shared = reshape(u_ith_shared,V_net,T_net,H_dim);
  
  v_ith_shared.slice(0) = v_ith_shared_mat;
  v_ith_shared = reshape(v_ith_shared,V_net,T_net,H_dim);
  
  // Sample GP prior mean 
  if(lat_mean){
    u_ith_shared_bar.randn();
    
    v_ith_shared_bar.randn();
    
    uv_ith_shared_var = 1/( accu(uv_t_sigma_prior_inv)+pow(sigma_uv_bar,-2) );
    
    aux_mat_1 = arma::ones<arma::mat>(1,T_net);
    aux_mat_1 = aux_mat_1*uv_t_sigma_prior_inv;
    for( i=0; i<V_net; i++ ) {
      // Rcpp::Rcout << " i=" << i << std::endl;
      for( h=0; h<H_dim; h++ ) {
        // Rcpp::Rcout << " h=" << h << std::endl;
        
        // Sender
        aux_mat_2 = u_ith_shared.subcube(i,0,h, i,T_net-1,h);
        u_ith_shared_mean = uv_ith_shared_var * (aux_mat_1*aux_mat_2.t());
        u_ith_shared_bar(i,h) = u_ith_shared_mean(0) + u_ith_shared_bar(i,h) * sqrt(uv_ith_shared_var);
        
        // Receiver
        aux_mat_2 = v_ith_shared.subcube(i,0,h, i,T_net-1,h);
        v_ith_shared_mean = uv_ith_shared_var * (aux_mat_1*aux_mat_2.t());
        v_ith_shared_bar(i,h) = v_ith_shared_mean(0) + v_ith_shared_bar(i,h) * sqrt(uv_ith_shared_var);
      }
    }
    
  }
  
  return Rcpp::List::create( Rcpp::Named("u_ith_shared") = u_ith_shared,
                             Rcpp::Named("v_ith_shared") = v_ith_shared,
                             Rcpp::Named("mu_ijtk") = mu_ijtk,
                             Rcpp::Named("u_ith_shared_bar") = u_ith_shared_bar,
                             Rcpp::Named("v_ith_shared_bar") = v_ith_shared_bar );
}


// [[Rcpp::export]]
Rcpp::List sample_add_eff_it_weight_cpp( arma::colvec sp_it,
                                         const arma::mat sp_t_cov_prior_inv,
                                         const arma::cube y_ijt,
                                         arma::cube mu_ijt,
                                         const double sigma_k,
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
  
  // Objects for the calculation
  arma::colvec Y = arma::zeros<arma::colvec>(1);
  arma::colvec linpred = arma::zeros<arma::colvec>(1);
  
  arma::colvec C = arma::zeros<arma::colvec>(1);
  arma::colvec C_valid = C;
  
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
      aux_mat_2 = mu_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
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
      aux_mat_2 = mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
      aux_mat_2.shed_row(i);
      aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
    }
    aux_mat_1.shed_row(0);
    aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
    linpred = aux_mat_1;
    
  }
  
  // identifies valid obs
  Y.replace(0, arma::datum::nan); // replace 0 with NA
  valid_obs = find_finite(Y); // find valid elements
  
  // Keep only valid cases
  Y = Y.rows(valid_obs);
  X_sp_valid = arma::sp_mat(X.rows(valid_obs));
  
  for( dir=0; dir<2; dir++ ) { // dir=0 samples p ; dir=1 samples s
    if( (dir==1) | (directed&(dir==0)) ) { // if directed, sample s and p. if undirected sample only s
      // Rcpp::Rcout << "dir=" << dir << std::endl;
      
      if( dir==0 ){
        // Marginal Posterior Covariance
        sp_it_cov_inv = (1/pow(sigma_k,2)) * X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1).t() * X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1) ;
        sp_it_cov_inv += sp_it_cov_prior_inv;
        sp_it_cov = arma::inv_sympd(sp_it_cov_inv);
        // Marginal Posterior Mean
        C = linpred - (X_sp.cols(T_net*V_net,2*T_net*V_net-1) * sp_it.rows(T_net*V_net,2*T_net*V_net-1));
        C_valid = C.rows(valid_obs);
        sp_it_mean = sp_it_cov * ( (1/pow(sigma_k,2)) * X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1).t() * (Y-C_valid) );
        
        // Sampling sp_it
        sp_it.rows(T_net*V_net,2*T_net*V_net-1) = arma::mvnrnd( sp_it_mean , sp_it_cov );
        
        // Recalculate linpred with the new values of sp_it
        linpred = X_sp.cols(T_net*V_net,2*T_net*V_net-1) * sp_it.rows(T_net*V_net,2*T_net*V_net-1) + C;
        
      } else if( dir==1 ){
        // Marginal Posterior Covariance
        sp_it_cov_inv = (1/pow(sigma_k,2)) * X_sp_valid.cols(0,T_net*V_net-1).t() * X_sp_valid.cols(0,T_net*V_net-1);
        sp_it_cov_inv += sp_it_cov_prior_inv;
        sp_it_cov = arma::inv_sympd(sp_it_cov_inv);
        // Marginal Posterior Mean
        C = linpred - (X_sp.cols(0,T_net*V_net-1) * sp_it.rows(0,T_net*V_net-1));
        C_valid = C.rows(valid_obs);
        sp_it_mean = sp_it_cov * ( (1/pow(sigma_k,2)) * X_sp_valid.cols(0,T_net*V_net-1).t() * (Y-C_valid) );
        
        // Sampling sp_it
        sp_it.rows(0,T_net*V_net-1) = arma::mvnrnd( sp_it_mean , sp_it_cov );
        
        // Recalculate linpred with the new values of sp_it
        linpred = X_sp.cols(0,T_net*V_net-1) * sp_it.rows(0,T_net*V_net-1) + C;
        
      }
    }
  }
  // return sp_it;
  
  // Redefine mu_ijt with the new values of linpred
  if( !directed ){
    aux_mat_1 = linpred;
    aux_mat_1.reshape(V_net*(V_net-1)/2,T_net);
    for( i=1; i<V_net; i++ ) {
      aux_mat_2 = aux_mat_1.rows((i-1)*V_net-((i-1)*i)/2,i*V_net-(i*(i+1))/2-1);
      mu_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1) = aux_mat_2;
    }
  } else {
    aux_mat_1 = linpred;
    aux_mat_1.reshape(V_net*(V_net-1),T_net);
    for( i=0; i<V_net; i++ ) {
      aux_mat_2 = aux_mat_1.rows(i*(V_net-1),(i+1)*(V_net-1)-1);
      aux_mat_2.insert_rows(i,arma::zeros<arma::mat>(1,T_net));
      mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_2;
    }
  }
  
  return Rcpp::List::create( Rcpp::Named("sp_it") = sp_it,
                             Rcpp::Named("mu_ijt") = mu_ijt );
}


// [[Rcpp::export]]
Rcpp::List sample_add_eff_it_shared_weight_cpp( arma::colvec sp_it,
                                                const arma::mat sp_t_cov_prior_inv,
                                                const arma::field<arma::cube> y_ijtk,
                                                arma::field<arma::cube> mu_ijtk,
                                                const arma::colvec sigma_k,
                                                const bool directed=false ) {
  
  arma::cube y_ijt = y_ijtk(0);
  arma::cube mu_ijt = mu_ijtk(0);
  
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
  
  // Objects for the calculation
  arma::colvec Y = arma::zeros<arma::colvec>(1);
  arma::colvec linpred = arma::zeros<arma::colvec>(1);
  
  arma::colvec C = arma::zeros<arma::colvec>(1);
  arma::colvec C_valid = C;
  
  // Variance associated with each observation in Y
  arma::colvec sigma_Y_inv = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  if(!directed){
    sigma_Y_inv = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
    for( k=0; k<K_net; k++ ) {
      sigma_Y_inv.rows( k*T_net*V_net*(V_net-1)/2,(k+1)*T_net*V_net*(V_net-1)/2-1 ).fill( 1/sigma_k(k) );
    }
  } else {
    sigma_Y_inv = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
    for( k=0; k<K_net; k++ ) {
      sigma_Y_inv.rows( k*T_net*V_net*(V_net-1),(k+1)*T_net*V_net*(V_net-1)-1 ).fill( 1/sigma_k(k) );
    }
  }
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
  arma::colvec sigma_Y_inv_valid = sigma_Y_inv;
  arma::sp_mat sigma_Y_inv_valid_mat = arma::speye<arma::sp_mat>(1,1);
  
  // Marginal posterior covariance matrix for sp_it
  arma::mat sp_it_cov_inv=arma::zeros<arma::mat>(T_net*V_net,T_net*V_net);
  arma::mat sp_it_cov;
  arma::colvec sp_it_mean = arma::zeros<arma::colvec>(2*V_net*T_net);
  
  
  // Create vectors Y and linpred tht serve as response and  //
  if( !directed ){
    Y = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
    linpred = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
  } else {
    Y = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
    linpred = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
  }
  
  
  for( k=0; k<K_net; k++ ) {
    // Rcpp::Rcout << "k=" << k << std::endl;
    y_ijt = y_ijtk(k);
    mu_ijt = mu_ijtk(k);
    
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
        aux_mat_2 = mu_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
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
        aux_mat_2 = mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_2.shed_row(i);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
      linpred.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1) = aux_mat_1;
      
    }
  }
  
  // identifies valid obs
  Y.replace(0, arma::datum::nan); // replace 0 with NA
  valid_obs = find_finite(Y); // find valid elements
  
  // Keep only valid cases
  Y = Y.rows(valid_obs);
  X_sp_valid = arma::sp_mat(X.rows(valid_obs));
  sigma_Y_inv_valid = sigma_Y_inv.rows(valid_obs);
  sigma_Y_inv_valid_mat = arma::speye<arma::sp_mat>(sigma_Y_inv_valid.n_rows,sigma_Y_inv_valid.n_rows);
  sigma_Y_inv_valid_mat.diag() = sigma_Y_inv_valid;
  
  // SAMPLING sp_it //
  for( dir=0; dir<2; dir++ ) { // dir=0 samples p ; dir=1 samples s
    if( (dir==1) | (directed&(dir==0)) ) { // if directed, sample s and p. if undirected sample only s
      // Rcpp::Rcout << "dir=" << dir << std::endl;
      
      if( dir==0 ){
        // Marginal Posterior Covariance
        sp_it_cov_inv = X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1).t() * (sigma_Y_inv_valid_mat * X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1));
        sp_it_cov_inv += sp_it_cov_prior_inv;
        sp_it_cov = arma::inv_sympd(sp_it_cov_inv);
        
        // Marginal Posterior Mean
        C = linpred - (X_sp.cols(T_net*V_net,2*T_net*V_net-1) * sp_it.rows(T_net*V_net,2*T_net*V_net-1));
        C_valid = C.rows(valid_obs);
        sp_it_mean = sp_it_cov * ( X_sp_valid.cols(T_net*V_net,2*T_net*V_net-1).t() * (sigma_Y_inv_valid % (Y-C_valid)) );
        
        // Sampling sp_it
        sp_it.rows(T_net*V_net,2*T_net*V_net-1) = arma::mvnrnd( sp_it_mean , sp_it_cov );
        
        // Recalculate linpred with the new values of sp_it
        linpred = X_sp.cols(T_net*V_net,2*T_net*V_net-1) * sp_it.rows(T_net*V_net,2*T_net*V_net-1) + C;
        
      } else if( dir==1 ){
        // Marginal Posterior Covariance
        sp_it_cov_inv = X_sp_valid.cols(0,T_net*V_net-1).t() * (sigma_Y_inv_valid_mat * X_sp_valid.cols(0,T_net*V_net-1));
        sp_it_cov_inv += sp_it_cov_prior_inv;
        sp_it_cov = arma::inv_sympd(sp_it_cov_inv);
        // Marginal Posterior Mean
        C = linpred - (X_sp.cols(0,T_net*V_net-1) * sp_it.rows(0,T_net*V_net-1));
        C_valid = C.rows(valid_obs);
        sp_it_mean = sp_it_cov * ( X_sp_valid.cols(0,T_net*V_net-1).t() * (sigma_Y_inv_valid % (Y-C_valid)) );
        
        // Sampling sp_it
        sp_it.rows(0,T_net*V_net-1) = arma::mvnrnd( sp_it_mean , sp_it_cov );
        
        // Recalculate linpred with the new values of sp_it
        linpred = X_sp.cols(0,T_net*V_net-1) * sp_it.rows(0,T_net*V_net-1) + C;
        
      }
    }
  }
  // return sp_it;
  
  // Redefine mu_ijt with the new values of linpred
  for( k=0; k<K_net; k++ ) {
    // Rcpp::Rcout << "k=" << k << std::endl;
    mu_ijt = mu_ijtk(k);
    
    if( !directed ){
      aux_mat_1 = linpred.rows(k*T_net*V_net*(V_net-1)/2, (k+1)*T_net*V_net*(V_net-1)/2-1);
      aux_mat_1.reshape(V_net*(V_net-1)/2,T_net);
      for( i=1; i<V_net; i++ ) {
        aux_mat_2 = aux_mat_1.rows((i-1)*V_net-((i-1)*i)/2,i*V_net-(i*(i+1))/2-1);
        mu_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1) = aux_mat_2;
      }
    } else {
      aux_mat_1 = linpred.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1);
      aux_mat_1.reshape(V_net*(V_net-1),T_net);
      for( i=0; i<V_net; i++ ) {
        aux_mat_2 = aux_mat_1.rows(i*(V_net-1),(i+1)*(V_net-1)-1);
        aux_mat_2.insert_rows(i,arma::zeros<arma::mat>(1,T_net));
        mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_2;
      }
    }
    
    mu_ijtk(k) = mu_ijt;
  }
  
  return Rcpp::List::create( Rcpp::Named("sp_it") = sp_it,
                             Rcpp::Named("mu_ijtk") = mu_ijtk );
}


// [[Rcpp::export]]
Rcpp::List sample_coeff_tp_weight_cpp( arma::mat beta_tp,
                                       const arma::mat beta_t_cov_prior_inv,
                                       const arma::field<arma::cube> y_ijtk,
                                       arma::field<arma::cube> mu_ijtk,
                                       const arma::colvec sigma_k,
                                       const arma::mat x_ijtkp_mat,
                                       const bool directed=false ) {
  
  // This function block sample
  // beta_tp = (beta_tp[1,1],...,beta_tp[T_net,1],...,beta_tp[1,P_pred],...,beta_tp[T_net,P_pred])
  
  arma::cube y_ijt = y_ijtk(0);
  arma::cube mu_ijt = mu_ijtk(0);
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int K_net = y_ijtk.n_rows;
  unsigned int P_pred = beta_tp.n_cols;
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int k=0;
  arma::uvec aux_uvec_1;
  arma::colvec aux_vec;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::mat aux_mat_3;
  arma::cube aux_cube_1;
  
  // reshape beta
  beta_tp.reshape(T_net*P_pred,1);
  
  //  Design matrix
  arma::sp_mat X_sp = arma::sp_mat(x_ijtkp_mat);
  
  // Objects for the calculation
  arma::colvec Y = arma::zeros<arma::colvec>(1);
  arma::colvec linpred = arma::zeros<arma::colvec>(1);
  
  arma::colvec C = arma::zeros<arma::colvec>(1);
  arma::colvec C_valid = C;
  
  // Variance associated with each observation in Y
  arma::colvec sigma_Y_inv = arma::zeros<arma::colvec>((V_net-1)*T_net*K_net);
  if(!directed){
    sigma_Y_inv = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
    for( k=0; k<K_net; k++ ) {
      sigma_Y_inv.rows( k*T_net*V_net*(V_net-1)/2,(k+1)*T_net*V_net*(V_net-1)/2-1 ).fill( 1/sigma_k(k) );
    }
  } else {
    sigma_Y_inv = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
    for( k=0; k<K_net; k++ ) {
      sigma_Y_inv.rows( k*T_net*V_net*(V_net-1),(k+1)*T_net*V_net*(V_net-1)-1 ).fill( 1/sigma_k(k) );
    }
  }
  // arma::sp_mat Omega_sp;
  
  // Prior covariance
  arma::mat beta_tp_cov_prior_inv = kron( beta_t_cov_prior_inv, arma::eye<arma::mat>(P_pred,P_pred) );
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  arma::colvec sigma_Y_inv_valid = sigma_Y_inv;
  arma::sp_mat sigma_Y_inv_valid_mat = arma::speye<arma::sp_mat>(1,1);
  arma::sp_mat X_sp_valid = arma::sp_mat(x_ijtkp_mat);
  
  // Marginal posterior covariance matrix for sp_it
  arma::mat beta_tp_cov_inv = arma::zeros<arma::mat>(T_net*V_net,T_net*V_net);
  arma::mat beta_tp_cov;
  arma::colvec beta_tp_mean = arma::zeros<arma::colvec>(P_pred*T_net);
  
  // Create vectors Y and linpred that serve as response and linear predictor //
  if( !directed ){
    Y = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
    linpred = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1)/2);
  } else {
    Y = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
    linpred = arma::zeros<arma::colvec>(K_net*T_net*V_net*(V_net-1));
  }
  
  // Extract info for Y and linpred from cubes //
  for( k=0; k<K_net; k++ ) {
    // Rcpp::Rcout << "k=" << k << std::endl;
    y_ijt = y_ijtk(k);
    mu_ijt = mu_ijtk(k);
    
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
        aux_mat_2 = mu_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1);
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
        aux_mat_2 = mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_2.shed_row(i);
        aux_mat_1.insert_rows( aux_mat_1.n_rows, aux_mat_2 );
      }
      aux_mat_1.shed_row(0);
      aux_mat_1.reshape(T_net*V_net*(V_net-1),1);
      linpred.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1) = aux_mat_1;
      
    }
  }
  
  // identifies valid obs
  Y.replace(0, arma::datum::nan); // replace 0 with NA
  valid_obs = find_finite(Y); // find valid elements
  
  // Keep only valid cases
  Y = Y.rows(valid_obs);
  X_sp_valid = arma::sp_mat(x_ijtkp_mat.rows(valid_obs));
  sigma_Y_inv_valid = sigma_Y_inv.rows(valid_obs);
  sigma_Y_inv_valid_mat = arma::speye<arma::sp_mat>(sigma_Y_inv_valid.n_rows,sigma_Y_inv_valid.n_rows);
  sigma_Y_inv_valid_mat.diag() = sigma_Y_inv_valid;
  
  // SAMPLING beta //
  
  // Marginal Posterior Covariance
  beta_tp_cov_inv = X_sp_valid.t() * (sigma_Y_inv_valid_mat * X_sp_valid);
  beta_tp_cov_inv += beta_tp_cov_prior_inv;
  beta_tp_cov = arma::inv_sympd(beta_tp_cov_inv);
  // Marginal Posterior Mean
  C = linpred - (X_sp * beta_tp);
  C_valid = C.rows(valid_obs);
  beta_tp_mean = beta_tp_cov * ( X_sp_valid.t() * (sigma_Y_inv_valid % (Y-C_valid)) );
  
  // Sampling sp_it
  beta_tp = arma::mvnrnd( beta_tp_mean , beta_tp_cov );
  
  // Recalculate linpred with the new values of sp_it
  linpred = X_sp * beta_tp + C;
  
  // Redefine mu_ijt with the new values of linpred
  for( k=0; k<K_net; k++ ) {
    // Rcpp::Rcout << "k=" << k << std::endl;
    mu_ijt = mu_ijtk(k);
    
    if( !directed ){
      aux_mat_1 = linpred.rows(k*T_net*V_net*(V_net-1)/2, (k+1)*T_net*V_net*(V_net-1)/2-1);
      aux_mat_1.reshape(V_net*(V_net-1)/2,T_net);
      for( i=1; i<V_net; i++ ) {
        aux_mat_2 = aux_mat_1.rows((i-1)*V_net-((i-1)*i)/2,i*V_net-(i*(i+1))/2-1);
        mu_ijt.subcube(i,i-1,0, V_net-1,i-1,T_net-1) = aux_mat_2;
      }
    } else {
      aux_mat_1 = linpred.rows(k*T_net*V_net*(V_net-1), (k+1)*T_net*V_net*(V_net-1)-1);
      aux_mat_1.reshape(V_net*(V_net-1),T_net);
      for( i=0; i<V_net; i++ ) {
        aux_mat_2 = aux_mat_1.rows(i*(V_net-1),(i+1)*(V_net-1)-1);
        aux_mat_2.insert_rows(i,arma::zeros<arma::mat>(1,T_net));
        mu_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_2;
      }
    }
    
    mu_ijtk(k) = mu_ijt;
  }
  
  // reshape beta
  beta_tp.reshape(T_net,P_pred);
  
  return Rcpp::List::create( Rcpp::Named("beta_tp") = beta_tp,
                             Rcpp::Named("mu_ijtk") = mu_ijtk );
}


// [[Rcpp::export]]
double sample_var_weight_cpp( double sigma_k,
                              double sigma_k_prop_int,
                              const arma::cube y_ijt,
                              const arma::cube mu_ijt,
                              const bool directed=false ) {
  
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  
  // column matrix for the continuous response
  arma::colvec Y = arma::zeros<arma::colvec>(1);
  
  // column matrix with the linear predictor
  arma::colvec linpred = arma::zeros<arma::colvec>(1);
  
  // vector that identifies valid observation for the model
  // in this case, where Y!=0 and is not NA
  arma::uvec valid_obs;
  arma::colvec Y_valid = Y;
  arma::colvec linpred_valid = linpred;
  
  // Auxiliar objects
  unsigned int i=0;
  double aux_double;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  // Variables for MH
  double sigma_k_prop = sigma_k; // proposal value
  double rho = 0; // auxiliar value for the proposal random walk
  double prop_probs=0;
  double cur_probs=0;
  double MH_ratio;
  
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
  
  // identifies valid obs
  Y.replace(0, arma::datum::nan); // replace 0 with NA
  valid_obs = find_finite(Y); // find valid elements
  Y_valid = Y.rows(valid_obs);
  linpred_valid = linpred.rows(valid_obs);
  
  // Propose a new value for sigma_k
  rho = R::runif(1/sigma_k_prop_int,sigma_k_prop_int);
  sigma_k_prop = rho * sigma_k;
  
  for( i=0; i<Y_valid.n_rows; i++ ) {
    prop_probs += R::dnorm( Y_valid(i), linpred_valid(i), sigma_k_prop,1);
    cur_probs += R::dnorm( Y_valid(i), linpred_valid(i), sigma_k,1);
  }
  
  // Metropolis-Hastings ratio
  MH_ratio = prop_probs - cur_probs + log( sigma_k_prop / sigma_k );
  //Rcpp::Rcout << MH_ratio << std::endl ;
  aux_double = log( R::runif(0,1) );
  if( aux_double < MH_ratio ) {
    // New value accepted!
    // Rcpp::Rcout << "accept!" << std::endl ;
    
    // Updating current values of the MCMC
    sigma_k = sigma_k_prop;
  }
  
  return( sigma_k );
}


