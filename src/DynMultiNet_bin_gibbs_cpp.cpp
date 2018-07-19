#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "DynMultiNet_shared.h"

// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
arma::mat sample_x_iht_mat_DynMultiNet_bin_cpp( arma::mat x_iht_mat,
                                                const arma::mat x_t_sigma_prior_inv,
                                                const arma::mat tau_h,
                                                const arma::cube y_ijt,
                                                const arma::cube w_ijt,
                                                const arma::cube s_ijt,
                                                const arma::mat mu_t ) {
  unsigned int i=0;
  unsigned int j=0;
  unsigned int t=0;
  
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = x_iht_mat.n_cols / T_net;
  
  // Rcpp::Rcout << "V_net:" << V_net << std::endl ;
  // Rcpp::Rcout << "T_net:" << T_net << std::endl ;
  // Rcpp::Rcout << "H_dim:" << H_dim << std::endl ;
  
  arma::mat y_i = arma::zeros<arma::mat>((V_net-1)*T_net,1);
  arma::mat w_i = arma::zeros<arma::mat>((V_net-1)*T_net,1);
  arma::mat s_i = arma::zeros<arma::mat>((V_net-1)*T_net,1);
  
  arma::mat w_diag_i = arma::zeros<arma::mat>((V_net-1)*T_net,(V_net-1)*T_net);
  
  arma::mat x_tilde_i = arma::zeros<arma::mat>((V_net-1)*T_net,T_net*H_dim);
  
  arma::mat x_i = arma::zeros<arma::mat>(T_net*H_dim,1);
  arma::mat x_i_cov = arma::zeros<arma::mat>(T_net*H_dim,T_net*H_dim);
  
  arma::cube aux_cube_1;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  
  // for( i=0; i<V_net; i++ ) {
  i=3;
  
  aux_mat_1 = y_ijt.subcube(i,0,0, i,i,T_net-1);
  aux_mat_2 = y_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
  y_i = join_vert(aux_mat_1,aux_mat_2);
  y_i.shed_rows(i,i+1);
  y_i = y_i.t();
  y_i.reshape((V_net-1)*T_net,1);
  
  aux_mat_1 = w_ijt.subcube(i,0,0, i,i,T_net-1);
  aux_mat_2 = w_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
  w_i = join_vert(aux_mat_1,aux_mat_2);
  w_i.shed_rows(i,i+1);
  w_i = w_i.t();
  w_i.reshape((V_net-1)*T_net,1);
  
  aux_mat_1 = s_ijt.subcube(i,0,0, i,i,T_net-1);
  aux_mat_2 = s_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
  s_i = join_vert(aux_mat_1,aux_mat_2);
  s_i.shed_rows(i,i+1);
  s_i = s_i.t();
  s_i.reshape((V_net-1)*T_net,1);
  
  // }
  
  return x_iht_mat;
}

// sample_x_iht_mat_DynMultiNet_bin <- function( x_iht_mat,
//                                               x_t_sigma_prior_inv, tau_h,
//                                               y_ijt, w_ijt, s_ijt, mu_t ) {
//   
//   ### For each unit, block-sample the set of time-varying latent coordinates x_iht ###
//   
//   V_net <- dim(y_ijt)[1]
//   T_net <- dim(y_ijt)[3]
//   K_net <- dim(y_ijt)[4]
//   
//   H_dim <- dim(x_iht_mat)[2]/T_net
//   
//     for(i in 1:V_net) {
//       y_i <- c(t(rbind( matrix(y_ijt[i,1:i,],i,T_net)[-i,],
//                         matrix(y_ijt[i:V_net,i,],V_net-i+1,T_net)[-1,] ) ))
//       y_i <- matrix(y_i)
//       w_i <- c(t(rbind( matrix(w_ijt[i,1:i,],i,T_net)[-i,],
//                         matrix(w_ijt[i:V_net,i,],V_net-i+1,T_net)[-1,] ) ))
//       w_diag_i <- diag(w_i)
//       
//       x_tilde_i <- x_iht_mat[rep((1:V_net)[-i],each=T_net),]
//       
//       x_tilde_i_rm <- matrix(TRUE,T_net,T_net*H_dim)
//       for(t in 1:T_net){
//         x_tilde_i_rm[t,seq(t,T_net*H_dim,by=T_net)] <- F
//       }
//       x_tilde_i_rm <- do.call(rbind, replicate(V_net-1, x_tilde_i_rm, simplify=FALSE))
//       
//       x_tilde_i[x_tilde_i_rm] <- 0
//       
//       x_i_cov_inv <- t(x_tilde_i) %*% w_diag_i %*% x_tilde_i + kronecker( diag(as.numeric(tau_h)), x_t_sigma_prior_inv )
//       #browser()
//       x_i_cov <- solve( x_i_cov_inv )
//       # x_i_cov <- chol2inv(chol(x_i_cov_inv))
//       #all.equal( solve( x_i_cov_inv ) , chol2inv(chol(x_i_cov_inv)) )
//       
//       if(!isSymmetric(x_i_cov)) {x_i_cov[upper.tri(x_i_cov)] <- t(x_i_cov)[upper.tri(x_i_cov)]}
//       x_i_mean <- x_i_cov %*% ( t(x_tilde_i) %*% ( y_i - kronecker(matrix(1,V_net-1,1),matrix(1,T_net,1)*0.5) - w_diag_i %*% kronecker(matrix(1,V_net-1,1),mu_t) ) )
//       x_i <- mvtnorm::rmvnorm( n=1,
//                                mean=x_i_mean,
//                                sigma=x_i_cov )
//       x_iht_mat[i,] <- x_i
//     }
//   return(x_iht_mat)
// }
