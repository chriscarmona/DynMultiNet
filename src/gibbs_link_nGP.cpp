#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "DynMultiNet_shared.h"

// [[Rcpp::interfaces(r, cpp)]]



// [[Rcpp::export]]
Rcpp::List sample_baseline_t_link_nGP_cpp( const arma::colvec eta_t,
                                           arma::mat alpha_eta_t,
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
  
  // Checking dimensions
  if(eta_t.n_rows!=T_net){ throw std::range_error("Input size mismatch, check eta_t");}
  if((alpha_eta_t.n_rows!=3)|(alpha_eta_t.n_cols!=T_net)){ throw std::range_error("Input size mismatch, check alpha_eta_t");}
  
  
  // initialising objects
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
    C_t.col(t) = linpred_t.col(t) - ( X_t.slice(t) * alpha_eta_t.col(t) );
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
  arma::mat P1chol = arma::eye<arma::mat>(3,3); P1chol.diag().fill(10);
  
  alpha_eta_t = kfsim_cpp( Z_t,
                           
                           dd,
                           X_t,
                           Z_t_cov_chol,
                           
                           cc,
                           nGP_G_t,
                           nGP_H_t,
                           nGP_Wchol_t,
                           
                           a1,
                           P1chol );
  // Recalculate linpred with the new values of alpha_t
  for( t=0; t<T_net; t++ ) {
    linpred_t.col(t) = C_t.col(t) + ( X_t.slice(t) * alpha_eta_t.col(t) );
  }
  
  // updating value of eta
  // eta_t = trans(alpha_eta_t.row(0));
  
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
  
  return Rcpp::List::create( Rcpp::Named("eta_t") = trans(alpha_eta_t.row(0)),
                             Rcpp::Named("gamma_ijt") = gamma_ijt,
                             Rcpp::Named("alpha_eta_t") = alpha_eta_t );
}


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_link_nGP_cpp( const arma::cube ab_ith,
                                          arma::field<arma::cube> alpha_ab_ith,
                                          const arma::cube y_ijt,
                                          const arma::cube w_ijt,
                                          arma::cube gamma_ijt,
                                          const arma::field<arma::cube> nGP_G_t,
                                          const arma::field<arma::cube> nGP_H_t,
                                          const arma::field<arma::cube> nGP_Wchol_t,
                                          const bool directed=false ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  unsigned int h=0;
  unsigned int alpha_i=0;
  arma::uvec aux_uvec_1;
  arma::colvec aux_vec;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::cube aux_cube_1;
  
  // Rcpp::Rcout << "(0)" << std::endl;
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = ab_ith.n_slices;
  
  // initialising objects
  arma::mat Y_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat PG_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat linpred_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat C_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat Z_t = arma::zeros<arma::mat>(V_net-1,T_net);
  arma::cube X_t = arma::ones<arma::cube>(1,1,T_net); 
  
  arma::mat alpha_t = arma::zeros<arma::mat>(3,T_net);
  
  arma::cube Z_t_cov_chol = arma::zeros<arma::cube>(V_net-1,V_net-1,T_net);
  
  arma::mat dd = arma::zeros<arma::mat>(V_net-1,T_net);
  arma::mat cc = arma::zeros<arma::mat>(3,T_net);
  arma::colvec a1 = arma::zeros<arma::colvec>(3);
  arma::mat P1chol = arma::eye<arma::mat>(3,3); P1chol.diag().fill(10);
  
  // Repeat simulation for each agent
  for( i=0; i<V_net; i++ ) {
    // Network data
    aux_mat_1 = y_ijt.subcube(i,0,0, i,i,T_net-1);
    aux_mat_2 = y_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
    aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
    aux_mat_1.shed_rows(i,i+1);
    Y_t = aux_mat_1;
    
    // Polya-Gamma data
    aux_mat_1 = w_ijt.subcube(i,0,0, i,i,T_net-1);
    aux_mat_2 = w_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
    aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
    aux_mat_1.shed_rows(i,i+1);
    PG_t = aux_mat_1;
    
    // linear predictor
    aux_mat_1 = gamma_ijt.subcube(i,0,0, i,i,T_net-1);
    aux_mat_2 = gamma_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
    aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
    aux_mat_1.shed_rows(i,i+1);
    linpred_t = aux_mat_1;
    
    // Repeat simulation for each agent and each dimension in its coordinates
    for( h=0; h<H_dim; h++ ) {
      alpha_t = arma::zeros<arma::mat>(3,T_net);
      alpha_t.row(0) = alpha_ab_ith(0).slice(h).row(i);
      // Model matrix
      // Has three columns for: U, U', A, only the first has effect on Y
      X_t = arma::zeros<arma::cube>(V_net,3,T_net);
      X_t.col(0) = alpha_ab_ith(0).slice(h);
      X_t.shed_row(i);
      
      // Constant term for theta in the linear predictor
      C_t = linpred_t; C_t.zeros();
      for( t=0; t<T_net; t++ ) {
        // Rcpp::Rcout << "t=" << t << std::endl;
        C_t.col(t) = linpred_t.col(t) - ( X_t.slice(t) * alpha_t.col(t) );
      }
      
      // Adjusted response
      Z_t = (Y_t-0.5)/PG_t - C_t;
      Z_t_cov_chol = arma::zeros<arma::cube>(Z_t.n_rows,Z_t.n_rows,T_net);
      for( t=0; t<T_net; t++ ) {
        Z_t_cov_chol.slice(t).diag() = pow( PG_t.col(t) , -1);
      }
      
      // Sampling coefficient using State-Space model and simulation smoother
      alpha_t = kfsim_cpp( Z_t,
                           
                           dd,
                           X_t,
                           Z_t_cov_chol,
                           
                           cc,
                           nGP_G_t(i),
                           nGP_H_t(i),
                           nGP_Wchol_t(i),
                           
                           a1,
                           P1chol );
      
      // store updated values of alpha_ab_ith
      for(alpha_i=0; alpha_i<3; alpha_i++ ) {
        alpha_ab_ith(alpha_i).slice(h).row(i) = alpha_t.row(alpha_i);
      }
      // store updated values of ab_ith
      // ab_ith.slice(h).row(i)=alpha_t.row(0);
      
      // Recalculate linpred with the new values of alpha_t
      for( t=0; t<T_net; t++ ) {
        linpred_t.col(t) = C_t.col(t) + ( X_t.slice(t) * alpha_t.col(t) );
      }
      
      // Update gamma_ijt with the new values of linpred
      aux_mat_1 = linpred_t;
      aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
      gamma_ijt.subcube(i,0,0, i,i,T_net-1) = aux_mat_1.rows(0,i);
      gamma_ijt.subcube(i,i,0, V_net-1,i,T_net-1) = aux_mat_1.rows(i,V_net-1);
    }
  }
  
  return Rcpp::List::create( Rcpp::Named("ab_ith") = alpha_ab_ith(0),
                             Rcpp::Named("gamma_ijt") = gamma_ijt,
                             Rcpp::Named("alpha_ab_ith") = alpha_ab_ith );
}


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_shared_link_nGP_cpp( const arma::cube ab_ith,
                                                 arma::field<arma::cube> alpha_ab_ith,
                                                 const arma::field<arma::cube> y_ijtk,
                                                 const arma::field<arma::cube> w_ijtk,
                                                 arma::field<arma::cube> gamma_ijtk,
                                                 const arma::field<arma::cube> nGP_G_t,
                                                 const arma::field<arma::cube> nGP_H_t,
                                                 const arma::field<arma::cube> nGP_Wchol_t,
                                                 const bool directed=false ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  unsigned int k=0;
  unsigned int h=0;
  unsigned int alpha_i=0;
  arma::uvec aux_uvec_1;
  arma::colvec aux_vec;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::cube aux_cube_1;
  
  arma::cube y_ijt = y_ijtk(0);
  arma::cube w_ijt = w_ijtk(0);
  arma::cube gamma_ijt = gamma_ijtk(0);
  
  // Rcpp::Rcout << "(0)" << std::endl;
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int K_net = y_ijtk.n_rows;
  unsigned int H_dim = ab_ith.n_slices;
  
  // initialising objects
  arma::mat Y_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::mat PG_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::mat linpred_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::mat C_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::mat Z_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::cube X_t = arma::ones<arma::cube>(1,1,T_net); 
  
  arma::mat alpha_t = arma::zeros<arma::mat>(3,T_net);
  
  arma::cube Z_t_cov_chol = arma::zeros<arma::cube>(V_net-1,V_net-1,T_net);
  
  arma::mat dd = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::mat cc = arma::zeros<arma::mat>(3,T_net);
  arma::colvec a1 = arma::zeros<arma::colvec>(3);
  arma::mat P1chol = arma::eye<arma::mat>(3,3); P1chol.diag().fill(10);
  
  // Repeat simulation for each agent
  for( i=0; i<V_net; i++ ) {
    Y_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
    PG_t = Y_t;
    linpred_t = Y_t;
    C_t = Y_t;
    Z_t = Y_t;
    
    for( k=0; k<K_net; k++ ) {
      // Network data
      y_ijt = y_ijtk(k);
      aux_mat_1 = y_ijt.subcube(i,0,0, i,i,T_net-1);
      aux_mat_2 = y_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
      aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
      aux_mat_1.shed_rows(i,i+1);
      Y_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1) = aux_mat_1;
      
      // Polya-Gamma data
      w_ijt = w_ijtk(k);
      aux_mat_1 = w_ijt.subcube(i,0,0, i,i,T_net-1);
      aux_mat_2 = w_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
      aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
      aux_mat_1.shed_rows(i,i+1);
      PG_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1) = aux_mat_1;
      
      // linear predictor
      gamma_ijt = gamma_ijtk(k);
      aux_mat_1 = gamma_ijt.subcube(i,0,0, i,i,T_net-1);
      aux_mat_2 = gamma_ijt.subcube(i,i,0, V_net-1,i,T_net-1);
      aux_mat_1.insert_rows(aux_mat_1.n_rows,aux_mat_2);
      aux_mat_1.shed_rows(i,i+1);
      linpred_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1) = aux_mat_1;
    }
    
    // Repeat simulation for each agent and each dimension in its coordinates
    for( h=0; h<H_dim; h++ ) {
      alpha_t = arma::zeros<arma::mat>(3,T_net);
      alpha_t.row(0) = alpha_ab_ith(0).slice(h).row(i);
      // Model matrix
      // Has three columns for: U, U', A, only the first has effect on Y
      X_t = arma::zeros<arma::cube>((V_net-1)*K_net,3,T_net);
      aux_mat_1 = alpha_ab_ith(0).slice(h); aux_mat_1.shed_row(i);
      X_t.col(0) = repmat(aux_mat_1,K_net,1);
      
      // Constant term for theta in the linear predictor
      C_t = linpred_t; C_t.zeros();
      for( t=0; t<T_net; t++ ) {
        // Rcpp::Rcout << "t=" << t << std::endl;
        C_t.col(t) = linpred_t.col(t) - ( X_t.slice(t) * alpha_t.col(t) );
      }
      
      // Adjusted response
      Z_t = (Y_t-0.5)/PG_t - C_t;
      Z_t_cov_chol = arma::zeros<arma::cube>(Z_t.n_rows,Z_t.n_rows,T_net);
      for( t=0; t<T_net; t++ ) {
        Z_t_cov_chol.slice(t).diag() = pow( PG_t.col(t) , -1);
      }
      
      // Sampling coefficient using State-Space model and simulation smoother
      alpha_t = kfsim_cpp( Z_t,
                           
                           dd,
                           X_t,
                           Z_t_cov_chol,
                           
                           cc,
                           nGP_G_t(i),
                           nGP_H_t(i),
                           nGP_Wchol_t(i),
                           
                           a1,
                           P1chol );
      
      // store updated values of alpha_ab_ith
      for(alpha_i=0; alpha_i<3; alpha_i++ ) {
        alpha_ab_ith(alpha_i).slice(h).row(i) = alpha_t.row(alpha_i);
      }
      
      // Recalculate linpred with the new values of alpha_t
      for( t=0; t<T_net; t++ ) {
        linpred_t.col(t) = C_t.col(t) + ( X_t.slice(t) * alpha_t.col(t) );
      }
      
      // Update gamma_ijt with the new values of linpred
      for( k=0; k<K_net; k++ ) {
        aux_mat_1 = linpred_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1);
        aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
        
        gamma_ijt = gamma_ijtk(k);
        gamma_ijt.subcube(i,0,0, i,i,T_net-1) = aux_mat_1.rows(0,i);
        gamma_ijt.subcube(i,i,0, V_net-1,i,T_net-1) = aux_mat_1.rows(i,V_net-1);
        gamma_ijtk(k)=gamma_ijt;
      }
    }
  }
  
  return Rcpp::List::create( Rcpp::Named("ab_ith") = alpha_ab_ith(0),
                             Rcpp::Named("gamma_ijtk") = gamma_ijtk,
                             Rcpp::Named("alpha_ab_ith") = alpha_ab_ith );
}


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_link_dir_nGP_cpp( const arma::cube ab_ith_send,
                                              const arma::cube ab_ith_receive,
                                              arma::field<arma::cube> alpha_ab_ith_send,
                                              arma::field<arma::cube> alpha_ab_ith_receive,
                                              const arma::cube y_ijt,
                                              const arma::cube w_ijt,
                                              arma::cube gamma_ijt,
                                              const arma::field<arma::cube> nGP_G_t,
                                              const arma::field<arma::cube> nGP_H_t,
                                              const arma::field<arma::cube> nGP_Wchol_t,
                                              const bool directed=false )  {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  unsigned int h=0;
  unsigned int alpha_i=0;
  unsigned int dir=0;
  arma::uvec aux_uvec_1;
  arma::colvec aux_vec;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::cube aux_cube_1;
  
  // Rcpp::Rcout << "(0)" << std::endl;
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int H_dim = ab_ith_send.n_slices;
  
  // initialising objects
  arma::mat Y_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat PG_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat linpred_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat C_t = arma::zeros<arma::mat>(1,T_net);
  arma::mat Z_t = arma::zeros<arma::mat>(V_net-1,T_net);
  arma::cube X_t = arma::ones<arma::cube>(1,1,T_net); 
  
  arma::mat alpha_t = arma::zeros<arma::mat>(3,T_net);
  
  arma::cube Z_t_cov_chol = arma::zeros<arma::cube>(V_net-1,V_net-1,T_net);
  
  arma::mat dd = arma::zeros<arma::mat>(V_net-1,T_net);
  arma::mat cc = arma::zeros<arma::mat>(3,T_net);
  arma::colvec a1 = arma::zeros<arma::colvec>(3);
  arma::mat P1chol = arma::eye<arma::mat>(3,3); P1chol.diag().fill(10);
  
  // Repeat simulation for each agent
  for( dir=0; dir<2; dir++ ) {
    for( i=0; i<V_net; i++ ) {
      
      if( dir==0 ) {
        aux_mat_1 = y_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
        aux_mat_1.shed_row(i);
        Y_t = aux_mat_1;
        
        aux_mat_1 = w_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
        aux_mat_1.shed_row(i);
        PG_t = aux_mat_1;
        
        aux_mat_1 = gamma_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
        aux_mat_1.shed_row(i);
        linpred_t = aux_mat_1;
        
      } else if( dir==1 ) {
        
        aux_mat_1 = y_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_1.shed_row(i);
        Y_t = aux_mat_1;
        
        aux_mat_1 = w_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_1.shed_row(i);
        PG_t = aux_mat_1;
        
        aux_mat_1 = gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
        aux_mat_1.shed_row(i);
        linpred_t = aux_mat_1;
        
      } else {
        throw std::range_error("direction not supported");
      }
      
      // Repeat simulation for each agent and each dimension in its coordinates
      for( h=0; h<H_dim; h++ ) {
        alpha_t = arma::zeros<arma::mat>(3,T_net);
        if( dir==0 ) {
          alpha_t.row(0) = alpha_ab_ith_send(0).slice(h).row(i);
        } else if( dir==1 ) {
          alpha_t.row(0) = alpha_ab_ith_receive(0).slice(h).row(i);
        }
        // Model matrix
        // Has three columns for: U, U', A, only the first has effect on Y
        X_t = arma::zeros<arma::cube>(V_net,3,T_net);
        if( dir==0 ) {
          X_t.col(0) = alpha_ab_ith_receive(0).slice(h);
        } else if( dir==1 ) {
          X_t.col(0) = alpha_ab_ith_send(0).slice(h);
        }
        X_t.shed_row(i);
        
        // Constant term for theta in the linear predictor
        C_t = linpred_t; C_t.zeros();
        for( t=0; t<T_net; t++ ) {
          // Rcpp::Rcout << "t=" << t << std::endl;
          C_t.col(t) = linpred_t.col(t) - ( X_t.slice(t) * alpha_t.col(t) );
        }
        
        // Adjusted response
        Z_t = (Y_t-0.5)/PG_t - C_t;
        Z_t_cov_chol = arma::zeros<arma::cube>(Z_t.n_rows,Z_t.n_rows,T_net);
        for( t=0; t<T_net; t++ ) {
          Z_t_cov_chol.slice(t).diag() = pow( PG_t.col(t) , -1);
        }
        
        // Sampling coefficient using State-Space model and simulation smoother
        alpha_t = kfsim_cpp( Z_t,
                             
                             dd,
                             X_t,
                             Z_t_cov_chol,
                             
                             cc,
                             nGP_G_t(i),
                             nGP_H_t(i),
                             nGP_Wchol_t(i),
                             
                             a1,
                             P1chol );
        
        // store updated values of alpha_ab_ith
        for(alpha_i=0; alpha_i<3; alpha_i++ ) {
          if( dir==0 ) {
            alpha_ab_ith_send(alpha_i).slice(h).row(i) = alpha_t.row(alpha_i);
          } else if( dir==1 ) {
            alpha_ab_ith_receive(alpha_i).slice(h).row(i) = alpha_t.row(alpha_i);
          }
        }
        
        // Recalculate linpred with the new values of alpha_t
        for( t=0; t<T_net; t++ ) {
          linpred_t.col(t) = C_t.col(t) + ( X_t.slice(t) * alpha_t.col(t) );
        }
        
        // Update gamma_ijt with the new values of linpred
        if( dir==0 ) {
          aux_mat_1 = linpred_t;
          aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
          gamma_ijt.subcube(i,0,0, i,V_net-1,T_net-1) = aux_mat_1;
          
        } else if( dir==1 ) {
          aux_mat_1 = linpred_t;
          aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
          gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1) = aux_mat_1;
        } else {
          throw std::range_error("direction not supported");
        }
        
      }
    }
  }
  
  return Rcpp::List::create( Rcpp::Named("ab_ith_send") = alpha_ab_ith_send(0),
                             Rcpp::Named("ab_ith_receive") = alpha_ab_ith_receive(0),
                             Rcpp::Named("gamma_ijt") = gamma_ijt,
                             Rcpp::Named("alpha_ab_ith_send") = alpha_ab_ith_send,
                             Rcpp::Named("alpha_ab_ith_receive") = alpha_ab_ith_receive );
}


// [[Rcpp::export]]
Rcpp::List sample_coord_ith_shared_link_dir_nGP_cpp( const arma::cube ab_ith_send,
                                                     const arma::cube ab_ith_receive,
                                                     arma::field<arma::cube> alpha_ab_ith_send,
                                                     arma::field<arma::cube> alpha_ab_ith_receive,
                                                     const arma::field<arma::cube> y_ijtk,
                                                     const arma::field<arma::cube> w_ijtk,
                                                     arma::field<arma::cube> gamma_ijtk,
                                                     const arma::field<arma::cube> nGP_G_t,
                                                     const arma::field<arma::cube> nGP_H_t,
                                                     const arma::field<arma::cube> nGP_Wchol_t,
                                                     const bool directed=false ) {
  
  // Auxiliar objects
  unsigned int i=0;
  unsigned int t=0;
  unsigned int k=0;
  unsigned int h=0;
  unsigned int alpha_i=0;
  unsigned int dir=0;
  arma::uvec aux_uvec_1;
  arma::colvec aux_vec;
  arma::mat aux_mat_1;
  arma::mat aux_mat_2;
  arma::cube aux_cube_1;
  
  arma::cube y_ijt = y_ijtk(0);
  arma::cube w_ijt = w_ijtk(0);
  arma::cube gamma_ijt = gamma_ijtk(0);
  
  // Rcpp::Rcout << "(0)" << std::endl;
  // Network and latent space dimensions
  unsigned int V_net = y_ijt.n_rows;
  unsigned int T_net = y_ijt.n_slices;
  unsigned int K_net = y_ijtk.n_rows;
  unsigned int H_dim = ab_ith_send.n_slices;
  
  // initialising objects
  arma::mat Y_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::mat PG_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::mat linpred_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::mat C_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::mat Z_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::cube X_t = arma::ones<arma::cube>(1,1,T_net); 
  
  arma::mat alpha_t = arma::zeros<arma::mat>(3,T_net);
  
  arma::cube Z_t_cov_chol = arma::zeros<arma::cube>(V_net-1,V_net-1,T_net);
  
  arma::mat dd = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
  arma::mat cc = arma::zeros<arma::mat>(3,T_net);
  arma::colvec a1 = arma::zeros<arma::colvec>(3);
  arma::mat P1chol = arma::eye<arma::mat>(3,3); P1chol.diag().fill(10);
  
  // Repeat simulation for each direction
  for( dir=0; dir<2; dir++ ) {
    // Repeat simulation for each agent
    for( i=0; i<V_net; i++ ) {
      
      Y_t = arma::zeros<arma::mat>((V_net-1)*K_net,T_net);
      PG_t = Y_t;
      linpred_t = Y_t;
      C_t = Y_t;
      Z_t = Y_t;
      
      if( dir==0 ) {
        // Updating Y, PG, linpred as vectors, lower triangular adjacency
        for( k=0; k<K_net; k++ ) {
          y_ijt = y_ijtk(k);
          aux_mat_1 = y_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
          aux_mat_1.shed_row(i);
          Y_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1) = aux_mat_1;
          
          w_ijt = w_ijtk(k);
          aux_mat_1 = w_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
          aux_mat_1.shed_row(i);
          PG_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1) = aux_mat_1;
          
          gamma_ijt = gamma_ijtk(k);
          aux_mat_1 = gamma_ijt.subcube(i,0,0, i,V_net-1,T_net-1);
          aux_mat_1.shed_row(i);
          linpred_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1) = aux_mat_1;
        }
        
      } else if( dir==1 ) {
        // Updating Y, PG, linpred as vectors, upper triangular adjacency
        for( k=0; k<K_net; k++ ) {
          y_ijt = y_ijtk(k);
          aux_mat_1 = y_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
          aux_mat_1.shed_row(i);
          Y_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1) = aux_mat_1;
          
          w_ijt = w_ijtk(k);
          aux_mat_1 = w_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
          aux_mat_1.shed_row(i);
          PG_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1) = aux_mat_1;
          
          gamma_ijt = gamma_ijtk(k);
          aux_mat_1 = gamma_ijt.subcube(0,i,0, V_net-1,i,T_net-1);
          aux_mat_1.shed_row(i);
          linpred_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1) = aux_mat_1;
        }
        
      } else {
        throw std::range_error("direction not supported");
      }
      
      // Repeat simulation for each agent and each dimension in its coordinates
      for( h=0; h<H_dim; h++ ) {
        alpha_t = arma::zeros<arma::mat>(3,T_net);
        if( dir==0 ) {
          alpha_t.row(0) = alpha_ab_ith_send(0).slice(h).row(i);
        } else if( dir==1 ) {
          alpha_t.row(0) = alpha_ab_ith_receive(0).slice(h).row(i);
        }
        
        // Model matrix
        // Has three columns for: U, U', A, only the first has effect on Y
        X_t = arma::zeros<arma::cube>((V_net-1)*K_net,3,T_net);
        if( dir==0 ) {
          aux_mat_1 = alpha_ab_ith_receive(0).slice(h); aux_mat_1.shed_row(i);
        } else if( dir==1 ) {
          aux_mat_1 = alpha_ab_ith_send(0).slice(h); aux_mat_1.shed_row(i);
        }
        X_t.col(0) = repmat(aux_mat_1,K_net,1);
        
        // Constant term for theta in the linear predictor
        C_t = linpred_t; C_t.zeros();
        for( t=0; t<T_net; t++ ) {
          // Rcpp::Rcout << "t=" << t << std::endl;
          C_t.col(t) = linpred_t.col(t) - ( X_t.slice(t) * alpha_t.col(t) );
        }
        
        // Adjusted response
        Z_t = (Y_t-0.5)/PG_t - C_t;
        Z_t_cov_chol = arma::zeros<arma::cube>(Z_t.n_rows,Z_t.n_rows,T_net);
        for( t=0; t<T_net; t++ ) {
          Z_t_cov_chol.slice(t).diag() = pow( PG_t.col(t) , -1);
        }
        
        // Sampling coefficient using State-Space model and simulation smoother
        alpha_t = kfsim_cpp( Z_t,
                             
                             dd,
                             X_t,
                             Z_t_cov_chol,
                             
                             cc,
                             nGP_G_t(i),
                             nGP_H_t(i),
                             nGP_Wchol_t(i),
                             
                             a1,
                             P1chol );
        
        // store updated values of alpha_ab_ith
        for(alpha_i=0; alpha_i<3; alpha_i++ ) {
          if( dir==0 ) {
            alpha_ab_ith_send(alpha_i).slice(h).row(i) = alpha_t.row(alpha_i);
          } else if( dir==1 ) {
            alpha_ab_ith_receive(alpha_i).slice(h).row(i) = alpha_t.row(alpha_i);
          }
        }
        
        // Recalculate linpred with the new values of alpha_t
        for( t=0; t<T_net; t++ ) {
          linpred_t.col(t) = C_t.col(t) + ( X_t.slice(t) * alpha_t.col(t) );
        }
        
        // Update gamma_ijt with the new values of linpred
        if( dir==0 ) {
          for( k=0; k<K_net; k++ ) {
            aux_mat_1 = linpred_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1);
            aux_mat_1.insert_rows( i, arma::zeros<arma::mat>( 1, T_net) );
            
            gamma_ijt = gamma_ijtk(k);
            gamma_ijt.subcube(i,0,0, i,V_net-1,T_net-1) = aux_mat_1;
            gamma_ijtk(k)= gamma_ijt;
          }
        } else if( dir==1 ) {
          for( k=0; k<K_net; k++ ) {
            aux_mat_1 = linpred_t.rows((V_net-1)*k, (V_net-1)*(k+1)-1);
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
  }
  return Rcpp::List::create( Rcpp::Named("ab_ith_send") = alpha_ab_ith_send(0),
                             Rcpp::Named("ab_ith_receive") = alpha_ab_ith_receive(0),
                             Rcpp::Named("gamma_ijtk") = gamma_ijtk,
                             Rcpp::Named("alpha_ab_ith_send") = alpha_ab_ith_send,
                             Rcpp::Named("alpha_ab_ith_receive") = alpha_ab_ith_receive );
}
