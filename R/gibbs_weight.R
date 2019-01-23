
#' @import foreach
#' @keywords internal
sample_baseline_tk_weight <- function( theta_tk,
                                       y_ijtk, mu_ijtk,
                                       sigma_k,
                                       theta_t_cov_prior_inv,
                                       directed=FALSE ) {
  ### Sample theta_t from its conditional N-variate Gaussian posterior ###
  # V_net <- dim(y_ijtk)[1]
  # T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  for( k in 1:K_net ) { # k<-1
    # theta_tk[,k] <- sample_theta_t_DynMultiNet_bin_v2_cpp( theta_t=theta_tk[,k,drop=F],
    out_aux <- sample_baseline_tk_weight_cpp( theta_t=theta_tk[,k,drop=F],
                                              theta_t_cov_prior_inv=theta_t_cov_prior_inv,
                                              y_ijt=y_ijtk[,,,k],
                                              mu_ijt=mu_ijtk[,,,k],
                                              sigma_k=sigma_k[k],
                                              directed=directed )
    theta_tk[,k] <- out_aux$theta_t
    mu_ijtk[,,,k] <- out_aux$mu_ijt
  }
  
  return( list( theta_tk=theta_tk,
                mu_ijtk=mu_ijtk ) );
}


#' @keywords internal
sample_add_eff_itk_weight <- function( sp_itk,
                                       sp_t_cov_prior_inv,
                                       y_ijtk, mu_ijtk,
                                       sigma_k,
                                       directed=FALSE ) {
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  # This function only deals with positive weighted edges
  y_ijtk[y_ijtk<0] <- NA
  
  ### Sample sp_it ###
  for( k in 1:K_net ) { # k<-1
    if(!directed) {
      sp_it = c(sp_itk[,,k])
    } else {
      sp_it = c(sp_itk[,,k,])
    }
    out_aux <- sample_add_eff_it_weight_cpp( sp_it=sp_it, # transformed to vector
                                             sp_t_cov_prior_inv=sp_t_cov_prior_inv,
                                             y_ijt=y_ijtk[,,,k],
                                             mu_ijt=mu_ijtk[,,,k],
                                             sigma_k=sigma_k[k],
                                             directed=directed )
    if(!directed) {
      sp_itk[,,k] <- array(out_aux$sp_it,dim=c(V_net,T_net))
    } else {
      sp_itk[,,k,] <- array(out_aux$sp_it,dim=c(V_net,T_net,2))
    }
    mu_ijtk[,,,k] <- out_aux$mu_ijt
  }
  
  return( list( sp_itk=sp_itk,
                mu_ijtk=mu_ijtk ) );
}


#' @keywords internal
sample_add_eff_it_shared_weight <- function( sp_it_shared,
                                             sp_t_cov_prior_inv,
                                             y_ijtk, mu_ijtk,
                                             sigma_k,
                                             directed=FALSE ) {
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  # This function only deals with positive weighted edges
  y_ijtk[y_ijtk<0] <- NA
  
  # transform arrays to list, as armadillo fields are required as input
  y_ijtk_list <- mu_ijtk_list <- list(NULL)
  for(k in 1:K_net) {
    y_ijtk_list[[k]] <- y_ijtk[,,,k]
    mu_ijtk_list[[k]] <- mu_ijtk[,,,k]
  }
  
  ### Sample sp_it_shared ###
  out_aux <- sample_add_eff_it_shared_weight_cpp( sp_it=c(sp_it_shared), # transformed to vector
                                                  sp_t_cov_prior_inv=sp_t_cov_prior_inv,
                                                  y_ijt=y_ijtk_list,
                                                  mu_ijt=mu_ijtk_list,
                                                  sigma_k=sigma_k,
                                                  directed=directed )
  if(!directed) {
    sp_it_shared <- array(out_aux$sp_it,dim=c(V_net,T_net))
  } else {
    sp_it_shared <- array(out_aux$sp_it,dim=c(V_net,T_net,2))
  }
  for(k in 1:K_net) {mu_ijtk[,,,k] <- out_aux$mu_ijtk[k,1][[1]]}
  
  return( list( sp_it_shared=sp_it_shared,
                mu_ijtk=mu_ijtk ) );
}


#' @keywords internal
sample_coord_ith_shared_weight <- function( uv_ith_shared,
                                            uv_t_sigma_prior_inv,
                                            tau_h,
                                            y_ijtk, mu_ijtk,
                                            sigma_k,
                                            directed=FALSE ) {
  
  ### For each unit, block-sample the set of time-varying latent coordinates uv_ith ###
  V_net <- dim(uv_ith_shared)[1]
  T_net <- dim(uv_ith_shared)[2]
  K_net <- dim(y_ijtk)[4]
  H_dim <- dim(uv_ith_shared)[3]
  
  y_ijtk_list <- mu_ijtk_list <- list(NULL)
  
  for(k in 1:K_net) {
    y_ijtk_list[[k]] <- y_ijtk[,,,k]
    mu_ijtk_list[[k]] <- mu_ijtk[,,,k]
  }
  
  if( directed ) {
    out_aux <- sample_coord_ith_shared_weight_dir_cpp( u_ith_shared = uv_ith_shared[[1]],
                                                       v_ith_shared = uv_ith_shared[[2]],
                                                       uv_t_sigma_prior_inv = uv_t_sigma_prior_inv,
                                                       tau_h_shared_send = tau_h[[1]],
                                                       tau_h_shared_receive = tau_h[[2]],
                                                       y_ijtk = y_ijtk_list,
                                                       mu_ijtk = mu_ijtk_list,
                                                       sigma_k=sigma_k )
    
    uv_ith_shared[[1]] <- out_aux$u_ith_shared
    uv_ith_shared[[2]] <- out_aux$v_ith_shared
    for(k in 1:K_net) {mu_ijtk[,,,k] <- out_aux$mu_ijtk[k,1][[1]]}
  } else {
    out_aux <- sample_coord_ith_shared_weight_cpp( uv_ith_shared = uv_ith_shared,
                                                   uv_t_sigma_prior_inv = uv_t_sigma_prior_inv,
                                                   tau_h = tau_h,
                                                   y_ijtk = y_ijtk_list,
                                                   mu_ijtk = mu_ijtk_list,
                                                   sigma_k=sigma_k )
    uv_ith_shared <- out_aux$uv_ith_shared
    for(k in 1:K_net) {mu_ijtk[,,,k] <- out_aux$mu_ijtk[k,1][[1]]}
  }
  
  return( list(uv_ith_shared=uv_ith_shared,
               mu_ijtk=mu_ijtk ) )
}


#' @keywords internal
#' @importFrom abind abind
sample_coord_ithk_weight <- function( uv_ithk,
                                      uv_t_sigma_prior_inv,
                                      tau_h,
                                      y_ijtk, mu_ijtk,
                                      sigma_k,
                                      directed=FALSE ) {
  
  ### For each unit, block-sample the set of time-varying latent coordinates x_ith ###
  
  K_net <- dim(y_ijtk)[4]
  if( directed ) {
    for(k in 1:K_net){ # k<-1
      out_aux <- sample_coord_ith_weight_dir_cpp( u_ith = uv_ithk[[1]][,,,k],
                                                  v_ith = uv_ithk[[2]][,,,k],
                                                  uv_t_sigma_prior_inv = uv_t_sigma_prior_inv,
                                                  tau_h_send = tau_h[[1]][,k],
                                                  tau_h_receive = tau_h[[2]][,k],
                                                  y_ijt = y_ijtk[,,,k],
                                                  mu_ijt = mu_ijtk[,,,k],
                                                  sigma_k=sigma_k[k])
      uv_ithk[[1]][,,,k] <- out_aux$u_ith
      uv_ithk[[2]][,,,k] <- out_aux$v_ith
      mu_ijtk[,,,k] <- out_aux$mu_ijt
    }
  } else {
    for(k in 1:K_net){ # k<-1
      out_aux <- sample_coord_ith_weight_cpp( uv_ith = uv_ithk[,,,k],
                                              uv_t_sigma_prior_inv = uv_t_sigma_prior_inv,
                                              tau_h = tau_h[,k],
                                              y_ijt = y_ijtk[,,,k],
                                              mu_ijt = mu_ijtk[,,,k],
                                              sigma_k=sigma_k[k])
      uv_ithk[,,,k] <- out_aux$uv_ith
      mu_ijtk[,,,k] <- out_aux$mu_ijt
    }
    
  }
  return( list( uv_ithk=uv_ithk,
                mu_ijtk=mu_ijtk) )
}


#' @keywords internal
sample_coeff_tp_weight <- function( beta_tp,
                                    beta_t_cov_prior_inv,
                                    y_ijtk=y_ijtk, mu_ijtk=mu_ijtk,
                                    sigma_k=sigma_k,
                                    x_ijtkp_mat=x_ijtkp_mat,
                                    directed=FALSE ){
  
  # transform arrays to list, as armadillo fields are required as input
  K_net <- dim(y_ijtk)[4]
  y_ijtk_list <- mu_ijtk_list <- list(NULL)
  for(k in 1:K_net) {
    y_ijtk_list[[k]] <- y_ijtk[,,,k]
    mu_ijtk_list[[k]] <- mu_ijtk[,,,k]
  }
  
  out_aux <- sample_coeff_tp_weight_cpp( beta_tp=beta_tp,
                                         beta_t_cov_prior_inv=beta_t_cov_prior_inv,
                                         y_ijtk=y_ijtk_list, mu_ijtk=mu_ijtk_list,
                                         sigma_k=sigma_k,
                                         x_ijtkp_mat=x_ijtkp_mat,
                                         directed=directed )
  
  for(k in 1:K_net) {mu_ijtk[,,,k] <- out_aux$mu_ijtk[k,1][[1]]}
  
  return( list( beta_tp=out_aux$beta_tp,
                mu_ijtk=mu_ijtk ) );
  
  return(out_aux)
}


#' @keywords internal
sample_var_weight <- function( sigma_k,
                               sigma_k_prop_int,
                               y_ijtk, mu_ijtk,
                               directed=FALSE ) {
  
  K_net <- dim(y_ijtk)[4]
  
  for(k in 1:K_net){ # k<-1
    sigma_k[k] <- sample_var_weight_cpp( sigma_k=sigma_k[k],
                                         sigma_k_prop_int=sigma_k_prop_int[k],
                                         y_ijt=y_ijtk[,,,k], mu_ijt=mu_ijtk[,,,k],
                                         directed=TRUE )
  }
  
  return( sigma_k )
}

