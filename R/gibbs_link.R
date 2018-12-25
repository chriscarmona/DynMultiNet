
#' @import BayesLogit
#' @keywords internal
sample_pg_w_ijtk_link <- function( w_ijtk,
                                   gamma_ijtk,
                                   directed=FALSE ) {
  
  ### Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
  V_net <- dim(w_ijtk)[1]
  T_net <- dim(w_ijtk)[3]
  K_net <- dim(w_ijtk)[4]
  
  if(directed){
    idx_tmp <- matrix(TRUE,V_net,V_net); diag(idx_tmp)<-FALSE
    idx_tmp <- array(idx_tmp,dim=dim(gamma_ijtk))
    # s_aux <- gamma_ijtk[idx_tmp]
    # n_sim <- K_net*T_net*V_net*(V_net-1)
  } else {
    idx_tmp <- array(lower.tri(gamma_ijtk[,,1,1]),dim=dim(gamma_ijtk))
    # s_aux <- gamma_ijtk[idx_tmp]
    # n_sim <- K_net*T_net*V_net*(V_net-1)/2
  }
  # if( length(s_aux)!=n_sim ){ stop("There was an error sampling w_ijtk") }
  
  w_ijtk[idx_tmp] <- BayesLogit::rpg.devroye( num=sum(idx_tmp), n=1, z=gamma_ijtk[idx_tmp] )
  w_ijtk[!idx_tmp] <- NA
  
  return(w_ijtk)
  
}


#' @import foreach
#' @keywords internal
sample_baseline_tk_link <- function( eta_tk,
                                     y_ijtk, w_ijtk, gamma_ijtk,
                                     eta_t_cov_prior_inv,
                                     directed=FALSE ) {
  
  # This function only deals with binary edges (non-weighted)
  y_ijtk[y_ijtk>0] <- 1
  y_ijtk[y_ijtk<0] <- NA
  
  ### Sample eta_t from its conditional N-variate Gaussian posterior ###
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  for( k in 1:K_net ) { # k<-1
    out_aux <- sample_baseline_t_link_cpp( eta_t=eta_tk[,k,drop=F],
                                           eta_t_cov_prior_inv=eta_t_cov_prior_inv,
                                           y_ijt=y_ijtk[,,,k],
                                           w_ijt=w_ijtk[,,,k],
                                           gamma_ijt=gamma_ijtk[,,,k],
                                           directed=directed )
    eta_tk[,k] <- out_aux$eta_t
    gamma_ijtk[,,,k] <- out_aux$gamma_ijt
  }
  
  return( list( eta_tk=eta_tk,
                gamma_ijtk=gamma_ijtk ) );
}


#' @keywords internal
sample_coord_ith_shared_link <- function( ab_ith_shared,
                                          ab_t_sigma_prior_inv, tau_h,
                                          y_ijtk, w_ijtk, gamma_ijtk,
                                          directed=FALSE ) {
  
  # This function only deals with binary edges (non-weighted)
  y_ijtk[y_ijtk>0] <- 1
  y_ijtk[y_ijtk<0] <- NA
  
  ### For each unit, block-sample the set of time-varying latent coordinates ab_ith ###
  V_net <- dim(ab_ith_shared)[1]
  T_net <- dim(ab_ith_shared)[2]
  K_net <- dim(y_ijtk)[4]
  H_dim <- dim(ab_ith_shared)[3]
  
  y_ijtk_list <- w_ijtk_list <- gamma_ijtk_list <- list(NULL)
  
  for(k in 1:K_net) {
    y_ijtk_list[[k]] <- y_ijtk[,,,k]
    w_ijtk_list[[k]] <- w_ijtk[,,,k]
    gamma_ijtk_list[[k]] <- gamma_ijtk[,,,k]
  }
  
  if( directed ) {
    out_aux <- sample_coord_ith_shared_link_dir_cpp( ab_ith_shared_send = ab_ith_shared[[1]],
                                                     ab_ith_shared_receive = ab_ith_shared[[2]],
                                                     ab_t_sigma_prior_inv = ab_t_sigma_prior_inv,
                                                     tau_h_shared_send = tau_h[[1]],
                                                     tau_h_shared_receive = tau_h[[2]],
                                                     y_ijtk = y_ijtk_list,
                                                     w_ijtk = w_ijtk_list,
                                                     gamma_ijtk = gamma_ijtk_list )
    
    ab_ith_shared[[1]] <- out_aux$ab_ith_shared_send
    ab_ith_shared[[2]] <- out_aux$ab_ith_shared_receive
    for(k in 1:K_net) {gamma_ijtk[,,,k] <- out_aux$gamma_ijtk[k,1][[1]]}
  } else {
    out_aux <- sample_coord_ith_shared_link_cpp( ab_ith_shared = ab_ith_shared,
                                                 ab_t_sigma_prior_inv = ab_t_sigma_prior_inv,
                                                 tau_h = tau_h,
                                                 y_ijtk = y_ijtk_list,
                                                 w_ijtk = w_ijtk_list,
                                                 gamma_ijtk = gamma_ijtk_list )
    ab_ith_shared <- out_aux$ab_ith_shared
    for(k in 1:K_net) {gamma_ijtk[,,,k] <- out_aux$gamma_ijtk[k,1][[1]]}
  }
  
  return( list(ab_ith_shared=ab_ith_shared,
               gamma_ijtk=gamma_ijtk ) )
}


#' @keywords internal
#' @importFrom abind abind
sample_coord_ithk_link <- function( ab_ithk,
                                    ab_t_sigma_prior_inv, tau_h,
                                    y_ijtk, w_ijtk, gamma_ijtk,
                                    directed=FALSE ) {
  
  # This function only deals with binary edges (non-weighted)
  y_ijtk[y_ijtk>0] <- 1
  y_ijtk[y_ijtk<0] <- NA
  
  ### For each unit, block-sample the set of time-varying latent coordinates ab_ith ###
  
  if( directed ) {
    V_net <- dim(ab_ithk[[1]])[1]
    T_net <- dim(ab_ithk[[1]])[2]
    H_dim <- dim(ab_ithk[[1]])[3]
    K_net <- dim(ab_ithk[[1]])[4]
    
    for(k in 1:K_net){ # k<-1
      out_aux <- sample_coord_ith_link_dir_cpp( ab_ith_send = ab_ithk[[1]][,,,k],
                                                ab_ith_receive = ab_ithk[[2]][,,,k],
                                                ab_t_sigma_prior_inv = ab_t_sigma_prior_inv,
                                                tau_h_send = tau_h[[1]][,k],
                                                tau_h_receive = tau_h[[2]][,k],
                                                y_ijt = y_ijtk[,,,k],
                                                w_ijt = w_ijtk[,,,k],
                                                gamma_ijt = gamma_ijtk[,,,k] )
      ab_ithk[[1]][,,,k] <- out_aux$ab_ith_send
      ab_ithk[[2]][,,,k] <- out_aux$ab_ith_receive
      gamma_ijtk[,,,k] <- out_aux$gamma_ijt
    }
    
  } else {
    V_net <- dim(ab_ithk)[1]
    T_net <- dim(ab_ithk)[2]
    H_dim <- dim(ab_ithk)[3]
    K_net <- dim(ab_ithk)[4]
    
    for(k in 1:K_net){ # k<-1
      out_aux <- sample_coord_ith_link_cpp( ab_ith = ab_ithk[,,,k],
                                            ab_t_sigma_prior_inv = ab_t_sigma_prior_inv,
                                            tau_h = tau_h[,k],
                                            y_ijt = y_ijtk[,,,k],
                                            w_ijt = w_ijtk[,,,k],
                                            gamma_ijt = gamma_ijtk[,,,k] )
      ab_ithk[,,,k] <- out_aux$ab_ith
      gamma_ijtk[,,,k] <- out_aux$gamma_ijt
    }
    
  }
  return( list( ab_ithk=ab_ithk,
                gamma_ijtk=gamma_ijtk) )
}


#' @keywords internal
sample_v_shrink_DynMultiNet_bin <- function( v_shrink,
                                             a_1, a_2,
                                             ab_ith,
                                             ab_t_sigma_prior_inv ){
  ### Sample the global shrinkage hyperparameters from conditional gamma distributions ###
  
  V_net <- dim(ab_ith)[1]
  T_net <- dim(ab_ith)[2]
  H_dim <- nrow(v_shrink)
  
  for(h in 1:H_dim) { # h <- 3
    aux_tau_x <- vector(mode="numeric",length=H_dim)
    for( l in h:H_dim ) { # l <- 4
      if((l==1)&(h==1)){
        aux_tau <- 0 
      } else {
        aux_tau <- prod(v_shrink[setdiff(1:l,h),])
      }
      aux_x <- vector(mode="numeric",length=V_net)
      for( i in 1:V_net ) { # l <- 1
        aux_x[i] <- matrix(ab_ith[i,,l],nrow=1) %*% ab_t_sigma_prior_inv %*% matrix(ab_ith[i,,l],ncol=1)
      }
      aux_tau_x[l] <- aux_tau * sum(aux_x)
    }
    
    if(h==1){
      v_shrink[h,1] <- rgamma( n=1,
                               shape = a_1+0.5*(V_net*T_net*H_dim),
                               rate = 1+0.5*sum(aux_tau_x,na.rm=T) )
    }
    if(h>1){
      v_shrink[h,1] <- rgamma( n=1,
                               shape = a_2+0.5*(V_net*T_net*(H_dim-h+1)),
                               rate = 1+0.5*sum(aux_tau_x,na.rm=T) )
    }
  }
  return(v_shrink)
}


#' @keywords internal
sample_beta_z_edge_DynMultiNet_bin <- function( beta_z_edge,
                                                z_ijtkp, pred_id_edge, pred_all, layer_all,
                                                y_ijtk, w_ijtk, gamma_ijtk,
                                                beta_t_cov_prior_inv ) {
  ### Sample beta_z_edge from its conditional N-variate Gaussian posterior ###
  
  ### output ###
  # beta_z_edge <- matrix(0,nrow=T_net,ncol=nrow(pred_id_edge))
  
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  for( row_p in 1:nrow(pred_id_edge) ) {
    
    p <- match(pred_id_edge[row_p,1],pred_all)
    k <- match(pred_id_edge[row_p,2],layer_all)
    
    aux_sum_z2w_tk <- apply( z_ijtkp[,,,k,p]^2 * w_ijtk[,,,k], c(3), sum, na.rm=T )
    
    beta_z_edge_cov <- solve( diag(aux_sum_z2w_tk) + beta_t_cov_prior_inv )
    if(!isSymmetric(beta_z_edge_cov)) {beta_z_edge_cov[upper.tri(beta_z_edge_cov)] <- t(beta_z_edge_cov)[upper.tri(beta_z_edge_cov)]}
    
    beta_aux <- array( rep(as.numeric(beta_z_edge[,row_p]),each=V_net^2),
                       dim=c(V_net,V_net,T_net) )
    aux_vec_mean <- apply( z_ijtkp[,,,k,p]*(y_ijtk[,,,k] - 0.5 - w_ijtk[,,,k] * (gamma_ijtk[,,,k]-z_ijtkp[,,,k,p]*beta_aux)),3,sum,na.rm=T)
    aux_vec_mean <- matrix( aux_vec_mean, nrow=T_net, ncol=1 )
    
    beta_z_edge[,row_p] <- mvtnorm::rmvnorm( n=1,
                                             mean=beta_z_edge_cov %*% aux_vec_mean,
                                             sigma=beta_z_edge_cov )
  }
  
  return(beta_z_edge);
}


#' @keywords internal
sample_beta_z_layer_DynMultiNet_bin <- function( beta_z_layer,
                                                 z_tkp, pred_id_layer, pred_all, layer_all,
                                                 y_ijtk, w_ijtk, gamma_ijtk,
                                                 beta_t_cov_prior_inv,
                                                 directed = FALSE ) {
  
  # This function only deals with binary edges (non-weighted)
  y_ijtk[y_ijtk>0] <- 1
  y_ijtk[y_ijtk<0] <- NA
  
  ### Sample beta_z_layer from its conditional N-variate Gaussian posterior ###
  ### output ###
  # beta_z_layer <- matrix(0,nrow=T_net,ncol=nrow(pred_id_layer))
  
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  for( row_p in 1:nrow(pred_id_layer) ) {
    p <- match(pred_id_layer[row_p,1],pred_all)
    k <- match(pred_id_layer[row_p,2],layer_all)
    
    beta_z_layer[,row_p] <- sample_beta_z_layer_DynMultiNet_bin_cpp( beta_t=beta_z_layer[,row_p],
                                                                     z_t=z_tkp[,k,p],
                                                                     beta_t_cov_prior_inv=beta_t_cov_prior_inv,
                                                                     y_ijt=y_ijtk[,,,k],
                                                                     w_ijt=w_ijtk[,,,k],
                                                                     gamma_ijt=gamma_ijtk[,,,k],
                                                                     directed=directed )
    
  }
  return(beta_z_layer);
}

