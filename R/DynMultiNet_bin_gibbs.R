
#' @import BayesLogit
#' @keywords internal
sample_w_ijtk_DynMultiNet_bin <- function( w_ijtk,
                                           s_ijtk,
                                           directed=FALSE ) {
  
  ### Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
  V_net <- dim(w_ijtk)[1]
  T_net <- dim(w_ijtk)[3]
  K_net <- dim(w_ijtk)[4]
  
  if(directed){
    idx_tmp <- matrix(TRUE,V_net,V_net); diag(idx_tmp)<-FALSE
    idx_tmp <- array(idx_tmp,dim=dim(s_ijtk))
    # s_aux <- s_ijtk[idx_tmp]
    # n_sim <- K_net*T_net*V_net*(V_net-1)
  } else {
    idx_tmp <- array(lower.tri(s_ijtk[,,1,1]),dim=dim(s_ijtk))
    # s_aux <- s_ijtk[idx_tmp]
    # n_sim <- K_net*T_net*V_net*(V_net-1)/2
  }
  # if( length(s_aux)!=n_sim ){ stop("There was an error sampling w_ijtk") }
  
  w_ijtk[idx_tmp] <- BayesLogit::rpg.devroye( num=sum(idx_tmp), n=1, z=s_ijtk[idx_tmp] )
  w_ijtk[!idx_tmp] <- NA
  
  return(w_ijtk)
  
}


#' @import foreach
#' @keywords internal
sample_mu_tk_DynMultiNet_bin <- function( mu_tk,
                                          y_ijtk, w_ijtk, s_ijtk,
                                          mu_t_cov_prior_inv,
                                          directed=FALSE,
                                          use_cpp=TRUE,
                                          parallel_mcmc=FALSE,
                                          calc_method_R=1 ) {
  ### Sample mu_t from its conditional N-variate Gaussian posterior ###
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  if(parallel_mcmc) {
    mu_tk <- foreach( k=1:K_net, .combine=cbind, .inorder=TRUE ) %dopar% { # k<-1
      sample_mu_t_DynMultiNet_bin_cpp( mu_t=mu_tk[,k,drop=FALSE],
                                       mu_t_cov_prior_inv=mu_t_cov_prior_inv,
                                       y_ijt=y_ijtk[,,,k],
                                       w_ijt=w_ijtk[,,,k],
                                       s_ijt=s_ijtk[,,,k],
                                       directed=directed )
    }
  } else {
    for( k in 1:K_net ) { # k<-1
      # mu_tk[,k] <- sample_mu_t_DynMultiNet_bin_v2_cpp( mu_t=mu_tk[,k,drop=F],
      out_aux <- sample_mu_t_DynMultiNet_bin_v2_cpp( mu_t=mu_tk[,k,drop=F],
                                                     mu_t_cov_prior_inv=mu_t_cov_prior_inv,
                                                     y_ijt=y_ijtk[,,,k],
                                                     w_ijt=w_ijtk[,,,k],
                                                     s_ijt=s_ijtk[,,,k],
                                                     directed=directed )
      mu_tk[,k] <- out_aux$mu_t
      s_ijtk[,,,k] <- out_aux$s_ijt
    }
  }
  
  return( list( mu_tk=mu_tk,
                s_ijtk=s_ijtk ) );
}


#' @keywords internal
sample_beta_z_edge_DynMultiNet_bin <- function( beta_z_edge,
                                                z_ijtkp, pred_id_edge, pred_all, layer_all,
                                                y_ijtk, w_ijtk, s_ijtk,
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
    aux_vec_mean <- apply( z_ijtkp[,,,k,p]*(y_ijtk[,,,k] - 0.5 - w_ijtk[,,,k] * (s_ijtk[,,,k]-z_ijtkp[,,,k,p]*beta_aux)),3,sum,na.rm=T)
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
                                                 y_ijtk, w_ijtk, s_ijtk,
                                                 beta_t_cov_prior_inv,
                                                 use_cpp = TRUE,
                                                 directed = FALSE ) {
  ### Sample beta_z_layer from its conditional N-variate Gaussian posterior ###
  ### output ###
  # beta_z_layer <- matrix(0,nrow=T_net,ncol=nrow(pred_id_layer))
  
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  for( row_p in 1:nrow(pred_id_layer) ) {
    p <- match(pred_id_layer[row_p,1],pred_all)
    k <- match(pred_id_layer[row_p,2],layer_all)
    if(use_cpp) {
      beta_z_layer[,row_p] <- sample_beta_z_layer_DynMultiNet_bin_cpp( beta_t=beta_z_layer[,row_p],
                                                                       z_t=z_tkp[,k,p],
                                                                       beta_t_cov_prior_inv=beta_t_cov_prior_inv,
                                                                       y_ijt=y_ijtk[,,,k],
                                                                       w_ijt=w_ijtk[,,,k],
                                                                       s_ijt=s_ijtk[,,,k],
                                                                       directed=directed )
    } else {
      if(directed) {
        stop("beta_z_layer sampling not supported in R code for directed networks.")
      } else {
        # method 1: Linear equations
        Y<-NULL; W<-NULL; S<-NULL
        for(i in 2:V_net) {
          #i<-3
          Y <- rbind( Y,
                      matrix(c(t(matrix(y_ijtk[i,1:i,,k],nrow=i,ncol=T_net)[-i,,drop=F])),T_net*(i-1),ncol=1) )
          W <- rbind( W,
                      matrix(c(t(matrix(w_ijtk[i,1:i,,k],nrow=i,ncol=T_net)[-i,,drop=F])),T_net*(i-1),ncol=1) )
          S <- rbind( S,
                      matrix(c(t(matrix(s_ijtk[i,1:i,,k],nrow=i,ncol=T_net)[-i,,drop=F])),T_net*(i-1),ncol=1) )
        }
        W_diag <- diag(as.numeric(W))
        
        X <- kronecker( matrix(1,V_net*(V_net-1)/2,1) ,diag(z_tkp[,k,p]))
        
        C = S - X %*% matrix(beta_z_layer[,row_p],T_net,1);
        Z = (Y-0.5)/W - C;
        
        beta_cov_inv <- t(X) %*% W_diag %*%  X + beta_t_cov_prior_inv
        beta_cov<- solve( beta_cov_inv )
        aux_vec_mean <- t(X) %*% W_diag %*% Z
        
        beta_z_layer[,row_p] <- mvtnorm::rmvnorm( n=1,
                                                  mean=beta_cov %*% aux_vec_mean,
                                                  sigma=beta_cov )
      }
    }
  }
  return(beta_z_layer);
}


#' @keywords internal
sample_x_ith_shared_DynMultiNet_bin <- function( x_ith_shared,
                                                 x_t_sigma_prior_inv, tau_h,
                                                 y_ijtk, w_ijtk, s_ijtk,
                                                 directed=FALSE ) {
  
  ### For each unit, block-sample the set of time-varying latent coordinates x_ith ###
  V_net <- dim(x_ith_shared)[1]
  T_net <- dim(x_ith_shared)[2]
  K_net <- dim(y_ijtk)[4]
  H_dim <- dim(x_ith_shared)[3]
  
  y_ijtk_list <- w_ijtk_list <- s_ijtk_list <- list(NULL)
  
  for(k in 1:K_net) {
    y_ijtk_list[[k]] <- y_ijtk[,,,k]
    w_ijtk_list[[k]] <- w_ijtk[,,,k]
    s_ijtk_list[[k]] <- s_ijtk[,,,k]
  }
  
  if( directed ) {
    x_ith_shared <- sample_x_ith_shared_DynMultiNet_bin_dir_cpp( x_ith_shared_send = x_ith_shared[[1]],
                                                                 x_ith_shared_receive = x_ith_shared[[2]],
                                                                 x_t_sigma_prior_inv = x_t_sigma_prior_inv,
                                                                 tau_h_shared_send = tau_h[[1]],
                                                                 tau_h_shared_receive = tau_h[[2]],
                                                                 y_ijtk = y_ijtk_list,
                                                                 w_ijtk = w_ijtk_list,
                                                                 s_ijtk = s_ijtk_list )
  } else {
    x_ith_shared <- sample_x_ith_shared_DynMultiNet_bin_cpp( x_ith_shared = x_ith_shared,
                                                             x_t_sigma_prior_inv = x_t_sigma_prior_inv,
                                                             tau_h = tau_h,
                                                             y_ijtk = y_ijtk_list,
                                                             w_ijtk = w_ijtk_list,
                                                             s_ijtk = s_ijtk_list )
  }
  
  return( x_ith_shared )
}


#' @keywords internal
#' @importFrom abind abind
sample_x_ithk_DynMultiNet_bin <- function( x_ithk,
                                           x_t_sigma_prior_inv, tau_h,
                                           y_ijtk, w_ijtk, s_ijtk,
                                           directed=FALSE,
                                           parallel_mcmc=FALSE ) {
  
  ### For each unit, block-sample the set of time-varying latent coordinates x_ith ###
  
  if( directed ) {
    V_net <- dim(x_ithk[[1]])[1]
    T_net <- dim(x_ithk[[1]])[2]
    H_dim <- dim(x_ithk[[1]])[3]
    K_net <- dim(x_ithk[[1]])[4]
    
    if(parallel_mcmc) {
      x_ithk <- foreach( k=1:K_net,
                         .combine=function(x,y){ list( abind::abind(x[[1]],y[[1]],along=4),
                                                       abind::abind(x[[2]],y[[2]],along=4) ) }
      ) %dopar% { # k<-1
        sample_x_ith_DynMultiNet_bin_dir_cpp( x_ith_send = x_ithk[[1]][,,,k],
                                              x_ith_receive = x_ithk[[2]][,,,k],
                                              x_t_sigma_prior_inv = x_t_sigma_prior_inv,
                                              tau_h_send = tau_h[[1]][,k],
                                              tau_h_receive = tau_h[[2]][,k],
                                              y_ijt = y_ijtk[,,,k],
                                              w_ijt = w_ijtk[,,,k],
                                              s_ijt = s_ijtk[,,,k] )
      }
    } else {
      for(k in 1:K_net){ # k<-1
        x_ithk_aux <- sample_x_ith_DynMultiNet_bin_dir_cpp( x_ith_send = x_ithk[[1]][,,,k],
                                                            x_ith_receive = x_ithk[[2]][,,,k],
                                                            x_t_sigma_prior_inv = x_t_sigma_prior_inv,
                                                            tau_h_send = tau_h[[1]][,k],
                                                            tau_h_receive = tau_h[[2]][,k],
                                                            y_ijt = y_ijtk[,,,k],
                                                            w_ijt = w_ijtk[,,,k],
                                                            s_ijt = s_ijtk[,,,k] )
        x_ithk[[1]][,,,k] <- x_ithk_aux[[1]]
        x_ithk[[2]][,,,k] <- x_ithk_aux[[2]]
      }
    }
  } else {
    V_net <- dim(x_ithk)[1]
    T_net <- dim(x_ithk)[2]
    H_dim <- dim(x_ithk)[3]
    K_net <- dim(x_ithk)[4]
    
    if(parallel_mcmc) {
      x_ithk <- foreach( k=1:K_net, .combine=function(x,y){abind::abind(x,y,along=4)} ) %dopar% {
        sample_x_ith_DynMultiNet_bin_cpp( x_ith = x_ithk[,,,k],
                                          x_t_sigma_prior_inv = x_t_sigma_prior_inv,
                                          tau_h = tau_h[,k],
                                          y_ijt = y_ijtk[,,,k],
                                          w_ijt = w_ijtk[,,,k],
                                          s_ijt = s_ijtk[,,,k] )
      }
    } else {
      for(k in 1:K_net){ # k<-1
        x_ithk[,,,k] <- sample_x_ith_DynMultiNet_bin_cpp( x_ith = x_ithk[,,,k],
                                                          x_t_sigma_prior_inv = x_t_sigma_prior_inv,
                                                          tau_h = tau_h[,k],
                                                          y_ijt = y_ijtk[,,,k],
                                                          w_ijt = w_ijtk[,,,k],
                                                          s_ijt = s_ijtk[,,,k] )
      }
    }
  }
  return( x_ithk )
}


#' @keywords internal
sample_v_shrink_DynMultiNet_bin <- function( v_shrink,
                                             a_1, a_2,
                                             x_ith,
                                             x_t_sigma_prior_inv ){
  ### Sample the global shrinkage hyperparameters from conditional gamma distributions ###
  
  V_net <- dim(x_ith)[1]
  T_net <- dim(x_ith)[2]
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
        aux_x[i] <- matrix(x_ith[i,,l],nrow=1) %*% x_t_sigma_prior_inv %*% matrix(x_ith[i,,l],ncol=1)
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
