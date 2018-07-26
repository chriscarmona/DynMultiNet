
#' @export
sample_w_ijtk_DynMultiNet_bin <- function( w_ijtk, s_ijtk ) {
  
  ### Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
  
  V_net <- dim(w_ijtk)[1]
  T_net <- dim(w_ijtk)[3]
  K_net <- dim(w_ijtk)[4]
  
  s_aux <- c(s_ijtk); s_aux <- s_aux[!is.na(s_aux)]
  if(length(s_aux)!=T_net*V_net*(V_net-1)/2){ stop("There was an error sampling w_ijtk") }
  w_ijtk[!is.na(s_ijtk)] <- BayesLogit::rpg.devroye( num=T_net*V_net*(V_net-1)/2, n=1, z=s_aux )
  
  return(w_ijtk)
  
}



#' @export
sample_mu_tk_DynMultiNet_bin <- function( mu_tk,
                                          y_ijtk, w_ijtk, s_ijtk,
                                          mu_t_cov_prior_inv,
                                          use_cpp=TRUE,
                                          calc_method_R=1 ) {
  ### Sample mu_t from its conditional N-variate Gaussian posterior ###
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  if(use_cpp) {
    for( k in 1:K_net ){ # k<-1
      sample_mu_tk_DynMultiNet_bin
      mu_tk[,k] <- sample_mu_t_DynMultiNet_bin_cpp( mu_t=mu_tk[,k,drop=F],
                                                    mu_t_cov_prior_inv=mu_t_cov_prior_inv,
                                                    y_ijt=y_ijtk[,,,k],
                                                    w_ijt=w_ijtk[,,,k],
                                                    s_ijt=s_ijtk[,,,k] )
    }
  } else {
    for( k in 1:K_net ){ # k<-1
      mu_t=mu_tk[,k,drop=F]
      y_ijt=y_ijtk[,,,k]
      w_ijt=w_ijtk[,,,k]
      s_ijt=s_ijtk[,,,k]
      
      if(calc_method_R==1){
        # Option 1: [Durante, 2014]
        aux_sum_w_t <- apply(w_ijt,3,sum,na.rm=T)
        aux_mat <- solve( diag(aux_sum_w_t) + mu_t_cov_prior_inv )
        if(!isSymmetric(aux_mat)) {aux_mat[upper.tri(aux_mat)] <- t(aux_mat)[upper.tri(aux_mat)]}
        mu_t_cov <- aux_mat
        
        mu_t_aux <- array( rep(as.numeric(mu_t),each=V_net^2),
                           dim=c(V_net,V_net,T_net) )
        aux_vec_mean <- apply(y_ijt[,,] - 0.5 - w_ijt[,,] * (s_ijt[,,]-mu_t_aux),3,sum,na.rm=T)
        aux_vec_mean <- matrix(aux_vec_mean,nrow=T_net,ncol=1)
      } else if(calc_method_R==2) {
        # Option 2: direct linear model
        X <- kronecker(matrix(1,V_net*(V_net-1)/2,1),diag(T_net))
        Y<-NULL; W<-NULL; S<-NULL
        for(i in 2:V_net) {
          #i<-3
          Y <- rbind( Y,
                      matrix(c(t(matrix(y_ijt[i,1:i,],nrow=i,ncol=T_net)[-i,,drop=F])),T_net*(i-1),ncol=1) )
          W <- rbind( W,
                      matrix(c(t(matrix(w_ijt[i,1:i,],nrow=i,ncol=T_net)[-i,,drop=F])),T_net*(i-1),ncol=1) )
          S <- rbind( S,
                      matrix(c(t(matrix(s_ijt[i,1:i,],nrow=i,ncol=T_net)[-i,,drop=F])),T_net*(i-1),ncol=1) )
        }
        
        S_minus <- kronecker(matrix(1,V_net*(V_net-1)/2,1),matrix(mu_t,T_net,1))
        W_diag <- diag(as.numeric(W))
        Z <- (Y - 0.5)/W-(S-S_minus)
        mu_t_cov <- solve( t(X) %*% W_diag %*%  X + mu_t_cov_prior_inv )
        aux_vec_mean <- t(X) %*% W_diag %*% Z
      }
      
      mu_tk[,k] <- mvtnorm::rmvnorm( n=1,
                                     mean=mu_t_cov %*% aux_vec_mean,
                                     sigma=mu_t_cov )
    }
  }
  
  return(mu_tk);
}



#' @export
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



#' @export
sample_beta_z_layer_DynMultiNet_bin <- function( beta_z_layer,
                                                 z_tkp, pred_id_layer, pred_all, layer_all,
                                                 y_ijtk, w_ijtk, s_ijtk,
                                                 beta_t_cov_prior_inv,
                                                 use_cpp=TRUE ) {
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
                                                                       s_ijt=s_ijtk[,,,k] )
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
  return(beta_z_layer);
}



#' @export
sample_x_iht_mat_DynMultiNet_bin <- function( x_iht_mat,
                                              x_t_sigma_prior_inv, tau_h,
                                              y_ijtk, w_ijtk, s_ijtk ) {
  
  ### For each unit, block-sample the set of time-varying latent coordinates x_iht ###
  K_net <- dim(y_ijtk)[4]
  
  y_ijtk_list <- w_ijtk_list <- s_ijtk_list <- list(NULL)
  for(k in 1:K_net) {
    y_ijtk_list[[k]] <- y_ijtk[,,,k]
    w_ijtk_list[[k]] <- w_ijtk[,,,k]
    s_ijtk_list[[k]] <- s_ijtk[,,,k]
  }
  
  x_iht_mat <- sample_x_iht_mat_DynMultiNet_bin_cpp( x_iht_mat = x_iht_mat,
                                                     x_t_sigma_prior_inv = x_t_sigma_prior_inv,
                                                     tau_h = tau_h,
                                                     y_ijtk = y_ijtk_list,
                                                     w_ijtk = w_ijtk_list,
                                                     s_ijtk = s_ijtk_list )
  
  return( x_iht_mat )
}



#' @export
sample_v_dim_DynMultiNet_bin <- function( v_dim, a_1, a_2,
                                          x_iht,
                                          x_t_sigma_prior_inv ){
  ### Sample the global shrinkage hyperparameters from conditional gamma distributions ###
  
  V_net <- dim(x_iht)[1]
  T_net <- dim(x_iht)[3]
  H_dim <- nrow(v_dim)
  
  for(h in 2:H_dim) {
    tau_h <- matrix(cumprod(v_dim), nrow=H_dim, ncol=1 )
    tau_minush_l <- matrix(tau_h, nrow=H_dim, ncol=H_dim )
    tau_minush_l[upper.tri(tau_minush_l)] <- NA ; tau_minush_l[1,1] <- NA
    tau_minush_l <- tau_minush_l / matrix(v_dim, nrow=H_dim, ncol=H_dim, byrow=T)
    aux_1 <- apply(tau_minush_l,2,sum,na.rm=TRUE)
    aux_2 <- vector(mode="numeric",length=V_net)
    for( i in 1:V_net ){
      aux_2[i]<-matrix(x_iht[i,h,],nrow=1) %*% x_t_sigma_prior_inv %*% matrix(x_iht[i,h,],ncol=1)
    }
    if(h==1){
      v_dim[h,1] <- rgamma( n=1,
                            shape = a_1+0.5*(V_net*T_net*H_dim),
                            rate = 1+0.5*aux_1[h]*sum(aux_2) )
    }
    if(h>1){
      v_dim[h,1] <- rgamma( n=1,
                            shape = a_2+0.5*(V_net*T_net*(H_dim-h+1)),
                            rate = 1+0.5*aux_1[h]*sum(aux_2) )
    }
  }
  return(v_dim)
}
