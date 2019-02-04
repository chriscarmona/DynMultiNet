
#' @import BayesLogit
#' @keywords internal
sample_pg_w_ijtk_link <- function( w_ijtk,
                                   gamma_ijtk,
                                   directed=FALSE ) {
  
  ### Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
  V_net <- dim(w_ijtk)[1]
  K_net <- dim(w_ijtk)[4]
  
  if(directed){
    idx_tmp <- matrix(TRUE,V_net,V_net); diag(idx_tmp)<-FALSE
    idx_tmp <- array(idx_tmp,dim=dim(gamma_ijtk))
  } else {
    idx_tmp <- array(lower.tri(gamma_ijtk[,,1,1]),dim=dim(gamma_ijtk))
  }
  
  w_ijtk[idx_tmp] <- BayesLogit::rpg.devroye( num=sum(idx_tmp), n=1, z=gamma_ijtk[idx_tmp] )
  w_ijtk[!idx_tmp] <- NA
  
  return(w_ijtk)
  
}


#' @keywords internal
sample_baseline_tk_link <- function( eta_tk,
                                     y_ijtk, w_ijtk, gamma_ijtk,
                                     class_dyn=c("GP","nGP")[1],
                                     eta_t_cov_prior_inv=NULL,
                                     nGP_mat=NULL,
                                     alpha_eta_tk=NULL,
                                     directed=FALSE ) {
  
  K_net <- dim(y_ijtk)[4]
  
  # This function only deals with binary edges (non-weighted)
  y_ijtk[y_ijtk>0] <- 1
  y_ijtk[y_ijtk<0] <- NA
  
  ### Sample eta_t from its conditional N-variate Gaussian posterior ###
  for( k in 1:K_net ) { # k<-1
    if( class_dyn==c("GP","nGP")[1] ){
      if(is.null(eta_t_cov_prior_inv)){stop("eta_t_cov_prior_inv not provided!")}
      out_aux <- sample_baseline_t_link_GP_cpp( eta_t=eta_tk[,k,drop=F],
                                                eta_t_cov_prior_inv=eta_t_cov_prior_inv,
                                                y_ijt=y_ijtk[,,,k],
                                                w_ijt=w_ijtk[,,,k],
                                                gamma_ijt=gamma_ijtk[,,,k],
                                                directed=directed )
    } else if( class_dyn==c("GP","nGP")[2] ){
      if(is.null(nGP_mat)){stop("nGP_mat not provided!")}
      out_aux <- sample_baseline_t_link_nGP_cpp( eta_t=eta_tk[,k,drop=F],
                                                 alpha_eta_t=t(alpha_eta_tk[,k,]),
                                                 
                                                 y_ijt=y_ijtk[,,,k],
                                                 w_ijt=w_ijtk[,,,k],
                                                 gamma_ijt=gamma_ijtk[,,,k],
                                                 
                                                 nGP_G_t = nGP_mat$G[k,,,],
                                                 nGP_H_t = nGP_mat$H[k,,,],
                                                 nGP_Wchol_t = nGP_mat$Wchol[k,,,],
                                                 
                                                 directed=directed )
      alpha_eta_tk[,k,] <- t(out_aux$alpha_eta_t)
    }
    eta_tk[,k] <- out_aux$eta_t
    gamma_ijtk[,,,k] <- out_aux$gamma_ijt
  }
  
  return( list( eta_tk=eta_tk,
                gamma_ijtk=gamma_ijtk,
                alpha_eta_tk=alpha_eta_tk ) );
}


#' @keywords internal
sample_coord_ith_shared_link <- function( ab_ith,
                                          y_ijtk, w_ijtk, gamma_ijtk,
                                          class_dyn=c("GP","nGP")[1],
                                          ab_t_sigma_prior_inv=NULL, tau_h=NULL,
                                          nGP_mat=NULL, alpha_ab_ith=NULL,
                                          directed=FALSE ) {
  
  # This function only deals with binary edges (non-weighted)
  y_ijtk[y_ijtk>0] <- 1
  y_ijtk[y_ijtk<0] <- NA
  
  ### For each unit, block-sample the set of time-varying latent coordinates ab_ith ###
  y_ijtk_list <- w_ijtk_list <- gamma_ijtk_list <- list(NULL)
  
  K_net <- dim(y_ijtk)[4]
  if( directed ) {
    V_net <- dim(ab_ith[[1]])[1]
    T_net <- dim(ab_ith[[1]])[2]
    H_dim <- dim(ab_ith[[1]])[3]
  } else {
    V_net <- dim(ab_ith)[1]
    T_net <- dim(ab_ith)[2]
    H_dim <- dim(ab_ith)[3]
  }
  
  for(k in 1:K_net) {
    y_ijtk_list[[k]] <- y_ijtk[,,,k]
    w_ijtk_list[[k]] <- w_ijtk[,,,k]
    gamma_ijtk_list[[k]] <- gamma_ijtk[,,,k]
  }
  
  if( directed & class_dyn=="GP" ){
    out_aux <- sample_coord_ith_shared_link_dir_GP_cpp( ab_ith_send = ab_ith[[1]],
                                                        ab_ith_receive = ab_ith[[2]],
                                                        y_ijtk = y_ijtk_list,
                                                        w_ijtk = w_ijtk_list,
                                                        gamma_ijtk = gamma_ijtk_list,
                                                        ab_t_sigma_prior_inv = ab_t_sigma_prior_inv,
                                                        tau_h_shared_send = tau_h[[1]],
                                                        tau_h_shared_receive = tau_h[[2]] )
    
    ab_ith[[1]] <- out_aux$ab_ith_send
    ab_ith[[2]] <- out_aux$ab_ith_receive
    for(k in 1:K_net) {gamma_ijtk[,,,k] <- out_aux$gamma_ijtk[k,1][[1]]}
  } else if(!directed & class_dyn=="GP"){
    out_aux <- sample_coord_ith_shared_link_GP_cpp( ab_ith = ab_ith,
                                                    y_ijtk = y_ijtk_list,
                                                    w_ijtk = w_ijtk_list,
                                                    gamma_ijtk = gamma_ijtk_list,
                                                    ab_t_sigma_prior_inv = ab_t_sigma_prior_inv,
                                                    tau_h = tau_h )
    ab_ith <- out_aux$ab_ith
    for(k in 1:K_net) {gamma_ijtk[,,,k] <- out_aux$gamma_ijtk[k,1][[1]]}
  } else if(directed & class_dyn=="GP"){
    stop("Not implemented")
  } else if(!directed & class_dyn=="nGP"){
    alpha_ab_ith_list <- nGP_G_t <- nGP_H_t <- nGP_Wchol_t <- list(NULL)
    for(i in 1:V_net){
      nGP_G_t[[i]] <- nGP_mat$G[i,,,]
      nGP_H_t[[i]] <- nGP_mat$H[i,,,]
      nGP_Wchol_t[[i]] <- nGP_mat$Wchol[i,,,]
    }
    for(alpha_i in 1:3) { alpha_ab_ith_list[[alpha_i]] <- alpha_ab_ith[,,,alpha_i] }
    out_aux <- sample_coord_ith_shared_link_nGP_cpp( ab_ith = ab_ith,
                                                     alpha_ab_ith = alpha_ab_ith_list,

                                                     y_ijtk = y_ijtk_list,
                                                     w_ijtk = w_ijtk_list,
                                                     gamma_ijtk = gamma_ijtk_list,

                                                     nGP_G_t = nGP_G_t,
                                                     nGP_H_t = nGP_H_t,
                                                     nGP_Wchol_t = nGP_Wchol_t )
    ab_ith <- out_aux$ab_ith
    for(k in 1:K_net) {gamma_ijtk[,,,k] <- out_aux$gamma_ijtk[k,1][[1]]}
    for(alpha_i in 1:3) { alpha_ab_ith[,,,alpha_i] <- out_aux$alpha_ab_ith[[alpha_i]] }
  }
  
  return( list( ab_ith=ab_ith,
                gamma_ijtk=gamma_ijtk,
                alpha_ab_ith=alpha_ab_ith ) )
}


#' @keywords internal
#' @importFrom abind abind
sample_coord_ithk_link <- function( ab_ithk,
                                    y_ijtk, w_ijtk, gamma_ijtk,
                                    class_dyn=c("GP","nGP")[1],
                                    ab_t_sigma_prior_inv=NULL, tau_h=NULL,
                                    nGP_mat=NULL, alpha_ab_ithk=NULL,
                                    directed=FALSE ) {
  
  # This function only deals with binary edges (non-weighted)
  y_ijtk[y_ijtk>0] <- 1
  y_ijtk[y_ijtk<0] <- NA
  
  if( directed ) {
    V_net <- dim(ab_ithk[[1]])[1]
    T_net <- dim(ab_ithk[[1]])[2]
    H_dim <- dim(ab_ithk[[1]])[3]
    K_net <- dim(ab_ithk[[1]])[4]
  } else {
    V_net <- dim(ab_ithk)[1]
    T_net <- dim(ab_ithk)[2]
    H_dim <- dim(ab_ithk)[3]
    K_net <- dim(ab_ithk)[4]
  }
  
  ### For each unit, block-sample the set of time-varying latent coordinates ab_ith ###
  if(directed & class_dyn=="GP"){
    
    for(k in 1:K_net){ # k<-1
      out_aux <- sample_coord_ith_link_dir_GP_cpp( ab_ith_send = ab_ithk[[1]][,,,k],
                                                   ab_ith_receive = ab_ithk[[2]][,,,k],
                                                   y_ijt = y_ijtk[,,,k],
                                                   w_ijt = w_ijtk[,,,k],
                                                   gamma_ijt = gamma_ijtk[,,,k],
                                                   ab_t_sigma_prior_inv = ab_t_sigma_prior_inv,
                                                   tau_h_send = tau_h[[1]][,k],
                                                   tau_h_receive = tau_h[[2]][,k] )
      ab_ithk[[1]][,,,k] <- out_aux$ab_ith_send
      ab_ithk[[2]][,,,k] <- out_aux$ab_ith_receive
      gamma_ijtk[,,,k] <- out_aux$gamma_ijt
    }
    
  } else if( !directed & class_dyn=="GP"){
    
    for(k in 1:K_net){ # k<-1
      out_aux <- sample_coord_ith_link_GP_cpp( ab_ith = ab_ithk[,,,k],
                                               y_ijt = y_ijtk[,,,k],
                                               w_ijt = w_ijtk[,,,k],
                                               gamma_ijt = gamma_ijtk[,,,k],
                                               ab_t_sigma_prior_inv = ab_t_sigma_prior_inv,
                                               tau_h = tau_h[,k] )
      ab_ithk[,,,k] <- out_aux$ab_ith
      gamma_ijtk[,,,k] <- out_aux$gamma_ijt
    }
    
  } else if(directed & class_dyn=="nGP"){
    stop("Not implemented...")
  } else if(!directed & class_dyn=="nGP"){
    
    for(k in 1:K_net){ # k<-1
      alpha_ab_ith_list <- nGP_G_t <- nGP_H_t <- nGP_Wchol_t <- list(NULL)
      for(i in 1:V_net){
        nGP_G_t[[i]] <- nGP_mat$G[i,k,,,]
        nGP_H_t[[i]] <- nGP_mat$H[i,k,,,]
        nGP_Wchol_t[[i]] <- nGP_mat$Wchol[i,k,,,]
      }
      for(alpha_i in 1:3) { alpha_ab_ith_list[[alpha_i]] <- alpha_ab_ithk[,,,k,alpha_i] }
      
      out_aux <- sample_coord_ith_link_nGP_cpp( ab_ith = ab_ithk[,,,k],
                                                alpha_ab_ith = alpha_ab_ith_list,
                                                
                                                y_ijt = y_ijtk[,,,k],
                                                w_ijt = w_ijtk[,,,k],
                                                gamma_ijt = gamma_ijtk[,,,k],
                                                
                                                nGP_G_t = nGP_G_t,
                                                nGP_H_t = nGP_H_t,
                                                nGP_Wchol_t = nGP_Wchol_t )
      ab_ithk[,,,k] <- out_aux$ab_ith
      gamma_ijtk[,,,k] <- out_aux$gamma_ijt
      for(alpha_i in 1:3) { alpha_ab_ithk[,,,k,alpha_i] <- out_aux$alpha_ab_ith[[alpha_i]] }
    }
  }
  return( list( ab_ithk=ab_ithk,
                gamma_ijtk=gamma_ijtk,
                alpha_ab_ithk=alpha_ab_ithk) )
}


#' @keywords internal
sample_add_eff_itk_link <- function( sp_itk,
                                     sp_t_cov_prior_inv,
                                     y_ijtk, w_ijtk, gamma_ijtk,
                                     directed=FALSE ) {
  
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  # This function only deals with binary edges (non-weighted)
  y_ijtk[y_ijtk>0] <- 1
  y_ijtk[y_ijtk<0] <- NA
  
  ### Sample sp_it ###
  for( k in 1:K_net ) { # k<-1
    if(!directed) {
      sp_it = c(sp_itk[,,k])
    } else {
      sp_it = c(sp_itk[,,k,])
    }
    out_aux <- sample_add_eff_it_link_cpp( sp_it=sp_it, # transformed to vector
                                           sp_t_cov_prior_inv=sp_t_cov_prior_inv,
                                           y_ijt=y_ijtk[,,,k],
                                           w_ijt=w_ijtk[,,,k],
                                           gamma_ijt=gamma_ijtk[,,,k],
                                           directed=directed )
    if(!directed) {
      sp_itk[,,k] <- array(out_aux$sp_it,dim=c(V_net,T_net))
    } else {
      sp_itk[,,k,] <- array(out_aux$sp_it,dim=c(V_net,T_net,2))
    }
    gamma_ijtk[,,,k] <- out_aux$gamma_ijt
  }
  
  return( list( sp_itk=sp_itk,
                gamma_ijtk=gamma_ijtk ) );
}


#' @keywords internal
sample_add_eff_it_shared_link <- function( sp_it_shared,
                                           sp_t_cov_prior_inv,
                                           y_ijtk, w_ijtk, gamma_ijtk,
                                           directed=FALSE ) {
  
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  # This function only deals with binary edges (non-weighted)
  y_ijtk[y_ijtk>0] <- 1
  y_ijtk[y_ijtk<0] <- NA
  
  # transform arrays to list, as armadillo fields are required as input
  y_ijtk_list <- w_ijtk_list <- gamma_ijtk_list <- list(NULL)
  for(k in 1:K_net) {
    y_ijtk_list[[k]] <- y_ijtk[,,,k]
    w_ijtk_list[[k]] <- w_ijtk[,,,k]
    gamma_ijtk_list[[k]] <- gamma_ijtk[,,,k]
  }
  
  ### Sample sp_it_shared ###
  out_aux <- sample_add_eff_it_shared_link_cpp( sp_it=c(sp_it_shared), # transformed to vector
                                                sp_t_cov_prior_inv=sp_t_cov_prior_inv,
                                                y_ijt=y_ijtk_list,
                                                w_ijt=w_ijtk_list,
                                                gamma_ijt=gamma_ijtk_list,
                                                directed=directed )
  
  if(!directed) {
    sp_it_shared <- array(out_aux$sp_it,dim=c(V_net,T_net))
  } else {
    sp_it_shared <- array(out_aux$sp_it,dim=c(V_net,T_net,2))
  }
  for(k in 1:K_net) {gamma_ijtk[,,,k] <- out_aux$gamma_ijtk[k,1][[1]]}
  
  return( list( sp_it_shared=sp_it_shared,
                gamma_ijtk=gamma_ijtk ) );
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


#' @keywords internal
sample_coeff_tp_link <- function( beta_tp,
                                  beta_t_cov_prior_inv,
                                  y_ijtk, w_ijtk, gamma_ijtk,
                                  x_ijtkp_mat,
                                  directed=FALSE ){
  
  # transform arrays to list, as armadillo fields are required as input
  K_net <- dim(y_ijtk)[4]
  y_ijtk_list <- w_ijtk_list <- gamma_ijtk_list <- list(NULL)
  for(k in 1:K_net) {
    y_ijtk_list[[k]] <- y_ijtk[,,,k]
    w_ijtk_list[[k]] <- w_ijtk[,,,k]
    gamma_ijtk_list[[k]] <- gamma_ijtk[,,,k]
  }
  
  out_aux <- sample_coeff_tp_link_cpp( beta_tp=beta_tp,
                                       beta_t_cov_prior_inv=beta_t_cov_prior_inv,
                                       y_ijtk=y_ijtk_list, w_ijtk=w_ijtk_list, gamma_ijtk=gamma_ijtk_list,
                                       x_ijtkp_mat=x_ijtkp_mat,
                                       directed=directed )
  
  for(k in 1:K_net) {gamma_ijtk[,,,k] <- out_aux$gamma_ijtk[k,1][[1]]}
  
  return( list( beta_tp=out_aux$beta_tp,
                gamma_ijtk=gamma_ijtk ) );
  
  return(out_aux)
}

#' @importFrom MCMCpack rinvgamma
#' @keywords internal
sample_nGP_sigma <- function( alpha_t,
                              a,b,
                              delta_t ) {
  # check dimensions #
  if(dim(alpha_t)[2]!=length(delta_t)){stop("Dimension inconsistency between alpha_t and delta_t")}
  n <- dim(alpha_t)[2]
  sigma_U <- MCMCpack::rinvgamma(1, a+(n+1)/2, b+0.5*sum((alpha_t[2,-1]-(alpha_t[2,-n]+alpha_t[3,-n]*delta_t[-n]))^2/delta_t[-n]) )
  sigma_A <- MCMCpack::rinvgamma(1, a+(n+1)/2, b+0.5*sum((alpha_t[3,-1]-alpha_t[3,-n])^2/delta_t[-n]) )
  return( c(sigma_U,sigma_A) )
}


#' @keywords internal
get_nGP_sigma_net <- function( a, b,
                               time_all,
                               alpha_baseline_tk,
                               alpha_coord_ith_shared,
                               alpha_coord_ithk,
                               alpha_add_eff_it_shared=NULL,
                               alpha_add_eff_itk=NULL,
                               directed=FALSE ) {
  
  # Network dimensions
  T_net <- length(time_all)
  K_net <- dim(alpha_baseline_tk)[2]
  if( !directed ){
    V_net <- dim(alpha_coord_ith_shared)[1]
    H_dim <- dim(alpha_coord_ith_shared)[3]
  } else {
    V_net <- dim(alpha_coord_ith_shared[[1]])[1]
    H_dim <- dim(alpha_coord_ith_shared[[1]])[3]
  }
  # lenght of time intervals
  diff_time_all <- c(diff(time_all),1)
  
  # Output #
  nGP_sigma_net <- list( baseline_k = array(NA,dim=c(K_net,2)),
                         coord_i = array(NA,dim=c(V_net,2)),
                         coord_ik = array(NA,dim=c(V_net,K_net,2)),
                         add_eff_i = NULL,
                         add_eff_ik = NULL,
                         time_all=time_all,
                         directed=directed,
                         V_net=V_net, T_net=T_net, K_net=K_net, H_dim=H_dim)
  
  # Computation #
  for(k in 1:K_net) {
    # variance of baseline process #
    alpha_t = t(alpha_baseline_tk[,k,]) # dim(alpha_t)=c(3,T_net)
    nGP_sigma_net$baseline_k[k,] <- sample_nGP_sigma( alpha_t=alpha_t,
                                                      a=a,b=b,
                                                      delta_t=diff_time_all )
  }
  # variance of latent coordinates #
  for(i in 1:V_net) {
    # variance of global coordinates #
    alpha_t = alpha_coord_ith_shared[i,,,] # dim(alpha_t)=c(T_net,H_dim,3)
    alpha_t = aperm(alpha_t,perm=c(3,1,2)) # dim(alpha_t)=c(3,T_net,H_dim)
    alpha_t = matrix(c(alpha_t),3,T_net*H_dim) # dim(alpha_t)=c(3,T_net*H_dim)
    nGP_sigma_net$coord_i[i,] <- sample_nGP_sigma( alpha_t=alpha_t,
                                                   a=a,b=b,
                                                   delta_t=rep(diff_time_all,H_dim) )
    # variance of layer-specific coordinates #
    if(K_net>1){
      for(k in 1:K_net) {
        alpha_t = alpha_coord_ithk[i,,,k,] # dim(alpha_t)=c(T_net,H_dim,3)
        alpha_t = aperm(alpha_t,perm=c(3,1,2)) # dim(alpha_t)=c(3,T_net,H_dim)
        alpha_t = matrix(c(alpha_t),3,T_net*H_dim) # dim(alpha_t)=c(3,T_net*H_dim)
        nGP_sigma_net$coord_ik[i,k,] <- sample_nGP_sigma( alpha_t=alpha_t,
                                                          a=a,b=b,
                                                          delta_t=rep(diff_time_all,H_dim) )
      }
    }
  }
  
  # variance of additive effects #
  if(!is.null(alpha_add_eff_it_shared)) {
    nGP_sigma_net$add_eff_i <- array(NA,dim=c(V_net,2))
    for(i in 1:V_net) {
      alpha_t = t(alpha_add_eff_it_shared[i,,]) # dim(alpha_t)=c(3,T_net)
      nGP_sigma_net$add_eff_i[i,] <- sample_nGP_sigma( alpha_t=alpha_t,
                                                       a=a,b=b,
                                                       delta_t=diff_time_all )
    }
  }
  # variance of layer-specific additive effects #
  if(K_net>1){
    if(!is.null(alpha_add_eff_itk)) {
      nGP_sigma_net$add_eff_ik <- array(NA,dim=c(V_net,K_net,2))
      for(k in 1:K_net) {
        alpha_t = t(alpha_add_eff_itk[i,,k,]) # dim(alpha_t)=c(3,T_net)
        nGP_sigma_net$add_eff_ik[i,k,] <- sample_nGP_sigma( alpha_t=alpha_t,
                                                            a=a,b=b,
                                                            delta_t=diff_time_all )
      }
    }
  }
  
  return(nGP_sigma_net)
}
