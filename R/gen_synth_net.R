#' @title
#'    Generates Synthetic Data of a Network with a latent space characterization.
#'
#' @description
#'    \code{gen_synth_net} Generates a Network.
#'
#' @param node_all Vector. Id's of nodes in the network.
#'
#' @details
#'    The model assumes a latent variable approach
#'
#' @return
#'    A list with the following components:
#' \describe{
#'     \item{\code{y_ijtk}}{Array. Edges in the net.}
#'     \item{\code{p_ijtk}}{Array. Associated probabilities for edge occurence.}
#'     \item{\code{s_ijtk}}{Array. Associated linear predictor.}
#'     \item{\code{mu_tk}}{Matrix. Baseline parameter.}
#'     \item{\code{x_iht_shared}}{Array. Shared latent coordinates.}
#'     \item{\code{x_ihtk}}{Array. Layer specific latent coordinates.}
#' }
#'
#'
#' @examples
#' 
#'    gen_synth_net( node_all=seq(1,10),
#'                   time_all=seq(1,30),
#'                   layer_all=seq(1,3),
#'                   H_dim=3, R_dim=3,
#'                   k_x=0.10, k_mu=0.10, k_p=0.10,
#'                   a_1=2, a_2=2.5,
#'                   
#'                   pred_all=NULL,
#'                   z_tkp=NULL, z_ijtkp=NULL,
#'                   beta_z_layer=NULL, beta_z_edge=NULL,
#'                   pred_id_layer=NULL, pred_id_edge=NULL,
#'                   
#'                   out_file=NULL, log_file=NULL )
#' 
#' @useDynLib DynMultiNet
#' 
#' @import mvtnorm
#' 
#' @export

gen_synth_net <- function( node_all=seq(1,10),
                           time_all=seq(1,30),
                           layer_all=seq(1,3),
                           H_dim=3, R_dim=3,
                           k_x=0.10, k_mu=0.10, k_p=0.10,
                           a_1=2, a_2=2.5,
                           
                           pred_all=NULL,
                           z_tkp=NULL, z_ijtkp=NULL,
                           beta_z_layer=NULL, beta_z_edge=NULL,
                           pred_id_layer=NULL, pred_id_edge=NULL,
                           
                           out_file=NULL, log_file=NULL ) {
  
  V_net <- length(node_all)
  T_net <- length(time_all)
  K_net <- length(layer_all)
  
  
  ### Baseline parameter ###
  # at time t for layer k
  
  # Covariance matrix prior for baseline mu_tk
  mu_t_cov_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_mu){ exp(-k*(x-y)^2) } )
  
  mu_tk <- matrix( NA,
                   nrow=T_net,
                   ncol=K_net )
  for(k in 1:K_net){ # k<-1
    mu_tk[,k] <- mvtnorm::rmvnorm( n = 1,
                                   mean = matrix(0,T_net,1),
                                   sigma = mu_t_cov_prior )
  }
  
  # Shrinkage Parameters for latent Coordinates
  v_shrink_shared <- matrix(NA, nrow=H_dim, ncol=1 )
  v_shrink_shared[1,1] <- rgamma(n=1,shape=a_1,rate=1); v_shrink_shared[-1,1] <- rgamma(n=H_dim-1,shape=a_2,rate=1)
  tau_h_shared <- matrix(cumprod(v_shrink_shared), nrow=H_dim, ncol=1 )
  # 1/tau_h
  if( K_net>1 ){
    v_shrink_k <- matrix(NA, nrow=R_dim, ncol=K_net )
    v_shrink_k[1,] <- rgamma(n=K_net,shape=a_1,rate=1); v_shrink_k[-1,] <- rgamma(n=K_net*(R_dim-1),shape=a_2,rate=1)
    tau_h_k <- matrix(apply(v_shrink_k,2,cumprod), nrow=R_dim, ncol=K_net )
  } else {
    v_shrink_k <- NULL
    tau_h_k <- NULL
  }
  
  
  
  ### Latent coordinates ###
  # Covariance matrix prior for coordinates x_t
  x_t_sigma_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_x){ exp(-k*(x-y)^2) } )
  # shared: hth coordinate of actor v at time t shared across the different layers
  x_iht_shared <- array( NA,
                         dim=c(V_net,H_dim,T_net) )
  for( i in 1:V_net){
    for( h in 1:H_dim){
      x_iht_shared[i,h,] <- mvtnorm::rmvnorm( n = 1,
                                              mean = matrix(0,T_net,1),
                                              sigma = (1/tau_h_shared[h,1]) * x_t_sigma_prior )
    }
  }; rm(i,h)
  # i<-1;h<-1;plot(y=x_iht_shared[i,h,],x=time_all,type="l")
  
  if( K_net>1 ){
    # by layer: hth coordinate of actor v at time t specific to layer k
    x_ihtk <- array( NA,
                     dim=c(V_net,R_dim,T_net,K_net) )
    for( i in 1:V_net){
      for( h in 1:R_dim){
        for( k in 1:K_net){
          
          x_ihtk[i,h,,k] <- mvtnorm::rmvnorm( n = 1,
                                              mean = matrix(0,T_net,1),
                                              sigma = (1/tau_h_k[h,k]) * x_t_sigma_prior )
        }
      }
    }; rm(i,h,k)
    # i<-1;h<-1;k<-1;plot(y=x_ihtk[i,h,,k],x=time_all,type="l")
  } else {
    x_ihtk <- NULL
  }
  
  ### Edges ###
  # Initialization #
  y_ijtk <- array( data=0,
                   dim=c(V_net,V_net,T_net,K_net) )
  for( k in 1:K_net) {
    for( t in 1:T_net) { #t<-1;k<-1
      diag(y_ijtk[,,t,k]) <- NA
    }
  }; rm(k,t)
  for( k in 1:K_net) {
    for( t in 1:T_net) { #t<-1;k<-1
      y_ijtk[,,t,k][upper.tri(y_ijtk[,,t,k])] <- NA
    }
  }; rm(k,t)
  
  ### Linear Predictor ###
  s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                pred_all=pred_all, layer_all=layer_all,
                                z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge )
  
  ### Edges probabilities ###
  pi_ijtk <- plogis(s_ijtk)
  
  ### Edges ###
  if(any(is.na(pi_ijtk[!is.na(y_ijtk)]))) {stop("There's a problem generating pi_ijtk.")}
  y_ijtk[!is.na(y_ijtk)] <- rbinom(n=K_net*T_net*V_net*(V_net-1)/2,size=1,prob=pi_ijtk[!is.na(pi_ijtk)])
  # sum(y_ijtk,na.rm=T)
  
  edge_data <- data.frame( source=NA,
                           target=NA,
                           time=NA,
                           layer=NA,
                           weight=NA )
  
  for(i in 2:V_net) {
    for(j in 1:(i-1)) {
      for(t in 1:T_net) {
        for(k in 1:K_net) {
          if(y_ijtk[i,j,t,k]==1){
            edge_data_aux <- data.frame( source=node_all[i],
                                         target=node_all[j],
                                         time=time_all[t],
                                         layer=layer_all[k],
                                         weight=1 )
            edge_data <- rbind(edge_data,edge_data_aux)
          }
        }
      }
    }
  }
  edge_data <- edge_data[-1,]
  
  synth_net <- list( edge_data=edge_data,
                      y_ijtk = y_ijtk,
                      pi_ijtk = pi_ijtk,
                      s_ijtk = s_ijtk,
                      mu_tk = mu_tk,
                      x_iht_shared=x_iht_shared,
                      x_ihtk=x_ihtk,
                      tau_h_shared=tau_h_shared,
                      tau_h_k=tau_h_k )
  if(!is.null(out_file)){
    save( synth_net, file=out_file )
  }
  return( synth_net )
}
