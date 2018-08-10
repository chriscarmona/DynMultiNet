#' @title
#'    Generates Synthetic Data of a Network with a latent space characterization.
#'
#' @description
#'    \code{gen_synth_net} Generates a Network.
#'
#' @param node_all Vector. Id's of nodes in the network.
#' @param time_all Vector. Timestamps of network's observations.
#' @param layer_all Vector. Id's of layers in the network.
#' @param directed Boolean. Indicates if the generated network must be directed, i.e. assymetric adjacency matrix.
#' @param H_dim Integer. Latent space dimension.
#' @param R_dim Integer. Latent space dimension, for layer specific latent vectors.
#' @param k_x Positive scalar. Hyperparameter controlling for the smoothness in the dynamic of latent coordinates. Smaller=smoother.
#' @param k_mu Positive scalar. Hyperparameter controlling for the smoothness in the dynamic of the baseline process. Smaller=smoother.
#' @param k_p Positive scalar. Hyperparameter controlling for the smoothness in the dynamic of the predictor coefficients. Smaller=smoother.
#' @param a_1 Positive scalar. Hyperparameter controlling for number of effective dimensions in the latent space.
#' @param a_2 Positive scalar. Hyperparameter controlling for number of effective dimensions in the latent space.
#' @param pred_all Vector. Id's of predictors.
#' @param pred_id_layer Data Frame. Catalogue matching predictors id's and layers, for layer-specific predictors.
#' @param pred_id_edge Data Frame. Catalogue matching predictors id's and layers, for edge-specific predictors.
#' @param z_tkp Numeric array. Data for layer-specific predictor p, for time t and layer k.
#' @param z_ijtkp Numeric array. Data for edge-specific predictor p, for edge i-j at time t and layer k.
#' @param beta_z_layer Numeric array. Layer specific predictor p, for edge i-j at time t and layer k.
#' @param beta_z_edge Numeric array. Layer specific predictor p, for edge i-j at time t and layer k.
#' @param out_file String. Indicates a file (.RData) where the output will be saved.
#' 
#' @details
#'    The model assumes a latent variable approach
#'
#' @return
#'    A list with the following components:
#' \describe{
#'     \item{\code{edge_data}}{}
#'     \item{\code{y_ijtk}}{Numeric array. Network data, weight of edge between nodes i and j at time t in layer k.}
#'     \item{\code{pi_ijtk}}{Numeric array. Underlying probability of edge existence.}
#'     \item{\code{s_ijtk}}{Numeric array. Associated linear predictor in the logit model.}
#'     \item{\code{mu_tk}}{Numeric matrix. Baseline process at time t for layer k.}
#'     \item{\code{x_ith_shared}}{Numeric array. Global latent coordinates.}
#'     \item{\code{x_ithk}}{Numeric array. Latent specific latent coordinates.}
#'     \item{\code{tau_h_shared}}{Numeric matrix. Shrinkage parameter for global latent coordinates.}
#'     \item{\code{tau_h_k}}{Numeric matrix. Shrinkage parameter for layer-specific latent coordinates.}
#' }
#'
#'
#' @examples
#' 
#' synth_net <- gen_synth_net( node_all=seq(1,5),
#'                             time_all=seq(1,10),
#'                             layer_all=seq(1,3),
#'                             H_dim=3, R_dim=3,
#'                             k_x=0.10, k_mu=0.10, k_p=0.10,
#'                             a_1=2, a_2=2.5 )
#' 
#' head(synth_net$edge_data)
#' 
#' @useDynLib DynMultiNet
#' 
#' @import mvtnorm
#' @importFrom stats plogis rbinom rgamma runif
#' 
#' @export
#' 

gen_synth_net <- function( node_all,
                           time_all,
                           layer_all,
                           directed=FALSE,
                           H_dim=3, R_dim=3,
                           k_x=0.10, k_mu=0.10, k_p=0.10,
                           a_1=2, a_2=2.5,
                           
                           pred_all=NULL,
                           pred_id_layer=NULL, pred_id_edge=NULL,
                           z_tkp=NULL, z_ijtkp=NULL,
                           beta_z_layer=NULL, beta_z_edge=NULL,
                           
                           out_file=NULL ) {
  
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
  if(directed){
    v_shrink_shared <- list( sender=v_shrink_shared,
                             receiver=v_shrink_shared)
    v_shrink_shared[[2]][1,1] <- rgamma(n=1,shape=a_1,rate=1); v_shrink_shared[[2]][-1,1] <- rgamma(n=H_dim-1,shape=a_2,rate=1)
    tau_h_shared <- list( sender=tau_h_shared,
                          receiver=tau_h_shared )
    tau_h_shared[[2]] <- matrix(cumprod(v_shrink_shared[[2]]), nrow=H_dim, ncol=1 )
  }
  
  if( K_net>1 ){
    v_shrink_k <- matrix(NA, nrow=R_dim, ncol=K_net )
    v_shrink_k[1,] <- rgamma(n=K_net,shape=a_1,rate=1); v_shrink_k[-1,] <- rgamma(n=K_net*(R_dim-1),shape=a_2,rate=1)
    tau_h_k <- matrix(apply(v_shrink_k,2,cumprod), nrow=R_dim, ncol=K_net )
    if(directed){
      v_shrink_k <- list( sender=v_shrink_k,
                          receiver=v_shrink_k)
      v_shrink_k[[2]][1,] <- rgamma(n=K_net,shape=a_1,rate=1); v_shrink_k[[2]][-1,] <- rgamma(n=K_net*(R_dim-1),shape=a_2,rate=1)
      tau_h_k <- list( sender=tau_h_k,
                       receiver=tau_h_k )
      tau_h_k[[2]] <- matrix(apply(v_shrink_k[[2]],2,cumprod), nrow=R_dim, ncol=K_net )
    }
  } else {
    v_shrink_k <- NULL
    tau_h_k <- NULL
    if(directed){
      v_shrink_k <- list( sender=v_shrink_k,
                          receiver=v_shrink_k)
      tau_h_k <- list( sender=tau_h_k,
                       receiver=tau_h_k )
    }
  }
  
  
  ### Latent coordinates ###
  # Covariance matrix prior for coordinates x_t
  x_t_sigma_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_x){ exp(-k*(x-y)^2) } )
  # shared: hth coordinate of actor v at time t shared across the different layers
  x_ith_shared <- array( NA,
                         dim=c(V_net,T_net,H_dim) )
  # i<-1;h<-1;plot(y=x_ith_shared[i,,h],x=time_all,type="l")
  if(directed){
    x_ith_shared <- list( sender=x_ith_shared,
                          receiver=x_ith_shared )
    for( dir_i in 1:2 ) {
      for( i in 1:V_net){
        for( h in 1:H_dim){
          x_ith_shared[[dir_i]][i,,h] <- mvtnorm::rmvnorm( n = 1,
                                                           mean = matrix(0,T_net,1),
                                                           sigma = (1/tau_h_shared[[dir_i]][h,1]) * x_t_sigma_prior )
        }
      }
    }; rm(dir_i,i,h)
  } else {
    for( i in 1:V_net){
      for( h in 1:H_dim){
        x_ith_shared[i,,h] <- mvtnorm::rmvnorm( n = 1,
                                                mean = matrix(0,T_net,1),
                                                sigma = (1/tau_h_shared[h,1]) * x_t_sigma_prior )
      }
    }; rm(i,h)
  }
  
  if( K_net>1 ){
    # by layer: hth coordinate of actor v at time t specific to layer k
    x_ithk <- array( NA,
                     dim=c(V_net,T_net,R_dim,K_net) )
    
    # i<-1;h<-1;k<-1;plot(y=x_ithk[i,,h,k],x=time_all,type="l")
    if(directed){
      x_ithk <- list( sender=x_ithk,
                      receiver=x_ithk )
      for( dir_i in 1:2){
        for( k in 1:K_net){
          for( i in 1:V_net){
            for( h in 1:H_dim){
              x_ithk[[dir_i]][i,,h,k] <- mvtnorm::rmvnorm( n = 1,
                                                           mean = matrix(0,T_net,1),
                                                           sigma = (1/tau_h_k[[dir_i]][h,k]) * x_t_sigma_prior )
            }
          }
        }
      }; rm(dir_i,i,h)
    } else {
      for( k in 1:K_net){
        for( i in 1:V_net){
          for( h in 1:R_dim){
            x_ithk[i,,h,k] <- mvtnorm::rmvnorm( n = 1,
                                                mean = matrix(0,T_net,1),
                                                sigma = (1/tau_h_k[h,k]) * x_t_sigma_prior )
          }
        }
      }; rm(i,h,k)
    }
  } else {
    x_ithk <- NULL
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
  if(!directed){
    for( k in 1:K_net) {
      for( t in 1:T_net) { #t<-1;k<-1
        y_ijtk[,,t,k][upper.tri(y_ijtk[,,t,k])] <- NA
      }
    }; rm(k,t)
  }
  
  ### Linear Predictor ###
  s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                x_ith_shared=x_ith_shared, x_ithk=x_ithk,
                                pred_all=pred_all, layer_all=layer_all,
                                z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                directed=directed )
  
  ### Edges probabilities ###
  pi_ijtk <- plogis(s_ijtk)
  
  ### Edges ###
  probs <- pi_ijtk[!is.na(y_ijtk)]
  if(any(is.na(probs))) {stop("There's a problem generating pi_ijtk.")}
  if( directed ){
    if(length(probs)!=K_net*T_net*V_net*(V_net-1)) {stop("There's a problem generating y_ijtk")}
    y_ijtk[!is.na(y_ijtk)] <- rbinom(n=length(probs),size=1,prob=probs)
  } else {
    if(length(probs)!=K_net*T_net*V_net*(V_net-1)/2) {stop("There's a problem generating y_ijtk")}
    y_ijtk[!is.na(y_ijtk)] <- rbinom(n=length(probs),size=1,prob=probs)
  }
  # sum(y_ijtk,na.rm=T)
  
  
  # Edge list data frame #
  edge_data <- data.frame( source=NA,
                           target=NA,
                           time=NA,
                           layer=NA,
                           weight=NA )
  
  for(t in 1:T_net) {
    for(k in 1:K_net) {
      for(i in 1:V_net) {
        for(j in 1:V_net) {
          if(!is.na(y_ijtk[i,j,t,k])){
            if(y_ijtk[i,j,t,k]!=0){
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
  }
  edge_data <- edge_data[-1,]
  
  synth_net <- list( edge_data=edge_data,
                     y_ijtk = y_ijtk,
                     pi_ijtk = pi_ijtk,
                     s_ijtk = s_ijtk,
                     mu_tk = mu_tk,
                     x_ith_shared=x_ith_shared,
                     x_ithk=x_ithk,
                     tau_h_shared=tau_h_shared,
                     tau_h_k=tau_h_k )
  if(!is.null(out_file)){
    save( synth_net, file=out_file )
  }
  return( synth_net )
}
