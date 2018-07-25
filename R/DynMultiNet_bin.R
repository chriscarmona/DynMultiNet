#' @title
#'    Bayesian Learning of Dynamic Multilayer Networks with binary data
#'
#' @description
#'    \code{DynMultiNet_bin} Implements model from Durante and Dunson, 2018
#'
#' @param net_data Data frame with network information
#' @param pred_data Data frame with linked predictors information
#' @param H_dim Latent space dimension
#' @param n_iter_mcmc number of iterations for the MCMC
#' @param out_file Indicates a file (.RData) where the output should be saved
#' @param log_file Indicates a file (.txt) where the log of the process should be saved
#' @param quiet_mcmc Silent mode, no progress update
#' @param use_cpp Boolean, indicates wheter to use Cpp routines for efficient computation
#'
#' @details
#'    The model assumes a latent variable approach
#'    
#'    \code{net_data} must be a data frame with the following columns
#'        \describe{
#'        \item{\code{source}}{Start Node}
#'        \item{\code{target}}{End node}
#'        \item{\code{weight}}{Edge weight}
#'        \item{\code{time}}{time associated with the edge}
#'        \item{\code{layer}}{Layer associated with the edge}
#'        }
#'
#' @return
#'    A list with the following components:
#' \describe{
#'     \item{\code{theta_mcmc}}{Matrix with the chain of the parameters in the model.}
#' }
#'
#'
#' @examples
#' 
#'    DynMultiNet_bin( net_data,
#'                     pred_data=NULL,
#'                     H_dim=10,
#'                     k_x=0.10, k_mu=0.10, k_p=0.10,
#'                     a_1=2, a_2=2.5,
#'                     n_iter_mcmc=100000,
#'                     out_file=NULL, log_file=NULL,
#'                     quiet_mcmc=FALSE )
#' 
#' @useDynLib DynMultiNet
#' 
#' @import BayesLogit
#' @import dplyr
#' 
#' @export
#' 

DynMultiNet_bin <- function( net_data,
                             pred_data=NULL,
                             H_dim=10,
                             k_x=0.10, k_mu=0.10, k_p=0.10,
                             a_1=2, a_2=2.5,
                             n_iter_mcmc=100000,
                             out_file=NULL, log_file=NULL,
                             quiet_mcmc=FALSE,
                             use_cpp=TRUE ) {
  
  if(!is.null(pred_data)) {
    if( !all( is.element(unique(pred_data[,"layer"]),c(NA,unique(net_data$layer))) ) ) {
      stop('Layers in "pred_data" must be one of layers in "net_data"')
    }
  }
  
  #### Start: Processing data ####
  ### Network data ###
  y_ijtk <- get_y_ijtk_from_edges( net_data,
                                   directed=FALSE,
                                   weighted=FALSE,
                                   self_edges=FALSE )
  node_all <- sort(unique(unlist(net_data[,c("source","target")])))
  V_net <- length(node_all)
  time_all <- sort(unique(unlist(net_data$time)))
  T_net <- length(time_all)
  layer_all <- sort(unique(unlist(net_data$layer)))
  K_net <- length(layer_all)
  
  ### Predictors data ###
  pred_net <- get_z_pred( pred_data,
                          node_all, time_all, layer_all,
                          quiet=FALSE )
  
  pred_all <- pred_net$pred_all
  
  pred_id_layer <- pred_net$pred_id_layer
  pred_id_edge <- pred_net$pred_id_edge
  
  z_tkp<-pred_net$z_tkp
  z_ijtkp<-pred_net$z_ijtkp
  
  beta_z_layer<-pred_net$beta_z_layer
  beta_z_layer_mcmc<-NULL
  beta_z_edge<-pred_net$beta_z_edge
  beta_z_edge_mcmc<-NULL
  #### End: Processing data ####
  
  
  
  #### Start: MCMC initialization ####
  # Edge between actors i and j at time t in layer k
  
  # Augmented Polya-gamma data
  w_ijtk <- y_ijtk
  w_ijtk[!is.na(w_ijtk)] <- 0
  
  # Baseline parameter #
  # at time t for layer k
  mu_tk <- matrix( #data=0,
                   data=runif(T_net*K_net),
                   nrow=T_net,
                   ncol=K_net )
  mu_tk_mcmc <- array( NA, dim=c(T_net,K_net,n_iter_mcmc) )
  
  # Covariance matrix prior for baseline mu_tk
  mu_t_cov_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_mu){ exp(-k*(x-y)^2) } )
  mu_t_cov_prior_inv <- solve(mu_t_cov_prior)
  
  # Covariance matrix prior for parameters beta
  beta_t_cov_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_p){ exp(-k*(x-y)^2) } )
  beta_t_cov_prior_inv <- solve(beta_t_cov_prior)
  
  # Latent coordinates #
  # shared: hth coordinate of actor v at time t shared across the different layers
  x_iht_shared <- array( data=runif(V_net*H_dim*T_net,-1,1),
                         dim=c(V_net,H_dim,T_net) )
  x_iht_mat_shared <- aperm(a=x_iht_shared,perm=c(1,3,2))
  dim(x_iht_mat_shared) <- c(V_net,T_net*H_dim)
  if( !all(x_iht_mat_shared[1,1:T_net]==x_iht_shared[1,1,]) ){stop("there is a problem arranging x_iht_shared into x_iht_mat_shared")}
  
  x_iht_mat_shared_mcmc <- array(NA,c(V_net,T_net*H_dim,n_iter_mcmc))
  
  if( K_net>1 ){
    # by layer: hth coordinate of actor v at time t specific to layer k
    x_ihtk <- array( data=runif(V_net*H_dim*T_net*K_net),
                     dim=c(V_net,H_dim,T_net,K_net) )
    for( k in 1:K_net ){
      x_iht_mat_k_aux <- aperm(a=x_ihtk[,,,k],perm=c(1,3,2))
      dim(x_iht_mat_k_aux) <- c(V_net,T_net*H_dim)
      x_iht_mat_k[,,k] <- x_iht_mat_k_aux
      if( !all(x_iht_mat_k[1,1:T_net,k]==x_ihtk[1,1,,k]) ){stop("there is a problem arranging x_ihtk into x_iht_mat_k")}
    }
    rm(x_iht_mat_k_aux)
  } else {
    x_ihtk <- NULL
    x_iht_mat_k <- NULL
  }
  # Covariance matrix prior for coordinates x_t
  x_t_sigma_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_x){ exp(-k*(x-y)^2) } )
  x_t_sigma_prior_inv <- solve(x_t_sigma_prior)
  
  # Predictor coefficients
  if(!is.null(pred_data)) {
    if(!is.null(pred_id_layer)){
      beta_z_layer <- matrix(0,nrow=T_net,ncol=nrow(pred_id_layer))
      beta_z_layer_mcmc <- array(NA, dim=c(T_net,nrow(pred_id_layer),n_iter_mcmc) )
    }
    if(!is.null(pred_id_edge)){
      beta_z_edge <- matrix(0,nrow=T_net,ncol=nrow(pred_id_edge))
      beta_z_edge_mcmc <- array(NA, dim=c(T_net,nrow(pred_id_edge),n_iter_mcmc) )
    }
  }
  
  # Linear predictor for the probability of an edge between actors i and j at time t in layer k
  s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                pred_all=pred_all, layer_all=layer_all,
                                z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                beta_z_tkp=beta_z_layer, beta_z_ijtkp=beta_z_edge,
                                pred_id_tkp=pred_id_layer, pred_id_ijtkp=pred_id_edge )
  
  # Probability of an edge between actors i and j at time t in layer k
  pi_ijtk <- plogis(s_ijtk)
  #pi_ijt[,,1]
  # all( abs(qlogis( pi_ijt ) - s_ijt)<1e-6,na.rm=T ) # TRUE
  
  # Shrinkage Parameters
  v_dim <- matrix( NA, nrow=H_dim, ncol=1 )
  v_dim[1,1] <- rgamma(n=1,shape=a_1,rate=1); v_dim[-1,1] <- rgamma(n=H_dim-1,shape=a_2,rate=1)
  tau_h <- matrix(cumprod(v_dim), nrow=H_dim, ncol=1 )
  # 1/tau_h
  #### End: MCMC initialization ####
  
  
  
  #### Start: MCMC Sampling ####
  if(!quiet_mcmc){ cat("Sampling MCMC ...\n") }
  for ( iter_i in 1:n_iter_mcmc) {
    #cat(iter_i,",")
    
    
    
    ### Step 1. Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
    w_ijtk <- sample_w_ijtk_DynMultiNet_bin( w_ijtk=w_ijtk,
                                             s_ijtk=s_ijtk )
    
    
    
    ### Step 2_mu. Sample mu_tk from its conditional N-variate Gaussian posterior ###
    mu_tk <- sample_mu_tk_DynMultiNet_bin( mu_tk=mu_tk,
                                           y_ijtk=y_ijtk, w_ijtk=w_ijtk, s_ijtk=s_ijtk,
                                           mu_t_cov_prior_inv=mu_t_cov_prior_inv,
                                           use_cpp=use_cpp )
    # MCMC chain #
    mu_tk_mcmc[,,iter_i] <- mu_tk
    
    # update linear predictor
    s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                  x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                  pred_all=pred_all, layer_all=layer_all,
                                  z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                  beta_z_tkp=beta_z_layer, beta_z_ijtkp=beta_z_edge,
                                  pred_id_tkp=pred_id_layer, pred_id_ijtkp=pred_id_edge )
    
    
    
    ### Step 2_beta. Sample beta_z_layer and beta_z_edge from its conditional N-variate Gaussian posterior ###
    if(!is.null(pred_all)){
      # Layer specific
      if(!is.null(beta_z_layer)&!is.null(pred_id_layer)){
        beta_z_layer <- sample_beta_z_layer_DynMultiNet_bin( beta_z_layer,
                                                             z_tkp, pred_id_layer, pred_all, layer_all,
                                                             y_ijtk, w_ijtk, s_ijtk,
                                                             beta_t_cov_prior_inv,
                                                             use_cpp=use_cpp )
        beta_z_layer_mcmc[,,iter_i] <- beta_z_layer
        # update linear predictor
        s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                      x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                      pred_all=pred_all, layer_all=layer_all,
                                      z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                      beta_z_tkp=beta_z_layer, beta_z_ijtkp=beta_z_edge,
                                      pred_id_tkp=pred_id_layer, pred_id_ijtkp=pred_id_edge )
      }
      # Edge specific
      if(!is.null(beta_z_edge)&!is.null(pred_id_edge)){
        beta_z_edge <- sample_beta_z_edge_DynMultiNet_bin( beta_z_edge,
                                                           z_ijtkp, pred_id_edge, pred_all, layer_all,
                                                           y_ijtk, w_ijtk, s_ijtk,
                                                           beta_t_cov_prior_inv )
        beta_z_edge_mcmc[,,iter_i] <- beta_z_edge
        # update linear predictor
        s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                      x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                      pred_all=pred_all, layer_all=layer_all,
                                      z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                      beta_z_tkp=beta_z_layer, beta_z_ijtkp=beta_z_edge,
                                      pred_id_tkp=pred_id_layer, pred_id_ijtkp=pred_id_edge )
      }
    }
    
    ### Step 3. For each unit, block-sample the set of time-varying latent coordinates x_iht ###
    browser()
    ### SHARED Latent Coordinates ### 
    x_iht_mat_shared <- sample_x_iht_mat_DynMultiNet_bin( x_iht_mat=x_iht_mat_shared,
                                                          x_t_sigma_prior_inv=x_t_sigma_prior_inv, tau_h=tau_h,
                                                          y_ijtk=y_ijtk, w_ijtk=w_ijtk, s_ijtk=s_ijtk )
    # redefine x_ihtk with the new sampled values in x_iht_mat_k
    x_iht_k_aux <- x_iht_mat_shared
    dim(x_iht_k_aux) <- c(V_net,T_net,H_dim)
    x_iht_shared <- aperm(a=x_iht_k_aux,perm=c(1,3,2))
    if( !all(x_iht_mat_shared[1,1:T_net]==x_iht_shared[1,1,]) ){stop("there is a problem arranging x_iht_shared from x_iht_mat_shared")}
    
    # MCMC chain #
    x_iht_mat_shared_mcmc[,,iter_i] <- x_iht_mat_shared
    
    # update linear predictor
    s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                  x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                  pred_all=pred_all, layer_all=layer_all,
                                  z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                  beta_z_tkp=beta_z_layer, beta_z_ijtkp=beta_z_edge,
                                  pred_id_tkp=pred_id_layer, pred_id_ijtkp=pred_id_edge )
    
    
    ### LAYER SPECIFIC Latent Coordinates ###
    if( K_net > 1 ) {
      ### Step 3A. For each unit, block-sample the EDGE SPECIFIC set of time-varying latent coordinates x_ihtk ###
      for(k in 1:K_net) { # k<-1
        x_iht_mat_k[,,k] <- sample_x_iht_mat_DynMultiNet_bin( x_iht_mat=x_iht_mat_k[,,k],
                                                              x_t_sigma_prior_inv=x_t_sigma_prior_inv, tau_h=tau_h,
                                                              y_ijtk=y_ijtk[,,,k,drop=F], w_ijtk=w_ijtk[,,,k,drop=F], s_ijtk=s_ijtk[,,,k,drop=F] )
        # redefine x_ihtk with the new sampled values in x_iht_mat_k
        x_iht_k_aux <- x_iht_mat_k[,,k]
        dim(x_iht_k_aux) <- c(V_net,T_net,H_dim)
        x_ihtk[,,,k] <- aperm(a=x_iht_k_aux,perm=c(1,3,2))
        if( !all(x_iht_mat_k[1,1:T_net,k]==x_ihtk[1,1,,k]) ){stop("there is a problem arranging x_ihtk from x_iht_mat_k, k=",k)}
      }; rm(k,x_iht_k_aux)
      
      # update linear predictor
      s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                    x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                    pred_all=pred_all, layer_all=layer_all,
                                    z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                    beta_z_tkp=beta_z_layer, beta_z_ijtkp=beta_z_edge,
                                    pred_id_tkp=pred_id_layer, pred_id_ijtkp=pred_id_edge )
    }
    
    
    
    # Edge probabilities #
    pi_ijtk <- plogis(s_ijtk)
    
    
    
    ### Step 4. Sample the global shrinkage hyperparameters from conditional gamma distributions ###
    v_dim <- sample_v_dim_DynMultiNet_bin( v_dim, a_1, a_2,
                                           x_iht,
                                           x_t_sigma_prior_inv )
    tau_h <- matrix(cumprod(v_dim), nrow=H_dim, ncol=1 )
    
    
    
    # display MCMC progress #
    if( is.element(iter_i, floor(n_iter_mcmc*seq(0,1,0.05)[-1]) ) ) {
      if(!quiet_mcmc){
        cat(round(100*iter_i/n_iter_mcmc),"% ",sep="")
      }
    }
    # save MCMC progress #
    if( is.element(iter_i, floor(n_iter_mcmc*seq(0,1,0.25)[-1]) ) ) {
      if(!is.null(out_file)){
        DynMultiNet_mcmc <- list( node_all=node_all, time_all=time_all, layer_all=layer_all,
                                  mu_tk_mcmc=mu_tk_mcmc,
                                  x_iht_mat_mcmc=x_iht_mat_mcmc,
                                  beta_z_layer_mcmc=beta_z_layer_mcmc,
                                  beta_z_edge_mcmc=beta_z_edge_mcmc,
                                  pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge )
        save(DynMultiNet_mcmc,file=out_file)
      }
    }
  }
  if(!quiet_mcmc){cat("\nMCMC finished!\n")}
  #### End: MCMC Sampling ####
  
  DynMultiNet_mcmc <- list( node_all=node_all, time_all=time_all, layer_all=layer_all,
                            mu_tk_mcmc=mu_tk_mcmc,
                            x_iht_mat_mcmc=x_iht_mat_mcmc,
                            beta_z_layer_mcmc=beta_z_layer_mcmc,
                            beta_z_edge_mcmc=beta_z_edge_mcmc,
                            pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge )
  return( DynMultiNet_mcmc )
  
}
