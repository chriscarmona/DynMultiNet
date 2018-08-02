#' @title
#'    Bayesian Learning of Dynamic Multilayer Networks with binary directed edges
#'
#' @description
#'    \code{DynMultiNet_bin} Implements model from Durante and Dunson, 2018
#'
#' @param net_data Data frame.Network information.
#' @param pred_data Data frame. Linked predictors information.
#' @param directed Boolean. Indicates if the provided network is directed, i.e. the adjacency matrix is assymetrical.
#' @param H_dim Integer. Latent space dimension.
#' @param R_dim Integer. Latent space dimension, for layer specific latent vectors.
#' @param n_iter_mcmc Integer. Number of iterations for the MCMC.
#' @param n_burn Integer. Number of iterations discarded as part of the MCMC warming up period at the beginning of the chain.
#' @param n_thin Integer. Number of iterations discarded for thining the chain (reducing the autocorrelation). We keep 1 of every n_thin iterations.
#' @param out_file String. Indicates a file (.RData) where the output should be saved.
#' @param log_file String. Indicates a file (.txt) where the log of the process should be saved.
#' @param quiet_mcmc Boolean. Indicates if silent mode is preferes, if \code{FALSE} progress update is displayed.
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
#'                     directed=TRUE,
#'                     H_dim=10, R_dim=10,
#'                     k_x=0.10, k_mu=0.10, k_p=0.10,
#'                     a_1=2, a_2=2.5,
#'                     n_iter_mcmc=10000, n_burn=n_iter_mcmc/2, n_thin=3,
#'                     out_file=NULL, log_file=NULL,
#'                     quiet_mcmc=FALSE )
#' 
#' @useDynLib DynMultiNet
#' 
#' @import foreach
#' @import BayesLogit
#' @import dplyr
#' 
#' @export
#' 

DynMultiNet_bin_dir <- function( net_data,
                                 pred_data=NULL,
                                 directed=TRUE,
                                 H_dim=10, R_dim=10,
                                 k_x=0.10, k_mu=0.10, k_p=0.10,
                                 a_1=2, a_2=2.5,
                                 n_iter_mcmc=10000, n_burn=n_iter_mcmc/2, n_thin=3,
                                 out_file=NULL, log_file=NULL,
                                 quiet_mcmc=FALSE ) {
  
  mcmc_clock <- Sys.time()
  ### iterations that will be reported ###
  # after burn-in period and thinning
  iter_out_mcmc <- seq(from=n_burn+1,to=n_iter_mcmc,by=n_thin)
  n_iter_mcmc_out <- length(iter_out_mcmc)
  
  if(!is.null(pred_data)) {
    if( !all( is.element(unique(pred_data[,"layer"]),c(NA,unique(net_data$layer))) ) ) {
      stop('Layers in "pred_data" must be one of layers in "net_data"')
    }
  }
  
  #### Start: Processing data ####
  ### Network data ###
  y_ijtk <- get_y_ijtk_from_edges( net_data,
                                   directed=directed,
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
  mu_tk_mcmc <- array( NA, dim=c(T_net,K_net,n_iter_mcmc_out) )
  
  # Covariance matrix prior for baseline mu_tk
  mu_t_cov_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_mu){ exp(-k*(x-y)^2) } )
  mu_t_cov_prior_inv <- solve(mu_t_cov_prior)
  
  # Covariance matrix prior for parameters beta
  beta_t_cov_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_p){ exp(-k*(x-y)^2) } )
  beta_t_cov_prior_inv <- solve(beta_t_cov_prior)
  
  # Latent coordinates #
  # shared: hth coordinate of actor v at time t shared across the different layers
  x_iht_shared <- array( #data=0,
    data=runif(V_net*H_dim*T_net,-1,1),
    dim=c(V_net,H_dim,T_net) )
  x_iht_mat_shared <- aperm(a=x_iht_shared,perm=c(1,3,2))
  dim(x_iht_mat_shared) <- c(V_net,T_net*H_dim)
  if( !all(x_iht_mat_shared[1,1:T_net]==x_iht_shared[1,1,]) ){stop("there is a problem arranging x_iht_shared into x_iht_mat_shared")}
  x_iht_shared_mcmc <- array(NA,c(V_net,H_dim,T_net,n_iter_mcmc_out))
  
  if( K_net>1 ){
    # by layer: hth coordinate of actor v at time t specific to layer k
    x_ihtk <- array( #data=0,
      data=runif(V_net*R_dim*T_net*K_net),
      dim=c(V_net,R_dim,T_net,K_net) )
    x_iht_mat_k <- array(NA,dim=c(V_net,R_dim*T_net,K_net))
    for( k in 1:K_net ){
      x_iht_mat_k_aux <- aperm(a=x_ihtk[,,,k],perm=c(1,3,2))
      dim(x_iht_mat_k_aux) <- c(V_net,T_net*R_dim)
      x_iht_mat_k[,,k] <- x_iht_mat_k_aux
      if( !all(x_iht_mat_k[1,1:T_net,k]==x_ihtk[1,1,,k]) ){stop("there is a problem arranging x_ihtk into x_iht_mat_k")}
    }
    rm(x_iht_mat_k_aux)
    x_ihtk_mcmc <- array(NA,c(V_net,R_dim,T_net,K_net,n_iter_mcmc_out))
  } else {
    x_ihtk <- NULL
    x_iht_mat_k <- NULL
    x_ihtk_mcmc <- NULL
  }
  if( directed ){
    x_iht_shared <- list( send=x_iht_shared,
                          receive=x_iht_shared )
    x_iht_mat_shared <- list( send=x_iht_mat_shared,
                              receive=x_iht_mat_shared )
    x_iht_shared_mcmc <- list( send=x_iht_shared_mcmc,
                               receive=x_iht_shared_mcmc )
    if( K_net>1 ){
      x_ihtk <- list( send=x_ihtk,
                      receive=x_ihtk )
      x_iht_mat_k <- list( send=x_iht_mat_k,
                           receive=x_iht_mat_k )
      x_ihtk_mcmc <- list( send=x_ihtk_mcmc,
                           receive=x_ihtk_mcmc )
    }
  }
  
  # Covariance matrix prior for coordinates x_t
  x_t_sigma_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_x){ exp(-k*(x-y)^2) } )
  x_t_sigma_prior_inv <- solve(x_t_sigma_prior)
  
  # Predictor coefficients
  if(!is.null(pred_data)) {
    if(!is.null(pred_id_layer)){
      beta_z_layer <- matrix(0,nrow=T_net,ncol=nrow(pred_id_layer))
      beta_z_layer_mcmc <- array(NA, dim=c(T_net,nrow(pred_id_layer),n_iter_mcmc_out) )
    }
    if(!is.null(pred_id_edge)){
      beta_z_edge <- matrix(0,nrow=T_net,ncol=nrow(pred_id_edge))
      beta_z_edge_mcmc <- array(NA, dim=c(T_net,nrow(pred_id_edge),n_iter_mcmc_out) )
    }
  }
  
  # Linear predictor for the probability of an edge between actors i and j at time t in layer k
  s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                pred_all=pred_all, layer_all=layer_all,
                                z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                directed=directed )
  
  # Probability of an edge between actors i and j at time t in layer k
  pi_ijtk_mcmc <- array(NA, dim=c(V_net,V_net,T_net,K_net,n_iter_mcmc_out))
  
  #pi_ijt[,,1]
  # all( abs(qlogis( pi_ijt ) - s_ijt)<1e-6,na.rm=T ) # TRUE
  
  # Shrinkage Parameters
  v_shrink_shared <- matrix(NA, nrow=H_dim, ncol=1 )
  v_shrink_shared[1,1] <- rgamma(n=1,shape=a_1,rate=1); v_shrink_shared[-1,1] <- rgamma(n=H_dim-1,shape=a_2,rate=1)
  tau_h_shared <- matrix(cumprod(v_shrink_shared), nrow=H_dim, ncol=1 )
  tau_h_shared_mcmc <- matrix(NA, nrow=H_dim, ncol=n_iter_mcmc_out )
  # 1/tau_h
  if( K_net>1 ){
    v_shrink_k <- matrix(NA, nrow=R_dim, ncol=K_net )
    v_shrink_k[1,] <- rgamma(n=K_net,shape=a_1,rate=1); v_shrink_k[-1,] <- rgamma(n=K_net*(R_dim-1),shape=a_2,rate=1)
    tau_h_k <- matrix(apply(v_shrink_k,2,cumprod), nrow=R_dim, ncol=K_net )
    tau_h_k_mcmc <- array( NA, dim=c(R_dim,K_net,n_iter_mcmc_out) )
  } else {
    v_shrink_k <- NULL
    tau_h_k <- NULL
    tau_h_k_mcmc <- NULL
  }
  if( directed ){
    v_shrink_shared <- list( send=v_shrink_shared,
                             receive=v_shrink_shared )
    tau_h_shared <- list( send=tau_h_shared,
                          receive=tau_h_shared )
    tau_h_shared_mcmc <- list( send=tau_h_shared_mcmc,
                               receive=tau_h_shared_mcmc )
    if( K_net>1 ){
      v_shrink_k <- list( send=v_shrink_k,
                          receive=v_shrink_k )
      tau_h_k <- list( send=tau_h_k,
                       receive=tau_h_k )
      tau_h_k_mcmc <- list( send=tau_h_k_mcmc,
                            receive=tau_h_k_mcmc )
    }
  }
  #### End: MCMC initialization ####
  
  if( !is.null(log_file) ) {
    model_des <- "Dynamic "
    if(K_net==1) {
      model_des <- paste(model_des,"single-layer network, ",collapse="")
    } else if(K_net>1) {
      model_des <- paste(model_des,"multi-layer network, ",collapse="")
    }
    if(directed) {
      model_des <- paste(model_des,"directed unweighted edges",collapse="")
    } else {
      model_des <- paste(model_des,"undirected unweighted edges",collapse="")
    }
    cat("**** DynMultiNet_bin *****\n\n",
        "----- Network topology -----\n",
        "Nodes = ",V_net,"\n",
        "Layers = ",K_net,"\n",
        "Times steps = ",T_net,"\n",
        "directed = ",directed,"\n",
        "weighted = ",FALSE,"\n",
        "----- Inferential parameters -----\n",
        "H_dim = ",H_dim,"\n",
        "R_dim = ",R_dim,"\n",
        "k_x = ",k_x,"\n",
        "k_mu = ",k_mu,"\n",
        "k_p = ",k_p,"\n",
        "a_1 = ",a_1,"\n",
        "a_2 = ",a_2,"\n",
        "----- MCMC parameters -----\n",
        "n_iter_mcmc = ",n_iter_mcmc,"\n",
        "n_burn = ",n_burn,"\n",
        "n_thin = ",n_thin,"\n",
        "----- Storage and processing -----\n",
        "out_file = ",out_file,"\n",
        "log_file = ",log_file,"\n",
        "---------------------------\n\n",
        "Process starting time:\n",as.character(mcmc_clock),"\n\n",
        "---------------------------\n\n",
        "MCMC Starting time:\n",as.character(Sys.time()),"\n\n",
        "---------------------------\n\n\n",
        "iter_i , mcmc_acum_minutes , Sys.time \n",
        file=log_file )
  }
  mcmc_clock <- Sys.time()
  
  
  
  #### Start: MCMC Sampling ####
  if(!quiet_mcmc){ cat("Sampling MCMC ...\n") }
  for ( iter_i in 1:n_iter_mcmc) { # iter_i <- 1
    #cat(iter_i,",")
    
    
    
    ### Step 1. Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
    w_ijtk <- sample_w_ijtk_DynMultiNet_bin( w_ijtk=w_ijtk,
                                             s_ijtk=s_ijtk,
                                             directed=directed )
    
    
    
    ### Step 2_mu. Sample mu_tk from its conditional N-variate Gaussian posterior ###
    mu_tk <- sample_mu_tk_DynMultiNet_bin( mu_tk=mu_tk,
                                           y_ijtk=y_ijtk, w_ijtk=w_ijtk, s_ijtk=s_ijtk,
                                           mu_t_cov_prior_inv=mu_t_cov_prior_inv,
                                           directed=directed,
                                           use_cpp=TRUE )
    # MCMC chain #
    if(is.element(iter_i,iter_out_mcmc)){
      mu_tk_mcmc[,,match(iter_i,iter_out_mcmc)] <- mu_tk
    }
    
    # update linear predictor
    s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                  x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                  pred_all=pred_all, layer_all=layer_all,
                                  z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                  beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                  pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                  directed=directed )
    
    
    
    ### Step 2_beta. Sample beta_z_layer and beta_z_edge from its conditional N-variate Gaussian posterior ###
    if(!is.null(pred_all)){
      # Layer specific
      if(!is.null(beta_z_layer)&!is.null(pred_id_layer)){
        beta_z_layer <- sample_beta_z_layer_DynMultiNet_bin( beta_z_layer,
                                                             z_tkp, pred_id_layer, pred_all, layer_all,
                                                             y_ijtk, w_ijtk, s_ijtk,
                                                             beta_t_cov_prior_inv,
                                                             use_cpp=TRUE )
        if(is.element(iter_i,iter_out_mcmc)){
          beta_z_layer_mcmc[,,match(iter_i,iter_out_mcmc)] <- beta_z_layer
        }
        
        # update linear predictor
        s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                      x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                      pred_all=pred_all, layer_all=layer_all,
                                      z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                      beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                      pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                      directed=directed )
      }
      # Edge specific
      if(!is.null(beta_z_edge)&!is.null(pred_id_edge)){
        beta_z_edge <- sample_beta_z_edge_DynMultiNet_bin( beta_z_edge,
                                                           z_ijtkp, pred_id_edge, pred_all, layer_all,
                                                           y_ijtk, w_ijtk, s_ijtk,
                                                           beta_t_cov_prior_inv )
        if(is.element(iter_i,iter_out_mcmc)){
          beta_z_edge_mcmc[,,match(iter_i,iter_out_mcmc)] <- beta_z_edge
        }
        # update linear predictor
        s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                      x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                      pred_all=pred_all, layer_all=layer_all,
                                      z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                      beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                      pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                      directed=directed )
      }
    }
    
    ### Step 3. For each unit, block-sample the set of time-varying latent coordinates x_iht ###
    ### SHARED Latent Coordinates ###
    x_iht_mat_shared <- sample_x_iht_mat_DynMultiNet_bin( x_iht_mat=x_iht_mat_shared,
                                                          x_t_sigma_prior_inv=x_t_sigma_prior_inv,
                                                          tau_h=tau_h_shared,
                                                          y_ijtk=y_ijtk, w_ijtk=w_ijtk, s_ijtk=s_ijtk,
                                                          directed=directed )
    
    if( directed ) {
      for( dir_i in 1:2) { # dir_i <- 1
        # redefine x_iht_shared with the new sampled values in x_iht_mat_shared
        x_iht_aux <- x_iht_mat_shared[[dir_i]]
        dim(x_iht_aux) <- c(V_net,T_net,H_dim)
        x_iht_shared[[dir_i]] <- aperm(a=x_iht_aux,perm=c(1,3,2))
        rm(x_iht_aux)
        if( !all(x_iht_mat_shared[[dir_i]][1,1:T_net]==x_iht_shared[[dir_i]][1,1,]) ){ stop("there is a problem arranging x_iht_shared from x_iht_mat_shared") }
        
        # MCMC chain #
        if(is.element(iter_i,iter_out_mcmc)){
          x_iht_shared_mcmc[[dir_i]][,,,match(iter_i,iter_out_mcmc)] <- x_iht_shared[[dir_i]]
        }
      }
    } else {
      # redefine x_iht_shared with the new sampled values in x_iht_mat_shared
      x_iht_aux <- x_iht_mat_shared
      dim(x_iht_aux) <- c(V_net,T_net,H_dim)
      x_iht_shared <- aperm(a=x_iht_aux,perm=c(1,3,2))
      rm(x_iht_aux)
      if( !all(x_iht_mat_shared[1,1:T_net]==x_iht_shared[1,1,]) ){ stop("there is a problem arranging x_iht_shared from x_iht_mat_shared") }
      
      # MCMC chain #
      if(is.element(iter_i,iter_out_mcmc)){
        x_iht_shared_mcmc[,,,match(iter_i,iter_out_mcmc)] <- x_iht_shared
      }
    }
    
    # update linear predictor
    s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                  x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                  pred_all=pred_all, layer_all=layer_all,
                                  z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                  beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                  pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                  directed=directed )
    
    
    ### LAYER SPECIFIC Latent Coordinates ###
    if( K_net>1 ) {
      ### Step 3A. For each unit, block-sample the EDGE SPECIFIC set of time-varying latent coordinates x_ihtk ###
      for( k in 1:K_net ) {
        if(directed){
          x_iht_mat_k_aux <- sample_x_iht_mat_DynMultiNet_bin( x_iht_mat=list(x_iht_mat_k[[1]][,,k],x_iht_mat_k[[2]][,,k]),
                                                               x_t_sigma_prior_inv=x_t_sigma_prior_inv,
                                                               tau_h=list(tau_h_k[[1]][,k,drop=F],tau_h_k[[2]][,k,drop=F]),
                                                               y_ijtk=y_ijtk[,,,k,drop=F], w_ijtk=w_ijtk[,,,k,drop=F], s_ijtk=s_ijtk[,,,k,drop=F],
                                                               directed=directed )
          x_iht_mat_k[[1]][,,k] <- x_iht_mat_k_aux[[1]]
          x_iht_mat_k[[2]][,,k] <- x_iht_mat_k_aux[[2]]
          for( dir_i in 1:2) { # dir_i <- 1
            # redefine x_ihtk with the new sampled values in x_iht_mat_k
            x_iht_aux <- x_iht_mat_k[[dir_i]][,,k]
            dim(x_iht_aux) <- c(V_net,T_net,R_dim)
            x_ihtk[[dir_i]][,,,k] <- aperm(a=x_iht_aux,perm=c(1,3,2))
            rm(x_iht_aux)
            if( !all(x_iht_mat_k[[dir_i]][1,1:T_net,k]==x_ihtk[[dir_i]][1,1,,k]) ){stop("there is a problem arranging x_ihtk from x_iht_mat_k, k=",k)}
          }
        } else {
          x_iht_mat_k[,,k] <- sample_x_iht_mat_DynMultiNet_bin( x_iht_mat=x_iht_mat_k[,,k],
                                                                x_t_sigma_prior_inv=x_t_sigma_prior_inv,
                                                                tau_h=tau_h_k[,k,drop=F],
                                                                y_ijtk=y_ijtk[,,,k,drop=F], w_ijtk=w_ijtk[,,,k,drop=F], s_ijtk=s_ijtk[,,,k,drop=F],
                                                                directed=directed )
          
          # redefine x_ihtk with the new sampled values in x_iht_mat_k
          x_iht_aux <- x_iht_mat_k[,,k]
          dim(x_iht_aux) <- c(V_net,T_net,R_dim)
          x_ihtk[,,,k] <- aperm(a=x_iht_aux,perm=c(1,3,2))
          rm(x_iht_aux)
          if( !all(x_iht_mat_k[1,1:T_net,k]==x_ihtk[1,1,,k]) ){stop("there is a problem arranging x_ihtk from x_iht_mat_k, k=",k)}
        }
      }; rm(k)
      
      # MCMC chain for x_ihtk #
      if(is.element(iter_i,iter_out_mcmc)){
        if(directed){
          x_ihtk_mcmc[[1]][,,,,match(iter_i,iter_out_mcmc)] <- x_ihtk[[1]]
          x_ihtk_mcmc[[2]][,,,,match(iter_i,iter_out_mcmc)] <- x_ihtk[[2]]
        } else {
          x_ihtk_mcmc[,,,,match(iter_i,iter_out_mcmc)] <- x_ihtk
        }
      }
      
      # update linear predictor
      s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                    x_iht_shared=x_iht_shared, x_ihtk=x_ihtk,
                                    pred_all=pred_all, layer_all=layer_all,
                                    z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                    beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                    pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                    directed=directed )
    }
    
    
    
    # Edge probabilities #
    if(is.element(iter_i,iter_out_mcmc)){
      pi_ijtk_mcmc[,,,,match(iter_i,iter_out_mcmc)] <- plogis(s_ijtk)
    }
    
    
    
    ### Step 4. Sample the global shrinkage hyperparameters from conditional gamma distributions ###
    if( directed ) {
      for( dir_i in 1:2 ) {
        v_shrink_shared[[dir_i]] <- sample_v_shrink_DynMultiNet_bin( v_shrink_shared[[dir_i]],
                                                                     a_1, a_2,
                                                                     x_iht_shared[[dir_i]],
                                                                     x_t_sigma_prior_inv )
        tau_h_shared[[dir_i]] <- matrix(cumprod(v_shrink_shared[[dir_i]]), nrow=H_dim, ncol=1 )
        if(is.element(iter_i,iter_out_mcmc)){
          tau_h_shared_mcmc[[dir_i]][,match(iter_i,iter_out_mcmc)] <- tau_h_shared[[dir_i]]
        }
      }
    } else {
      v_shrink_shared <- sample_v_shrink_DynMultiNet_bin( v_shrink_shared,
                                                          a_1, a_2,
                                                          x_iht_shared,
                                                          x_t_sigma_prior_inv )
      tau_h_shared <- matrix(cumprod(v_shrink_shared), nrow=H_dim, ncol=1 )
      if(is.element(iter_i,iter_out_mcmc)){
        tau_h_shared_mcmc[,match(iter_i,iter_out_mcmc)] <- tau_h_shared
      }
    }
    
    if(K_net>1){
      if( directed ){
        for( dir_i in 1:2 ) {
          for(k in 1:K_net) {
            v_shrink_k[[dir_i]][,k] <- sample_v_shrink_DynMultiNet_bin( v_shrink_k[[dir_i]][,k,drop=F],
                                                                        a_1, a_2,
                                                                        x_ihtk[[dir_i]][,,,k],
                                                                        x_t_sigma_prior_inv )
          }
          tau_h_k[[dir_i]] <- matrix(apply(v_shrink_k[[dir_i]],2,cumprod), nrow=R_dim, ncol=K_net )
          if(is.element(iter_i,iter_out_mcmc)){
            tau_h_k_mcmc[[dir_i]][,,match(iter_i,iter_out_mcmc)] <- tau_h_k[[dir_i]]
          }
        }
      } else {
        for(k in 1:K_net) {
          v_shrink_k[,k] <- sample_v_shrink_DynMultiNet_bin( v_shrink_k[,k,drop=F],
                                                             a_1, a_2,
                                                             x_ihtk[,,,k],
                                                             x_t_sigma_prior_inv )
        }
        tau_h_k <- matrix(apply(v_shrink_k,2,cumprod), nrow=R_dim, ncol=K_net )
        if(is.element(iter_i,iter_out_mcmc)){
          tau_h_k_mcmc[,,match(iter_i,iter_out_mcmc)] <- tau_h_k
        }
      }
    }
    
    # display MCMC progress #
    if( is.element(iter_i, floor(n_iter_mcmc*seq(0,1,0.05)[-1]) ) ) {
      if(!quiet_mcmc){
        cat(round(100*iter_i/n_iter_mcmc),"% ",sep="")
      }
      if( !is.null(log_file) ) {
        cat(iter_i," , ", as.numeric(difftime(Sys.time(),mcmc_clock,units="mins"))," , ", as.character(Sys.time()),"\n",
            file=log_file,append=TRUE )
      }
    }
    
    # save MCMC progress #
    if( is.element(iter_i, floor(n_iter_mcmc*seq(0,1,0.25)[-1]) ) & iter_i>min(iter_out_mcmc,na.rm=T) ) {
      if(!is.null(out_file)){
        DynMultiNet_mcmc <- list( y_ijtk=y_ijtk,
                                  pi_ijtk_mcmc=pi_ijtk_mcmc,
                                  node_all=node_all, time_all=time_all, layer_all=layer_all,
                                  mu_tk_mcmc=mu_tk_mcmc,
                                  x_iht_shared_mcmc=x_iht_shared_mcmc,
                                  x_ihtk_mcmc=x_ihtk_mcmc,
                                  tau_h_shared_mcmc=tau_h_shared_mcmc,
                                  tau_h_k_mcmc=tau_h_k_mcmc,
                                  pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                  beta_z_layer_mcmc=beta_z_layer_mcmc,
                                  beta_z_edge_mcmc=beta_z_edge_mcmc )
        save(DynMultiNet_mcmc,file=out_file)
      }
    }
  }
  if(!quiet_mcmc){cat("\nMCMC finished!\n")}
  #### End: MCMC Sampling ####
  if( !is.null(log_file) ) {
    cat("\n\n---------------------------\n\n",
        "Finishing time:\n",as.character(Sys.time()),"\n\n",
        file=log_file, append=TRUE)
  }
  
  DynMultiNet_mcmc <- list( y_ijtk=y_ijtk,
                            pi_ijtk_mcmc=pi_ijtk_mcmc,
                            node_all=node_all, time_all=time_all, layer_all=layer_all,
                            mu_tk_mcmc=mu_tk_mcmc,
                            x_iht_shared_mcmc=x_iht_shared_mcmc,
                            x_ihtk_mcmc=x_ihtk_mcmc,
                            tau_h_shared_mcmc=tau_h_shared_mcmc,
                            tau_h_k_mcmc=tau_h_k_mcmc,
                            pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                            beta_z_layer_mcmc=beta_z_layer_mcmc,
                            beta_z_edge_mcmc=beta_z_edge_mcmc )
  return( DynMultiNet_mcmc )
  
}
