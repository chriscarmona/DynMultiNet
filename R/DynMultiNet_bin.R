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
#' @param use_cpp whether the process is executed in C++ or not
#'
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
#' @importFrom abind abind
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
  
  #### Start: Checking inputs ####
  colnames_net_data <- c("source","target","time","layer","weight")
  if( !all(is.element(colnames_net_data,colnames(net_data))) ) {
    stop('"net_data" must contain the following columns: ',paste(colnames_net_data,collapse=", ") )
  }
  net_data <- net_data[,colnames_net_data]
  rm(colnames_net_data)
  
  if(!is.null(pred_data)){
    colnames_pred_data <- c("source","target","time","layer","id","pred")
    if( !all(is.element(colnames_pred_data,colnames(pred_data))) ) {
      stop('"pred_data" must contain the following columns: ',paste(colnames_pred_data,collapse=", ") )
    }
    pred_data <- pred_data[,colnames_pred_data]
    if( any(is.na(pred_data[,c("time","id","pred")])) ) {
      stop('"pred_data" does NOT allow NAs in columns: ',paste(c("time","id","pred"),collapse=", ") )
    }
    
    rm(colnames_pred_data)
  }
  #### End: Checking inputs ####
  
  
  #### Start: Global parameters ####
  node_all <- sort(unique(unlist(net_data[,c("source","target")])))
  V_net <- length(node_all)
  
  time_all <- sort(unique(unlist(net_data$time)))
  T_net <- length(time_all)
  
  layer_all <- sort(unique(unlist(net_data$layer)))
  K_net <- length(layer_all)
  #### End: Global parameters ####
  
  #### Start: Processing data ####
  
  ### Network data ###
  y_ijtk <- get_y_ijtk_from_edges( edges_data=net_data,
                                   quiet=FALSE )
  
  
  
  ### Predictors data ###
  pred_all <- NULL
  pred_id_global<-NULL; pred_id_layer<-NULL; pred_id_edge<-NULL
  z_tp<-NULL; z_tkp<-NULL; z_ijtkp<-NULL
  beta_z_global<-NULL; beta_z_layer<-NULL; beta_z_edge<-NULL
  if(!is.null(pred_data)) {
    #colnames_pred_data <- c("source","target","time","layer","id","pred")
    
    if( !all(is.element(unique(pred_data[,"layer"]),c(NA,unique(net_data[,"layer"])))) ) {
      stop('Layers in "pred_data" must be one of layers in "net_data"')
    }
    
    ## Predictors ##
    pred_all <- sort(unique(unlist(pred_data$id)))
    P_pred <- length(pred_all)
    
    cat("Procesing predictors data...\n")
    ## Global ##
    cond_pred <- apply( matrix(!is.na(pred_data),nrow=nrow(pred_data)), 1, identical, c(F,F,T,F,T,T) )
    if( any(cond_pred) ) {
      pred_id_global <- pred_data[cond_pred,c("id","layer")]
      pred_id_global <- pred_id_global[!duplicated(pred_id_global),]
      z_tp <- matrix(NA,nrow=T_net,ncol=P_pred)
      for(row_i in 1:sum(cond_pred)){
        # row_i <- 1
        t <- match(pred_data[row_i,"time"],time_all)
        p <- match(pred_data[row_i,"id"],pred_all)
        z_tp[t,p] <- pred_data[cond_pred,][row_i,"pred"]
      }
      rm(row_i,t,p)
    }
    ## Layer specific ##
    cond_pred <- apply( matrix(!is.na(pred_data),nrow=nrow(pred_data)), 1, identical, c(F,F,T,T,T,T) )
    if( any(cond_pred) ) {
      pred_id_layer <- pred_data[cond_pred,c("id","layer")]
      pred_id_layer <- pred_id_layer[!duplicated(pred_id_layer),]
      z_tkp <- array( NA, dim=c(T_net,K_net,P_pred) )
      for(row_i in 1:sum(cond_pred)){
        # row_i <- 1
        t <- match(pred_data[row_i,"time"],time_all)
        k <- match(pred_data[row_i,"layer"],time_all)
        p <- match(pred_data[row_i,"id"],pred_all)
        z_tp[t,p] <- pred_data[cond_pred,][row_i,"pred"]
      }
      rm(row_i,t,k,p)
    }
    ## Edge specific ##
    # The edge specific predictors MUST specify also the asscociated layer.
    # there will be one coefficient for each predictor-layer combination
    cond_pred <- apply( matrix(!is.na(pred_data),nrow=nrow(pred_data)), 1, identical, c(T,T,T,T,T,T) )
    if( any(cond_pred) ) {
      pred_id_edge <- pred_data[cond_pred,c("id","layer")]
      pred_id_edge <- pred_id_edge[!duplicated(pred_id_edge),]
      z_ijtkp <- array( NA, dim=c(V_net,V_net,T_net,K_net,P_pred) )
      for(row_i in 1:sum(cond_pred)){
        # row_i <- 1
        i <- match(pred_data[row_i,c("source")],node_all)
        j <- match(pred_data[row_i,c("target")],node_all)
        t <- match(pred_data[row_i,"time"],time_all)
        k <- match(pred_data[row_i,"layer"],time_all)
        p <- match(pred_data[row_i,"id"],pred_all)
        z_ijtkp[i,j,t,k,p] <- pred_data[cond_pred,][row_i,"pred"]
      }
      rm(row_i,i,j,t,k,p)
    }
    rm(cond_pred)
    
    cat("done!\n")
  }
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
  mu_tk_mcmc <- array( NA, dim=c(T_net,K_net,1) )
  
  # Covariance matrix prior for baseline mu_tk
  mu_t_cov_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_mu){ exp(-k*(x-y)^2) } )
  mu_t_cov_prior_inv <- solve(mu_t_cov_prior)
  
  # Covariance matrix prior for parameters beta
  beta_t_cov_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_p){ exp(-k*(x-y)^2) } )
  beta_t_cov_prior_inv <- solve(beta_t_cov_prior)
  
  # Latent coordinates #
  # shared: hth coordinate of actor v at time t shared across the different layers
  x_iht <- array( #data=0.01,
                  data=runif(V_net*H_dim*T_net,-1,1),
                  dim=c(V_net,H_dim,T_net) )
  x_iht_mat <- aperm(a=x_iht,perm=c(1,3,2))
  dim(x_iht_mat) <- c(V_net,T_net*H_dim)
  if( !all(x_iht_mat[1,1:T_net]==x_iht[1,1,]) ){stop("there is a problem arranging x_iht into x_iht_mat")}
  x_iht_mat_mcmc <- array(NA,c(V_net,T_net*H_dim,1))
  
  if( K_net>1 ){
    # by layer: hth coordinate of actor v at time t specific to layer k
    x_ihtk <- array( data=runif(V_net*H_dim*T_net*K_net),
                     dim=c(V_net,H_dim,T_net,K_net) )
  }
  # Covariance matrix prior for coordinates x_t
  x_t_sigma_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_x){ exp(-k*(x-y)^2) } )
  x_t_sigma_prior_inv <- solve(x_t_sigma_prior)
  
  # Predictor coefficients
  if(!is.null(pred_data)) {
    if(!is.null(pred_id_global)){
      beta_z_global <- matrix(0,nrow=T_net,ncol=nrow(pred_id_global))
    }
    if(!is.null(pred_id_layer)){
      beta_z_layer <- matrix(0,nrow=T_net,ncol=nrow(pred_id_layer))
    }
    if(!is.null(pred_id_edge)){
      beta_z_edge <- matrix(0,nrow=T_net,ncol=nrow(pred_id_edge))
    }
  }
  
  # Linear predictor for the probability of an edge between actors i and j at time t in layer k
  s_ijtk <- get_linpred_s_ijtk( y_ijtk, mu_tk, x_iht,
                                pred_all, layer_all,
                                z_tp, z_tkp, z_ijtkp,
                                beta_z_global, beta_z_layer, beta_z_edge,
                                pred_id_global, pred_id_layer, pred_id_edge )
  
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
    w_ijtk <- sample_w_ijtk_DynMultiNet_bin( s_ijtk )
    
    
    
    ### Step 2_mu. Sample mu_tk from its conditional N-variate Gaussian posterior ###
    mu_tk <- sample_mu_tk_DynMultiNet_bin( mu_tk,
                                           y_ijtk, w_ijtk, s_ijtk,
                                           mu_t_cov_prior_inv )
    # MCMC chain #
    mu_tk_mcmc <- abind::abind(mu_tk_mcmc,mu_tk,along=3)
    
    # update linear predictor
    s_ijtk <- get_linpred_s_ijtk( y_ijtk, mu_tk, x_iht,
                                  pred_all, layer_all,
                                  z_tp, z_tkp, z_ijtkp,
                                  beta_z_global, beta_z_layer, beta_z_edge,
                                  pred_id_global, pred_id_layer, pred_id_edge )
    
    
    
    ### Step 2_beta. Sample beta_z_global, beta_z_layer and beta_z_edge from its conditional N-variate Gaussian posterior ###
    if(!is.null(pred_all)){
      if(F&!is.null(beta_z_global)&!is.null(pred_id_global)){
        beta_z_global <- sample_beta_z_layer_DynMultiNet_bin( beta_z_global,
                                                              z_tp, pred_id_global, pred_all,
                                                              y_ijtk, w_ijtk, s_ijtk,
                                                              beta_t_cov_prior_inv )
        # update linear predictor
        s_ijtk <- get_linpred_s_ijtk( y_ijtk, mu_tk, x_iht,
                                      pred_all, layer_all,
                                      z_tp, z_tkp, z_ijtkp,
                                      beta_z_global, beta_z_layer, beta_z_edge,
                                      pred_id_global, pred_id_layer, pred_id_edge )
      }
      if(F&!is.null(beta_z_layer)&!is.null(pred_id_layer)){
        beta_z_layer <- sample_beta_z_layer_DynMultiNet_bin( beta_z_layer,
                                                             z_tkp, pred_id_layer, pred_all,
                                                             y_ijtk, w_ijtk, s_ijtk,
                                                             beta_t_cov_prior_inv )
        # update linear predictor
        s_ijtk <- get_linpred_s_ijtk( y_ijtk, mu_tk, x_iht,
                                      pred_all, layer_all,
                                      z_tp, z_tkp, z_ijtkp,
                                      beta_z_global, beta_z_layer, beta_z_edge,
                                      pred_id_global, pred_id_layer, pred_id_edge )
      }
      if(!is.null(beta_z_edge)&!is.null(pred_id_edge)){
        beta_z_edge <- sample_beta_z_edge_DynMultiNet_bin( beta_z_edge,
                                                           z_ijtkp, pred_id_edge, pred_all,
                                                           y_ijtk, w_ijtk, s_ijtk,
                                                           beta_t_cov_prior_inv )
        # update linear predictor
        s_ijtk <- get_linpred_s_ijtk( y_ijtk, mu_tk, x_iht,
                                      pred_all, layer_all,
                                      z_tp, z_tkp, z_ijtkp,
                                      beta_z_global, beta_z_layer, beta_z_edge,
                                      pred_id_global, pred_id_layer, pred_id_edge )
      }
    }
    
    ### Step 3. For each unit, block-sample the set of time-varying latent coordinates x_iht ###
    x_iht_mat <- sample_x_iht_mat_DynMultiNet_bin( x_iht_mat,
                                                   x_t_sigma_prior_inv, tau_h,
                                                   y_ijtk, w_ijtk, s_ijtk, mu_tk,
                                                   use_cpp=use_cpp )
    # MCMC chain #
    x_iht_mat_mcmc <- abind::abind(x_iht_mat_mcmc,x_iht_mat,along=3)
    
    # redefine x_iht with the new sampled values in x_iht_mat
    x_iht <- x_iht_mat
    dim(x_iht) <- c(V_net,T_net,H_dim)
    x_iht <- aperm(a=x_iht,perm=c(1,3,2))
    if( !all(x_iht_mat[1,1:T_net]==x_iht[1,1,]) ){stop("there is a problem arranging x_iht into x_iht_mat")}
    
    # update linear predictor
    s_ijtk <- get_linpred_s_ijtk( y_ijtk, mu_tk, x_iht,
                                  pred_all, layer_all,
                                  z_tp, z_tkp, z_ijtkp,
                                  beta_z_global, beta_z_layer, beta_z_edge,
                                  pred_id_global, pred_id_layer, pred_id_edge )
    
    # Edge probabilities
    pi_ijtk <- plogis(s_ijtk)
    
    
    
    ### Step 4. Sample the global shrinkage hyperparameters from conditional gamma distributions ###
    v_dim <- sample_v_dim_DynMultiNet_bin( v_dim, a_1, a_2,
                                           x_iht,
                                           x_t_sigma_prior_inv )
    tau_h <- matrix(cumprod(v_dim), nrow=H_dim, ncol=1 )
    
    
    
    # display MCMC progress #
    if( is.element(iter_i, floor(n_iter_mcmc*seq(0.05,1,0.05)) ) ) {
      if(!is.null(out_file)){
        DynMultiNet_mcmc <- list( mu_tk_mcmc=mu_tk_mcmc[,,-1],
                                  x_iht_mat_mcmc=x_iht_mat_mcmc[,,-1],
                                  V_net=V_net, T_net=T_net, H_dim=H_dim )
        save(DynMultiNet_mcmc,file=out_file)
      }
      if(!quiet_mcmc){
        cat(round(100*iter_i/n_iter_mcmc),"% ",sep="")
      }
    }
    
  }
  if(!quiet_mcmc){cat("\nMCMC finished!\n")}
  #### End: MCMC Sampling ####
  
  
  
  DynMultiNet_mcmc <- list( mu_tk_mcmc=mu_tk_mcmc,
                            x_iht_mat_mcmc=x_iht_mat_mcmc,
                            V_net=V_net, T_net=T_net, H_dim=H_dim )
  return( DynMultiNet_mcmc )
  
}
