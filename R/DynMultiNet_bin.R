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
#'    DynMultiNet_bin( net_data, pred_data,
#'                     H_dim=10,
#'                     n_iter_mcmc=100000,
#'                     out_file=NULL, log_file=NULL,
#'                     quiet_mcmc=FALSE )
#' 
#' @useDynLib DynMultiNet
#' 
#' @import BayesLogit
#' @import igraph
#' 
#' @export

DynMultiNet_bin <- function( net_data,
                             pred_data=NULL,
                             H_dim=10, k_x=0.01, k_mu=0.10, a_1=2, a_2=2.5,
                             n_iter_mcmc=100000,
                             out_file=NULL, log_file=NULL,
                             quiet_mcmc=FALSE ) {
  
  #### Start: Checking inputs ####
  colnames_net_data <- c("source","target","weight","time","layer")
  if( !all(is.element(colnames_net_data,colnames(net_data))) ) {
    stop('"net_data" must contain the following columns: ',paste(colnames_net_data,collapse=", ") )
  }
  net_data <- net_data[,colnames_net_data]
  #### End: Checking inputs ####
  
  
  #### Start: Global parameters ####
  node_all <- sort(unique(unlist(net_data[,c("source","target")])))
  V_net <- length(node_all)
  
  time_all <- sort(unique(unlist(net_data$time)))
  T_net <- length(time_all)
  
  layer_all <- sort(unique(unlist(net_data$layer)))
  K_net <- length(layer_all)
  #### End: Global parameters ####
  
  
  #### Start: MCMC initialization ####
  # Edge between actors i and j at time t in layer k
  cat("Procesing data...\n")
  y_ijtk <- array( data=NA,
                   dim=c(V_net,V_net,T_net,K_net),
                   dimnames=list(node_all,node_all,time_all,layer_all))
  for( t in 1:T_net) {
    for( k in 1:K_net) {
      #t<-1;k<-1
      y_ijtk[,,t,k][lower.tri(y_ijtk[,,t,k])] <- 0
    }
  }
  
  for( row_i in 1:nrow(net_data) ){
    # row_i <- 1
    aux_ij <- match(net_data[row_i,c("source","target")],node_all)
    i <- max(aux_ij)
    j <- min(aux_ij)
    t <- match(net_data[row_i,"time"],time_all)
    k <- match(net_data[row_i,"layer"],layer_all)
    if(net_data[row_i,"weight"]>0){
      y_ijtk[i,j,t,k] <- 1
    }
  }
  cat("done!\n")
  
  # Probability of an edge between actors i and j at time t in layer k
  pi_ijtk <- y_ijtk
  pi_ijtk[!is.na(pi_ijtk)] <- 0.5
  
  # Linear predictor for the probability of an edge between actors i and j at time t in layer k
  # all( qlogis( pi_ijtk ) == s_ijtk ) # TRUE
  s_ijtk <- y_ijtk
  s_ijtk[!is.na(s_ijtk)] <- 0
  
  # Augmented Polya-gamma data
  w_ijtk <- y_ijtk
  w_ijtk[!is.na(w_ijtk)] <- 0
  
  # Baseline parameter #
  # at time t for layer k
  mu_tk <- matrix( data=0,
                   nrow=T_net,
                   ncol=K_net )
  
  # Covariance matrix for baseline mu_t
  mu_tk_sigma <- array( data=0,
                        dim=c(T_net,T_net,K_net) )
  for(k in 1:K_net) { mu_tk_sigma[,,k]<-diag(T_net) }
  rm(k)
  
  # Covariance matrix prior for baseline mu_t
  mu_tk_sigma_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_mu){ exp(-k_mu*(x-y)^2) } )
  
  # Latent coordinates #
  # shared: hth coordinate of actor v at time t shared across the different layers
  x_iht <- array( data=0.10,
                  dim=c(V_net,H_dim,T_net) )
  
  if( K_net>1 ){
    # by layer: hth coordinate of actor v at time t specific to layer k
    x_ihtk <- array( data=0.10,
                     dim=c(V_net,H_dim,T_net,K_net) )
  }
  #### End: MCMC initialization ####
  
  #### Start: MCMC Sampling ####
  if(!quiet_mcmc){ cat("Sampling MCMC ...\n") }
  for ( iter_i in 1:n_iter_mcmc) {
    
    ### Step 1. Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
    for( i in 2:V_net) {
      for( j in 1:i) {
        for( t in 1:T_net) {
          for( k in 1:K_net) {
            # i<-1;j<-1;t<-1;k<-1
            s_ijtk[i,j,t,k] <- mu_tk[t,k] + t(x_iht[i,,t]) %*% x_iht[j,,t]
            w_ijtk[i,j,t,k] <- rpg( num=1, h=1, z=s_ijtk[i,j,t,k] )
          }
        }
      }
    }
    rm(i,j,t,k)
    
    ### Step 2. Sample mu_t_k from its conditional N-variate Gaussian posterior ###
    aux_sum_w_tk <- apply(w_ijtk,c(3,4),sum,na.rm=T)
    for( k in 1:K_net) {
      #k <- 1
      mu_tk_sigma[,,k] <- solve( diag(aux_sum_w_tk[,k]) + solve(mu_tk_sigma_prior) )
      aux_vec1 <- matrix(NA,T_net,1)
      for( t in 1:T_net) {
        aux_vec1[t,] <- sum( y_ijtk[,,t,k] - 0.5 - x_iht[,,t] %*% t(x_iht[,,t]), na.rm=TRUE )
      }
      mu_tk[,k] <- mvtnorm::rmvnorm( n=1,
                                     mean=mu_tk_sigma[,,k] %*% aux_vec1,
                                     sigma=mu_tk_sigma[,,k] )
      # plot(y=mu_tk[,k],x=time_all,type="l")
      
    }
    rm(aux_sum_w_tk,k,aux_vec1)
    
    # Step 3. For each unit, block-sample the set of time-varying latent coordinates x_iht
    browser()
    for(i in 1:V_net) {
      1
    }
    
    if(!quiet_mcmc){
      if( is.element(iter_i, floor(n_iter_mcmc*seq(0.05,1,0.05)) ) ) {
        cat(round(100*iter_i/n_iter_mcmc),"% ",sep="")
      }
    }
  }
  if(!quiet_mcmc){cat("\nMCMC finished!\n")}
  #### End: MCMC Sampling ####
  
  
  
  out <- list(NULL)
  return( out )
  
}
