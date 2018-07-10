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
#' @param devel Development mode
#' @param out_file Indicates a file (.RData) where the output should be saved
#' @param log_file Indicates a file (.txt) where the log of the process should be saved
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
#'                     devel=FALSE,
#'                     out_file=NULL, log_file=NULL )
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
                             devel=FALSE,
                             out_file=NULL, log_file=NULL ) {
  
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
  y_ijtk <- array( data=0,
                   dim=c(V_net,V_net,T_net,K_net) )
  
  # Probability of an edge between actors i and j at time t in layer k
  pi_ijtk <- array( data=0.5,
                    dim=c(V_net,V_net,T_net,K_net) )
  
  # Linear predictor for the probability of an edge between actors i and j at time t in layer k
  # all( qlogis( pi_ijtk ) == s_ijtk ) # TRUE
  s_ijtk <- array( data=0,
                   dim=c(V_net,V_net,T_net,K_net) )
  
  # Augmented Polya-gamma data
  w_ijtk <- array( data=0,
                   dim=c(V_net,V_net,T_net,K_net) )
  
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
  for ( iter_i in 1:n_iter_mcmc) {
    ### Step 1. Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
    for( i in 1:V_net) {
      for( j in 1:V_net) {
        for( t in 1:T_net) {
          for( k in 1:K_net) {
            # i<-1;j<-1;t<-1;k<-1
            #browser()
            #cat("i=",i,", j=",j,", t=",t,", k=",k,"\n")
            s_ijtk[i,j,t,k] <- mu_tk[t,k] + t(x_iht[i,,t]) %*% x_iht[j,,t]
            w_ijtk[i,j,t,k] <- rpg( num=1, h=1, z=s_ijtk[i,j,t,k] )
          }
        }
      }
    }
    rm(i,j,t,k)
    
    
  }
  browser()
  
  #### End: MCMC Sampling ####
  
  
  
  out <- list(NULL)
  return( out )
  
}
