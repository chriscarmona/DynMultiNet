#' @title
#'    Generates Synthetic Data of a Network with a latent space characterization.
#'
#' @description
#'    \code{gen_synth_net} Generates a Network.
#'
#' @param nodes_net Vector. Id's of nodes in the network.
#' @param times_net Vector. Timestamps of network's observations.
#' @param layers_net Vector. Id's of layers in the network.
#' @param directed Boolean. Indicates if the provided network is directed, i.e. the adjacency matrix is assymetrical.
#' @param weighted Boolean. Indicates if the provided network is weighted, i.e. edges with values other that 0 and 1.
#' @param H_dim Integer. Latent space dimension.
#' @param R_dim Integer. Latent space dimension, for layer specific latent vectors.
#' @param add_eff_weight Boolean. Indicates if dynamic additive effects by node should be considered for edge weights.
#' @param add_eff_link Boolean. Indicates if dynamic additive effects by node should be considered for links.
#' @param nu_sigma Positive scalar. Hyperparameter controlling the prior of the weight variance.
#' @param tau_sigma Positive scalar. Hyperparameter controlling the prior of the weight variance.
#' @param nu_delta Positive scalar. Hyperparameter controlling the prior for smoothness in the dynamic of latent elements.
#' @param tau_delta Positive scalar. Hyperparameter controlling the prior for smoothness in the dynamic of latent elements.
#' @param x_ijtkp Numeric array. External covariates (edge specific) to inform network.
#' @param rds_file String. Indicates a file (.rds) where the output will be saved.
#' 
#' @details
#'    The model assumes a latent variable approach
#'
#' @return
#'    A list with the following components:
#' \describe{
#'     \item{\code{y_ijtk}}{Numeric array. Network data, weight of edge between nodes i and j at time t in layer k.}
#'     \item{\code{pi_ijtk}}{Numeric array. Underlying probability of link existence.}
#'     \item{\code{gamma_ijtk}}{Numeric array. Associated linear predictor in the logit model.}
#'     \item{\code{eta_tk}}{Numeric matrix. Baseline process at time t for layer k.}
#'     \item{\code{ab_ith_shared}}{Numeric array. Global latent coordinates.}
#'     \item{\code{ab_ithk}}{Numeric array. Latent specific latent coordinates.}
#'     \item{\code{tau_h_shared}}{Numeric matrix. Shrinkage parameter for global latent coordinates.}
#'     \item{\code{tau_h_k}}{Numeric matrix. Shrinkage parameter for layer-specific latent coordinates.}
#' }
#'
#'
#' @examples
#' 
#' synth_net <- gen_synth_net( nodes_net=seq(1,5),
#'                             times_net=seq(1,10),
#'                             layers_net=seq(1,3),
#'                             H_dim=3, R_dim=3,
#'                             k_x=0.10, k_mu=0.10, k_p=0.10,
#'                             a_1=2, a_2=2.5 )
#' 
#' @useDynLib DynMultiNet
#' 
#' @import mvtnorm
#' @importFrom MCMCpack rinvgamma
#' @importFrom stats plogis rbinom rgamma runif
#' @export
#' 

gen_synth_net <- function( nodes_net,
                           times_net,
                           layers_net,
                           directed=TRUE, weighted=TRUE,
                           H_dim=3, R_dim=3,
                           add_eff_weight=TRUE, add_eff_link=TRUE,
                           nu_sigma=1,tau_sigma=1,
                           nu_delta=1,tau_delta=1,
                           x_ijtkp=NULL,
                           rds_file=NULL ) {
  
  V_net <- length(nodes_net)
  T_net <- length(times_net)
  K_net <- length(layers_net)
  P_pred <- NULL
  if(!is.null(x_ijtkp)){P_pred<-dim(x_ijtkp)[5]}
  
  ### NETWORK ###
  y_ijtk <- array( data=NA,
                   dim=c(V_net,V_net,T_net,K_net) )
  dimnames(y_ijtk) <- list(nodes_net,nodes_net,times_net,layers_net)
  # indicator y_ijtk>0
  z_ijtk <- y_ijtk
  
  # Generating model parameters #
  # sigma_k <- 1/rgamma(4+P_pred,shape=nu_sigma,scale=tau_sigma)
  sigma_k <- MCMCpack::rinvgamma(K_net,shape=nu_sigma,scale=tau_sigma)
  
  delta <- setNames( MCMCpack::rinvgamma(4+P_pred,shape=nu_delta,scale=tau_delta),
                     paste("delta",c("mu","lambda","s","p",paste("beta",1:P_pred,sep="_")),sep="_") )
  if(!directed){delta["delta_p"]=delta["delta_s"]}
  
  ### Baseline parameter ###
  # at time t for layer k
  
  # Prior Covariance matrix for GPs
  GP_cov_prior <- vector("list",4+P_pred)
  names(GP_cov_prior) <- names(delta)
  for(i in 1:(4+P_pred)) { #i<-1
    GP_cov_prior[[i]] <- outer( times_net, times_net, FUN=function(x,y,k=delta[i]){ exp(-((x-y)/k)^2) } )
  }
  
  ## then we do:
  # GP_cov_prior_mu = GP_cov_prior^(-1/delta["delta_mu"]^2)
  ## ...and so on
  
  ### Start: Dynamics ###
  
  ## Baseline process ##
  theta_tk <- eta_tk <- matrix( NA, nrow=T_net, ncol=K_net )
  for(k in 1:K_net){ # k<-1
    theta_tk[,k] <- mvtnorm::rmvnorm( n = 1,
                                      mean = rep(0,T_net),
                                      sigma = GP_cov_prior[["delta_mu"]] )
    eta_tk[,k] <- mvtnorm::rmvnorm( n = 1,
                                    mean = rep(0,T_net),
                                    sigma = GP_cov_prior[["delta_lambda"]] )
  }
  
  ## Global Latent coordinates ##
  if(!directed){
    uv_ith_shared <- ab_ith_shared <- array( NA, dim=c(V_net,T_net,H_dim) )
    for( i in 1:V_net){
      for( h in 1:H_dim){
        uv_ith_shared[i,,h] <- mvtnorm::rmvnorm( n = 1,
                                                 mean = rep(0,T_net),
                                                 sigma = GP_cov_prior[["delta_mu"]] )
        ab_ith_shared[i,,h] <- mvtnorm::rmvnorm( n = 1,
                                                 mean = rep(0,T_net),
                                                 sigma = GP_cov_prior[["delta_lambda"]] )
      }
    }
  } else {
    uv_ith_shared <- ab_ith_shared <- list( sender=array( NA, dim=c(V_net,T_net,H_dim) ),
                                            receiver=array( NA, dim=c(V_net,T_net,H_dim) ) )
    for( dir_i in 1:2 ) {
      for( i in 1:V_net){
        for( h in 1:H_dim){
          uv_ith_shared[[dir_i]][i,,h] <- mvtnorm::rmvnorm( n = 1,
                                                            mean = rep(0,T_net),
                                                            sigma = GP_cov_prior[["delta_mu"]] )
          ab_ith_shared[[dir_i]][i,,h] <- mvtnorm::rmvnorm( n = 1,
                                                            mean = rep(0,T_net),
                                                            sigma = GP_cov_prior[["delta_lambda"]] )
        }
      }
    }
  }
  ## Layer-specific Latent coordinates ##
  if( K_net>1 ){
    if(!directed){
      uv_ithk <- ab_ithk <- array( NA, dim=c(V_net,T_net,R_dim,K_net) )
      for( k in 1:K_net){
        for( i in 1:V_net){
          for( h in 1:R_dim){
            # Link strength #
            uv_ithk[i,,h,k] <- mvtnorm::rmvnorm( n = 1,
                                                 mean = matrix(0,T_net,1),
                                                 sigma = GP_cov_prior[["delta_mu"]] )
            # Link incidence #
            ab_ithk[i,,h,k] <- mvtnorm::rmvnorm( n = 1,
                                                 mean = matrix(0,T_net,1),
                                                 sigma = GP_cov_prior[["delta_lambda"]] )
          }
        }
      }
    } else {
      uv_ithk <- ab_ithk <- list( sender=array( NA, dim=c(V_net,T_net,R_dim,K_net) ),
                                  receiver=array( NA, dim=c(V_net,T_net,R_dim,K_net) ) )
      for( dir_i in 1:2){
        for( k in 1:K_net){
          for( i in 1:V_net){
            for( h in 1:H_dim){
              # Link strength #
              uv_ithk[[dir_i]][i,,h,k] <- mvtnorm::rmvnorm( n = 1,
                                                            mean = matrix(0,T_net,1),
                                                            sigma = GP_cov_prior[["delta_mu"]] )
              # Link incidence #
              ab_ithk[[dir_i]][i,,h,k] <- mvtnorm::rmvnorm( n = 1,
                                                            mean = matrix(0,T_net,1),
                                                            sigma = GP_cov_prior[["delta_lambda"]] )
            }
          }
        }
      }
    }
  } else {
    uv_ithk <- ab_ithk <- NULL
  }
  
  ## Additive effects ##
  if(!directed){
    sp_link_it_shared <- sp_weight_it_shared <- array(NA,dim=c(V_net,T_net))
    for( i in 1:V_net){
      # Link incidence #
      sp_link_it_shared[i,] <- mvtnorm::rmvnorm( n = 1,
                                                 mean = rep(0,T_net),
                                                 sigma = GP_cov_prior[["delta_s"]] )
      # Link strength #
      sp_weight_it_shared[i,] <- mvtnorm::rmvnorm( n = 1,
                                                   mean = rep(0,T_net),
                                                   sigma = GP_cov_prior[["delta_s"]] )
      
    }
  } else {
    sp_link_it_shared <- sp_weight_it_shared <- array(NA,dim=c(V_net,T_net,2))
    for( dir_i in 1:2){
      for( i in 1:V_net){
        # Link incidence #
        sp_link_it_shared[i,,dir_i] <- mvtnorm::rmvnorm( n = 1,
                                                         mean = rep(0,T_net),
                                                         sigma = GP_cov_prior[[ c("delta_s","delta_p")[dir_i] ]] )
        # Link strength #
        sp_weight_it_shared[i,,dir_i] <- mvtnorm::rmvnorm( n = 1,
                                                           mean = rep(0,T_net),
                                                           sigma = GP_cov_prior[[ c("delta_s","delta_p")[dir_i] ]] )
      }
    }
  }
  if( !add_eff_link ){ sp_link_it_shared <- NULL }
  if( !add_eff_weight ){ sp_weight_it_shared <- NULL }
  
  ## External effects ##
  if(!is.null(x_ijtkp)) {
    beta_mu_tp <- beta_lambda_tp <- matrix( NA, nrow=T_net, ncol=P_pred )
    for(l in 1:P_pred){ # k<-1
      beta_mu_tp[,l] <- mvtnorm::rmvnorm( n = 1,
                                          mean = rep(0,T_net),
                                          sigma = GP_cov_prior[[ paste("delta_beta",l,sep="_") ]] )
      beta_lambda_tp[,l] <- mvtnorm::rmvnorm( n = 1,
                                           mean = rep(0,T_net),
                                           sigma = GP_cov_prior[[ paste("delta_beta",l,sep="_") ]] )
    }
  } else {
    beta_mu_tp <- beta_lambda_tp <- NULL
  }
  ### End: Latent dynamics ###
  
  ### Start: Linear predictors ###
  mu_ijtk <- get_linpred_ijtk( baseline_tk=theta_tk,
                               add_eff_it_shared=sp_weight_it_shared,
                               coord_ith_shared=uv_ith_shared, coord_ithk=uv_ithk,
                               beta_edge_tp=beta_mu_tp, x_ijtkp=x_ijtkp,
                               directed=directed )
  # Make mu more insteresting
  mu_ijtk[] <- scale(c(mu_ijtk)) + rep(rnorm(K_net,8,1),each=V_net*V_net*T_net)
  
  gamma_ijtk <- get_linpred_ijtk( baseline_tk=eta_tk,
                                  add_eff_it_shared=sp_link_it_shared,
                                  coord_ith_shared=ab_ith_shared, coord_ithk=ab_ithk,
                                  beta_edge_tp=beta_lambda_tp, x_ijtkp=x_ijtkp,
                                  directed=directed )
  ### End: Linear predictors ###
  
  ### links probabilities ###
  pi_ijtk <- plogis(gamma_ijtk)
  for( t in 1:T_net){
    for( k in 1:K_net){
      diag(pi_ijtk[,,t,k])<-0
    }
  }
  
  ### link incidence ###
  z_ijtk[] <- rbinom(n=length(c(gamma_ijtk)),size=1,prob=c(pi_ijtk))
  ### link weight, assuming there's a link ###
  y_ijtk[] <- rnorm(n=length(c(mu_ijtk)),mean=c(mu_ijtk),sd=rep(sigma_k,each=V_net*V_net*T_net))
  
  low_tri <-lower.tri(y_ijtk[,,t,k])
  
  if(!directed){
    for( t in 1:T_net) {
      for( k in 1:K_net) { #t<-1;k<-1
        z_ijtk[,,t,k][t(low_tri)] <- z_ijtk[,,t,k][low_tri]
        y_ijtk[,,t,k][t(low_tri)] <- y_ijtk[,,t,k][low_tri]
      }
    }
  }
  
  ### Weights ###
  if(!weighted){
    y_ijtk = z_ijtk
  } else {
    y_ijtk = z_ijtk * y_ijtk
  }
  
  synth_net <- list( y_ijtk = y_ijtk,
                     
                     x_ijtkp = x_ijtkp,
                     
                     nodes_net=nodes_net,
                     times_net=times_net,
                     layers_net=layers_net,
                     directed=directed,
                     weighted=weighted,
                     H_dim=H_dim,
                     R_dim=R_dim,
                     
                     mu_ijtk = mu_ijtk,
                     sigma_k=sigma_k,
                     theta_tk = theta_tk,
                     uv_ith_shared=uv_ith_shared,
                     uv_ithk=uv_ithk,
                     sp_weight_it_shared=sp_weight_it_shared,
                     beta_mu_tp=beta_mu_tp,
                     
                     pi_ijtk = pi_ijtk,
                     
                     gamma_ijtk = gamma_ijtk,
                     eta_tk = eta_tk,
                     ab_ith_shared=ab_ith_shared,
                     ab_ithk=ab_ithk,
                     sp_link_it_shared=sp_link_it_shared,
                     beta_lambda_tp=beta_lambda_tp )
  
  if(!is.null(rds_file)){
    if( substr(rds_file,nchar(rds_file)-3,nchar(rds_file))!=".rds" ) {rds_file<-paste(rds_file,".rds",sep="")}
    saveRDS( synth_net, file=rds_file )
  }
  return( synth_net )
}
