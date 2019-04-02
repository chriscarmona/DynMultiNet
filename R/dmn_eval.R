#' @title
#'    DMN Model assesment
#'
#' @description
#'    \code{dmc_eval} calcuates measures of model performance
#'
#' @param x object of class \code{dmn_mcmc}
#' @param y_ijtk_complete an array with the complete network with no missing links
#'
#' @details
#'    Calculates WAIC and AUC
#'
#' @return
#'    A \code{waic} object
#' 
#' @useDynLib DynMultiNet, .registration = TRUE
#' 
#' @importFrom loo waic
#' 
#' @export
#' 

dmn_eval <- function( x,
                      y_ijtk_complete=NULL ) {
  
  V_net <- dim(x$y_ijtk)[1]
  
  # assume no self-edges
  diag_y_idx <- matrix(FALSE,V_net,V_net); diag(diag_y_idx)<-TRUE
  diag_y_idx <- array(diag_y_idx,dim=dim(x$y_ijtk))
  x$y_ijtk[diag_y_idx] <- NA
  if(!is.null(y_ijtk_complete)){
    if(!all_equal(dim(x$y_ijtk),dim(y_ijtk_complete))){stop("The y_ijtk_complete provided is inconsistent")}
    y_ijtk_complete[diag_y_idx] <- NA
    
    # Indicator of observations considered in the out-of-sample test set
    y_oos <- y_ijtk_complete; y_oos[] <- 0
    y_oos[is.na(x$y_ijtk)&(!is.na(y_ijtk_complete))] <- 1
  }
  
  ### start: WAIC ###
  loglik_y_link_pi_post <- dbinom(1*(x$y_ijtk>0),1,x$pi_ijtk_mcmc,log=T)
  
  aux_log_f <- matrix( loglik_y_link_pi_post,
                       nrow = prod(dim(loglik_y_link_pi_post)[1:4]),
                       ncol = dim(loglik_y_link_pi_post)[5] )
  aux_log_f <- t(aux_log_f)
  if(!all.equal(is.na(c(x$y_ijtk)), apply(is.na(aux_log_f),2,all) )){stop("Problem calculating WAIC")}
  aux_log_f <- aux_log_f[,!is.na(c(x$y_ijtk))]
  dmn_waic_link <- loo::waic(aux_log_f)
  
  # log pointwise predictive density
  lppd_link <- mean( log( apply( exp(aux_log_f), 2, mean ) ) )
  
  dmn_waic_weight <- NULL
  lppd_weight <- NULL
  if(x$weighted){
    # Array with MCMC for sigma
    aux_sigma_k <- x$mu_ijtk_mcmc; aux_sigma_k[] <- NA
    aux_sigma_k <- aperm(aux_sigma_k, c(4,5,1,2,3))
    aux_sigma_k[] <- x$sigma_k_mcmc
    aux_sigma_k <- aperm(aux_sigma_k, c(3,4,5,1,2))
    dim(aux_sigma_k)
    
    # Gaussian probability for weight
    loglik_y_weight_pi_post <- dnorm(x$y_ijtk,mean=x$mu_ijtk_mcmc,sd=aux_sigma_k,log=T)
    
    aux_log_f_w <- matrix( loglik_y_weight_pi_post,
                           nrow = prod(dim(loglik_y_weight_pi_post)[1:4]),
                           ncol = dim(loglik_y_weight_pi_post)[5] )
    aux_log_f_w <- t(aux_log_f_w)
    if(!all.equal(is.na(c(x$y_ijtk)), apply(is.na(aux_log_f_w),2,all) )){stop("Problem calculating WAIC")}
    
    # dirac delta for non-linked pairs
    aux_log_f_w[,which(c(x$y_ijtk)==0)] <- 0
    
    aux_log_f_w <- aux_log_f_w[,!is.na(c(x$y_ijtk))]
    
    dmn_waic_weight <- loo::waic(aux_log_f+aux_log_f_w)
    
    # log pointwise predictive density
    lppd_weight <- mean( log( apply( exp(aux_log_f+aux_log_f_w), 2, mean ) ) )
  }
  ### end: WAIC ###
  
  ### start: Predictive out-of-sample ###
  lppd_link_oos <- NULL
  lppd_weight_oos <- NULL
  if(!is.null(y_ijtk_complete)){
    
    loglik_y_link_pi_post <- dbinom(1*(y_ijtk_complete>0),1,x$pi_ijtk_mcmc,log=T)
    
    aux_oos_log_f <- matrix( loglik_y_link_pi_post,
                             nrow = prod(dim(loglik_y_link_pi_post)[1:4]),
                             ncol = dim(loglik_y_link_pi_post)[5] )
    aux_oos_log_f <- t(aux_oos_log_f)
    aux_oos_log_f <- aux_oos_log_f[,which(y_oos==1)]
    
    # log pointwise predictive density
    lppd_link_oos <- mean( log( apply( exp(aux_oos_log_f), 2, mean ) ) )
    
    if(x$weighted){
      # Array with MCMC for sigma
      aux_sigma_k <- x$mu_ijtk_mcmc; aux_sigma_k[] <- NA
      aux_sigma_k <- aperm(aux_sigma_k, c(4,5,1,2,3))
      aux_sigma_k[] <- x$sigma_k_mcmc
      aux_sigma_k <- aperm(aux_sigma_k, c(3,4,5,1,2))
      dim(aux_sigma_k)
      
      # Gaussian probability for weight
      loglik_y_weight_pi_post <- dnorm(y_ijtk_complete,mean=x$mu_ijtk_mcmc,sd=aux_sigma_k,log=T)
      
      aux_oos_log_f_w <- matrix( loglik_y_weight_pi_post,
                                 nrow = prod(dim(loglik_y_weight_pi_post)[1:4]),
                                 ncol = dim(loglik_y_weight_pi_post)[5] )
      aux_oos_log_f_w <- t(aux_oos_log_f_w)
      
      # dirac delta for non-linked pairs
      aux_oos_log_f_w[,which(c(y_ijtk_complete)==0)] <- 0
      
      aux_oos_log_f_w <- aux_oos_log_f_w[,which(y_oos==1)]
      
      # log pointwise predictive density
      lppd_weight_oos <- mean( log( apply( exp(aux_oos_log_f+aux_oos_log_f_w), 2, mean ) ) )
    }
  }
  ### end: Predictive out-of-sample ###
  
  return( list( waic_binary=dmn_waic_link,
                waic_weighted=dmn_waic_weight,
                
                lppd_binary=lppd_link,
                lppd_weighted=lppd_weight,
                
                lppd_binary_oos=lppd_link_oos,
                lppd_weighted_oos=lppd_weight_oos ) )
}
