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
#' @importFrom reshape melt
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
  if(any(is.na(x$pi_ijtk_mcmc[!diag_y_idx]))){stop("There are NA values in pi_ijtk_mcmc")}
  
  ### start: WAIC ###
  loglik_y_link_pi_post <- dbinom(1*(x$y_ijtk>0),1,x$pi_ijtk_mcmc,log=T)
  
  # eliminates -Inf
  loglik_y_link_pi_post[is.infinite(loglik_y_link_pi_post)] <- -1000
  
  aux_log_f <- matrix( loglik_y_link_pi_post,
                       nrow = prod(dim(loglik_y_link_pi_post)[1:4]),
                       ncol = dim(loglik_y_link_pi_post)[5] )
  aux_log_f <- t(aux_log_f)
  if(!all.equal(is.na(c(x$y_ijtk)), apply(is.na(aux_log_f),2,all) )){stop("Problem calculating WAIC")}
  aux_log_f <- aux_log_f[,!is.na(c(x$y_ijtk))]
  
  dmn_waic_link <- loo::waic(aux_log_f)
  
  ### initialize output table ###
  out_table <- reshape::melt(dmn_waic_link$estimates)
  colnames(out_table)[1:2] <- c("criteria","estim_type")
  out_table$response <- 'binary'
  
  # Add log pointwise predictive density
  out_table_temp <- out_table[1,,drop=F];out_table_temp[]<-NA
  out_table_temp$criteria <- "lpd_mcmc"
  out_table_temp$estim_type <- "Estimate"
  out_table_temp$value <- sum( log( apply( exp(aux_log_f), 2, mean ) ) )
  out_table_temp$response <- 'binary'
  out_table <- rbind(out_table,out_table_temp)
  
  if(x$weighted){
    # Array with MCMC for sigma
    aux_sigma_k <- x$mu_ijtk_mcmc; aux_sigma_k[] <- NA
    aux_sigma_k <- aperm(aux_sigma_k, c(4,5,1,2,3))
    aux_sigma_k[] <- x$sigma_k_mcmc
    aux_sigma_k <- aperm(aux_sigma_k, c(3,4,5,1,2))
    dim(aux_sigma_k)
    
    # Gaussian probability for weight
    loglik_y_weight_pi_post <- dnorm(x$y_ijtk,mean=x$mu_ijtk_mcmc,sd=aux_sigma_k,log=T)
    # eliminates -Inf
    loglik_y_weight_pi_post[is.infinite(loglik_y_weight_pi_post)] <- -1000
    
    aux_log_f_w <- matrix( loglik_y_weight_pi_post,
                           nrow = prod(dim(loglik_y_weight_pi_post)[1:4]),
                           ncol = dim(loglik_y_weight_pi_post)[5] )
    aux_log_f_w <- t(aux_log_f_w)
    if(!all.equal(is.na(c(x$y_ijtk)), apply(is.na(aux_log_f_w),2,all) )){stop("Problem calculating WAIC")}
    
    # dirac delta for non-linked pairs
    aux_log_f_w[,which(c(x$y_ijtk)==0)] <- 0
    
    aux_log_f_w <- aux_log_f_w[,!is.na(c(x$y_ijtk))]
    
    dmn_waic_weight <- loo::waic(aux_log_f+aux_log_f_w)
    
    ### initialize output table ###
    out_table_temp <- reshape::melt(dmn_waic_weight$estimates)
    colnames(out_table_temp)[1:2] <- c("criteria","estim_type")
    out_table_temp$response <- 'weighted'
    out_table <- rbind(out_table,out_table_temp)
    
    # Add log pointwise predictive density
    out_table_temp <- out_table[1,,drop=F];out_table_temp[]<-NA
    out_table_temp$criteria <- "lpd_mcmc"
    out_table_temp$estim_type <- "Estimate"
    out_table_temp$value <- sum( log( apply( exp(aux_log_f+aux_log_f_w), 2, mean ) ) )
    out_table_temp$response <- 'weighted'
    out_table <- rbind(out_table,out_table_temp)
  }
  ### end: WAIC ###
  
  ### start: Predictive out-of-sample ###
  if(!is.null(y_ijtk_complete)){
    
    loglik_y_link_pi_post <- dbinom(1*(y_ijtk_complete>0),1,x$pi_ijtk_mcmc,log=T)
    
    aux_oos_log_f <- matrix( loglik_y_link_pi_post,
                             nrow = prod(dim(loglik_y_link_pi_post)[1:4]),
                             ncol = dim(loglik_y_link_pi_post)[5] )
    aux_oos_log_f <- t(aux_oos_log_f)
    aux_oos_log_f <- aux_oos_log_f[,which(y_oos==1)]
    
    # Add log pointwise predictive density, out of sample
    out_table_temp <- out_table[1,,drop=F];out_table_temp[]<-NA
    out_table_temp$criteria <- "lpd_mcmc_oos"
    out_table_temp$estim_type <- "Estimate"
    out_table_temp$value <- sum( log( apply( exp(aux_oos_log_f), 2, mean ) ) )
    out_table_temp$response <- 'binary'
    out_table <- rbind(out_table,out_table_temp)
    
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
      
      # Add log pointwise predictive density, out of sample
      out_table_temp <- out_table[1,,drop=F];out_table_temp[]<-NA
      out_table_temp$criteria <- "lpd_mcmc_oos"
      out_table_temp$estim_type <- "Estimate"
      out_table_temp$value <- sum( log( apply( exp(aux_oos_log_f+aux_oos_log_f_w), 2, mean ) ) )
      out_table_temp$response <- 'weighted'
      out_table <- rbind(out_table,out_table_temp)
    }
  }
  ### end: Predictive out-of-sample ###
  
  return( out_table )
}
