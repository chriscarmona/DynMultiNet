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
#' @importFrom foreach foreach
#' @importFrom pROC auc
#' @import ggplot2
#' @export
#' 

plot_loss <- function( x,
                       y_ijtk_complete=NULL,
                       FUN=median ) {
  
  V_net <- dim(x$y_ijtk)[1]
  T_net <- dim(x$y_ijtk)[3]
  K_net <- dim(x$y_ijtk)[4]
  
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
  
  mu_hat <- mu <- 0
  
  p_hat <- apply( x$pi_ijtk_mcmc, 1:4, FUN )
  if(x$weighted){
    mu_hat <- apply( x$mu_ijtk_mcmc, 1:4, FUN )
    E_y_hat <- apply( x$mu_ijtk_mcmc * x$pi_ijtk_mcmc, 1:4, FUN )
  }
  
  loss_table <- foreach::foreach( k = 1:K_net, .combine=rbind ) %:%
    foreach::foreach( t = 1:T_net, .combine=rbind ) %dopar% {
      if( is.null(y_ijtk_complete) ) {
        y <- x$y_ijtk[,,t,k]
      } else {
        y <- y_ijtk_complete[,,t,k]
      }
      z <- c(y>0)*1
      p <- c(p_hat[,,t,k])
      E_y <- c(E_y_hat[,,t,k])
      
      if(x$weighted){ mu <- c(mu_hat[,,t,k]) }
      auc <- ifelse( all(is.na(z)) | all(is.na(p)) ,
                     NA,
                     pROC::auc( response=z, predictor=p ) )
      rmse_cond <- ifelse( all(is.na(y)) | all(is.na(mu)) ,
                           NA,
                           sqrt(mean((y[!is.na(y)&y>0]-mu[!is.na(y)&y>0])^2)) )
      rmse <- ifelse( all(is.na(y)) | all(is.na(E_y)) ,
                      NA,
                      sqrt(mean((y[!is.na(y)]-E_y[!is.na(y)])^2)) )
      c(t,k,auc,rmse_cond,rmse)
    }
  rownames(loss_table) <- NULL
  colnames(loss_table) <- c("time","layer","AUC","RMSE_cond","RMSE")
  
  loss_table <- loss_table %>%
    as.data.frame() %>%
    dplyr::mutate( time=x$time_all[time],
                   layer=x$layer_all[layer] ) %>%
    gather("criteria","value",-c(time,layer))
    
  auc_plot <- loss_table %>%
    filter(criteria=="AUC") %>%
    ggplot() +
    theme_bw() +
    geom_line( aes(x=time,y=value,color=layer) ) +
    labs(y="AUC") +
    coord_cartesian(y=c(0.5,1))
  
  rmse_cond_plot <- loss_table %>%
    filter(criteria=="RMSE_cond") %>%
    ggplot() +
    theme_bw() +
    geom_line( aes(x=time,y=value,color=layer) ) +
    labs(y="RMSE | y>0")
  
  rmse_plot <- loss_table %>%
    filter(criteria=="RMSE") %>%
    ggplot() +
    theme_bw() +
    geom_line( aes(x=time,y=value,color=layer) ) +
    labs(y="RMSE")
  
  rmse_both_plot <- loss_table %>%
    filter(is.element(criteria,c("RMSE_cond","RMSE"))) %>%
    ggplot() +
    theme_bw() +
    geom_line( aes(x=time,y=value,color=layer,lty=criteria) ) +
    labs(y="RMSE",lty="")
    
  
  return( list( loss_table=loss_table,
                auc_plot=auc_plot,
                rmse_cond_plot=rmse_cond_plot,
                rmse_plot=rmse_plot,
                rmse_both_plot=rmse_both_plot ) )
}
