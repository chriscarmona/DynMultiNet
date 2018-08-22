#' @title
#'    Bayesian Learning of Dynamic Multilayer Networks
#'
#' @description
#'    \code{DynMultiNet_bin} Implements model from Durante and Dunson, 2018
#'
#' @param stan_fit Object of class "stanfit".
#' @param directed Boolean. Indicates if the network is directed.
#' @param weighted Boolean. Indicates if the network is weighted.
#' 
#' @details
#'    Transform the output from the stan implemetation of the Latent Network model into a dmn_mcmc class object.
#' 
#' @return
#'    An object of class \code{dmn_mcmc}.
#' 
#' @importFrom rstan extract
#' 
#' @keywords internal
#' 

dmn_mcmc_from_stan <- function( stan_fit,
                                directed=FALSE,
                                weighted=FALSE ) {
  
  posterior <- rstan::extract(stan_fit)
  
  if( !directed & !weighted ) {
    dmn_mcmc <- list( y_ijtk=NULL,
                      
                      directed = directed,
                      weighted = weighted,
                      
                      n_chains_mcmc=NULL,
                      n_iter_mcmc = NULL, n_burn = NULL, n_thin = NULL,
                      
                      a_1=NULL, a_2=NULL,
                      k_mu=NULL, k_x=NULL,
                      
                      pi_ijtk_mcmc = aperm(posterior[["pi_ij_tk"]],c(4,5,2,3,1)),
                      
                      node_all=NULL, time_all=NULL, layer_all=NULL,
                      time_all_idx_net=NULL,
                      
                      mu_tk_mcmc = posterior[["mu_tk"]],
                      
                      x_ith_shared_mcmc = aperm(posterior[["x_ti_h_shared"]],c(4,3,2,1)),
                      x_ithk_mcmc = aperm(posterior[["x_ti_hk"]],c(5,4,2,3,1)),
                      tau_h_shared_mcmc = aperm(posterior[["tau_h_shared"]],c(2,1)),
                      tau_h_k_mcmc = aperm(posterior[["tau_hk"]],c(2,3,1)),
                      
                      pred_id_layer=NULL, pred_id_edge=NULL,
                      beta_z_layer_mcmc=NULL,
                      beta_z_edge_mcmc=NULL )
  } else if( directed & !weighted ) {
    dmn_mcmc <- list( y_ijtk=NULL,
                      
                      directed = directed,
                      weighted = weighted,
                      
                      n_chains_mcmc=NULL,
                      n_iter_mcmc = NULL, n_burn = NULL, n_thin = NULL,
                      time_all_idx_net=NULL,
                      
                      a_1=NULL, a_2=NULL,
                      k_mu=NULL, k_x=NULL,
                      
                      pi_ijtk_mcmc = aperm(posterior[["pi_ij_tk"]],c(4,5,2,3,1)),
                      
                      node_all=NULL, time_all=NULL, layer_all=NULL,
                      
                      mu_tk_mcmc = aperm(posterior[["mu_tk"]],c(2,3,1)),
                      
                      x_ith_shared_mcmc = list( aperm(posterior[["x_ti_h_shared"]][,1,,,],c(4,3,2,1)),
                                                aperm(posterior[["x_ti_h_shared"]][,2,,,],c(4,3,2,1)) ),
                      x_ithk_mcmc = list( aperm(posterior[["x_ti_hk"]][,1,,,,],c(5,4,2,3,1)),
                                          aperm(posterior[["x_ti_hk"]][,2,,,,],c(5,4,2,3,1)) ),
                      tau_h_shared_mcmc = list( aperm(posterior[["tau_h_shared"]][,1,],c(2,1) ),
                                                aperm(posterior[["tau_h_shared"]][,1,],c(2,1)) ),
                      tau_h_k_mcmc = list( aperm(posterior[["tau_hk"]][,1,,],c(2,3,1)),
                                           aperm(posterior[["tau_hk"]][,1,,],c(2,3,1)) ),
                      
                      pred_id_layer=NULL, pred_id_edge=NULL,
                      beta_z_layer_mcmc=NULL,
                      beta_z_edge_mcmc=NULL )
  } else {
    stop("Not supported.")
  }
  dmn_mcmc <- structure( dmn_mcmc, class="dmn_mcmc" )
  
  return( dmn_mcmc )
}
