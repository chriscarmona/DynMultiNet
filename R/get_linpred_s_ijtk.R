#' @title
#'    Bayesian Learning of Dynamic Multilayer Networks with binary data
#'
#' @description
#'    \code{DynMultiNet_bin} Implements model from Durante and Dunson, 2018
#'
#' @param y_ijtk Array of dimension dim(y_ijtk)=c(V_net,V_net,T_net,K_net). Binary undirected links between nodes i and j, at time t, on layer k
#' @param mu_tk Matrix. Edge specific baseline factor of the linear predictor
#' @param x_iht Array. Latent coordinates of element i in dimension h
#' @param z_tp Matrix. Global predictors data, predictor p at time t
#' @param z_tkp Array. Layer specific predictors data, predictor p at time t for layer k
#' @param z_ijtkp Array. Edge specific predictors data, predictor p at time t for edge bewtween i and j in layer k
#' @param beta_z_tp Column Matrix. Coefficients associated with z_tp
#' @param beta_z_tkp Column Matrix. Coefficients associated with z_tkp
#' @param beta_z_ijtkp Column Matrix. Coefficients associated with z_ijtkp
#' @param pred_id_tp Matrix. List of global predictors
#' @param pred_id_tkp Matrix. List of layer specific predictors
#' @param pred_id_ijtkp Matrix. List of edge specific predictors
#' 
#' @details
#'    Linear predictor of response variable y_ijtk under a logistic link,
#'    logit( E(y_ijtk) ) = s_ijtk
#' 
#' @return
#'    An array s_ijtk of dimension dim(s_ijtk)=c(V_net,V_net,T_net,K_net)
#' 
#' @keywords internal
#' 

get_linpred_s_ijtk <- function( y_ijtk, mu_tk, x_iht,
                                pred_all=NULL, layer_all=NULL,
                                z_tp=NULL, z_tkp=NULL, z_ijtkp=NULL,
                                beta_z_tp=NULL, beta_z_tkp=NULL, beta_z_ijtkp=NULL,
                                pred_id_tp=NULL, pred_id_tkp=NULL, pred_id_ijtkp=NULL ) {
  
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  s_ijtk <- array( data=NA,
                   dim=c(V_net,V_net,T_net,K_net) )
  
  for( k in 1:K_net ){
    for( t in 1:T_net ){
      s_ijtk[,,t,k] <- mu_tk[t,k] + x_iht[,,t] %*% t(x_iht[,,t])
    }
  }
  rm(k,t)
  
  # Global predictors
  if(!is.null(z_tp) & !is.null(beta_z_tp) & !is.null(pred_id_tp) ){
    for( row_i in 1:nrow(pred_id_tp)) {
      p <- match(pred_id_tp[row_i,"id"],pred_all)
      for(t in 1:T_net){
        s_ijtk[,,t,] <- s_ijtk[,,t,] + z_tp[t,p] * beta_z_tp[t,row_i]
      }
    }
    rm(row_i,p,t)
  }
  # Layer specific predictors
  if(!is.null(z_tkp) & !is.null(beta_z_tkp) & !is.null(pred_id_tkp) ){
    for( row_i in 1:nrow(pred_id_tkp)) {
      k <- match(pred_id_tkp[row_i,"layer"],layer_all)
      p <- match(pred_id_tkp[row_i,"id"],pred_all)
      for(t in 1:T_net){
        s_ijtk[,,t,k] <- s_ijtk[,,t,k] + z_tkp[t,k,p] * beta_z_tkp[t,row_i]
      }
    }
    rm(row_i,k,p,t)
  }
  # Edge specific predictors
  if( !is.null(z_ijtkp) & !is.null(beta_z_ijtkp) & !is.null(pred_id_ijtkp) ){
    for( row_i in 1:nrow(pred_id_ijtkp)) {
      k <- match(pred_id_ijtkp[row_i,"layer"],layer_all)
      p <- match(pred_id_ijtkp[row_i,"id"],pred_all)
      for(t in 1:T_net){
        s_ijtk[,,t,k] <- s_ijtk[,,t,k] + z_ijtkp[,,t,k,p] * beta_z_ijtkp[t,row_i]
      }
    }
    rm(row_i,k,p,t)
  }
  
  s_ijtk[is.na(y_ijtk)]<-NA
  
  if(any(is.na(s_ijtk[!is.na(y_ijtk)]))){stop('There is a problem creating "s_ijtk" (perhaps related to predictors)')}
  
  return(s_ijtk)
}