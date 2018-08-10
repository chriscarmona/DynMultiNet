#' @title
#'    Bayesian Learning of Dynamic Multilayer Networks with binary data
#'
#' @description
#'    \code{DynMultiNet_bin} Implements model from Durante and Dunson, 2018
#'
#' @param y_ijtk Array of dimension dim(y_ijtk)=c(V_net,V_net,T_net,K_net). Binary undirected links between nodes i and j, at time t, on layer k
#' @param mu_tk Matrix. Edge specific baseline factor of the linear predictor
#' @param x_ith_shared Array. Global latent coordinates of element i in dimension h
#' @param x_ithk Array. Layer Specific Latent coordinates of element i in dimension h. Only for multilayer networks.
#' @param pred_all Vector. Names of all predictors.
#' @param layer_all Vector. Names of all layers.
#' @param z_tkp Array. Layer specific predictors data, predictor p at time t for layer k
#' @param z_ijtkp Array. Edge specific predictors data, predictor p at time t for edge bewtween i and j in layer k
#' @param beta_z_layer Column Matrix. Coefficients associated with z_tkp
#' @param beta_z_edge Column Matrix. Coefficients associated with z_ijtkp
#' @param pred_id_layer Matrix. List of layer specific predictors
#' @param pred_id_edge Matrix. List of edge specific predictors
#' @param directed Boolean. Indicates if the provided network is directed, i.e. the adjacency matrix is assymetrical.
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

get_linpred_s_ijtk <- function( y_ijtk, mu_tk,
                                x_ith_shared, x_ithk=NULL,
                                pred_all=NULL, layer_all=NULL,
                                z_tkp=NULL, z_ijtkp=NULL,
                                beta_z_layer=NULL, beta_z_edge=NULL,
                                pred_id_layer=NULL, pred_id_edge=NULL,
                                directed=FALSE ) {
  
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  s_ijtk <- array( data=NA,
                   dim=c(V_net,V_net,T_net,K_net) )
  
  # Baseline Process and Global latent coordinates
  if( directed ) {
    for( k in 1:K_net ){ # k<-1
      for( t in 1:T_net ){ # t<-1
        s_ijtk[,,t,k] <- mu_tk[t,k] + x_ith_shared[[1]][,t,] %*% t(x_ith_shared[[2]][,t,])
      }
    }; rm(k,t)
  } else {
    for( k in 1:K_net ){
      for( t in 1:T_net ){
        s_ijtk[,,t,k] <- mu_tk[t,k] + x_ith_shared[,t,] %*% t(x_ith_shared[,t,])
      }
    }; rm(k,t)
  }
  
  # Layer-specific latent coordinates
  if(!is.null(x_ithk)) {
    if( directed ) {
      for( k in 1:K_net ){
        for( t in 1:T_net ){
          s_ijtk[,,t,k] <- s_ijtk[,,t,k] + x_ithk[[1]][,t,,k] %*% t(x_ithk[[2]][,t,,k])
        }
      }; rm(k,t)
    } else {
      for( k in 1:K_net ){
        for( t in 1:T_net ){
          s_ijtk[,,t,k] <- s_ijtk[,,t,k] + x_ithk[,t,,k] %*% t(x_ithk[,t,,k])
        }
      }; rm(k,t)
    }
  }
  
  # Layer specific predictors
  if(!is.null(z_tkp) & !is.null(beta_z_layer) & !is.null(pred_id_layer) ){
    for( row_i in 1:nrow(pred_id_layer)) { # row_i<-1
      k <- match(pred_id_layer[row_i,"layer"],layer_all)
      p <- match(pred_id_layer[row_i,"id"],pred_all)
      for(t in 1:T_net){ # t<-1
        s_ijtk[,,t,k] <- s_ijtk[,,t,k] + z_tkp[t,k,p] * beta_z_layer[t,row_i]
      }
    }
    rm(row_i,k,p,t)
  }
  if(any(is.na(s_ijtk[!is.na(y_ijtk)]))){stop('There is a problem creating "s_ijtk" (perhaps related to layer predictors)')}
  
  # Edge specific predictors
  if( !is.null(z_ijtkp) & !is.null(beta_z_edge) & !is.null(pred_id_edge) ){
    for( row_i in 1:nrow(pred_id_edge)) {
      k <- match(pred_id_edge[row_i,"layer"],layer_all)
      p <- match(pred_id_edge[row_i,"id"],pred_all)
      for(t in 1:T_net){
        s_ijtk[,,t,k] <- s_ijtk[,,t,k] + z_ijtkp[,,t,k,p] * beta_z_edge[t,row_i]
      }
    }
    rm(row_i,k,p,t)
  }
  
  s_ijtk[is.na(y_ijtk)]<-NA
  
  if(any(is.na(s_ijtk[!is.na(y_ijtk)]))){stop('There is a problem creating "s_ijtk" (perhaps related to edge predictors)')}
  
  return(s_ijtk)
}
