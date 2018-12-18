#' @title
#'    Bayesian Learning of Dynamic Multilayer Networks with binary data
#'
#' @description
#'    \code{DynMultiNet_bin} Implements model from Durante and Dunson, 2018
#'
#' @param baseline_tk Matrix. Edge specific baseline factor of the linear predictor
#' @param coord_ith Array. Global latent coordinates of element i in dimension h
#' @param coord_ithk Array. Layer Specific Latent coordinates of element i in dimension h. Only for multilayer networks.
#' @param directed Boolean. Indicates if the provided network is directed, i.e. the adjacency matrix is assymetrical.
#' 
#' @details
#'    Linear predictor of response variable,
#' 
#' @return
#'    An array linpred_ijtk of dimension dim(linpred_ijtk)=c(V_net,V_net,T_net,K_net)
#' 
#' @keywords internal
#' 

get_linpred_ijtk <- function( baseline_tk,
                              coord_ith, coord_ithk=NULL,
                              directed=FALSE ) {
  
  T_net <- dim(baseline_tk)[1]
  K_net <- dim(baseline_tk)[2]
  if(directed){
    V_net <- dim(coord_ith[[1]])[1]
  } else {
    V_net <- dim(coord_ith)[1]
  }
  
  linpred_ijtk <- array( data=NA,
                         dim=c(V_net,V_net,T_net,K_net) )
  
  
  # Baseline Process and Global latent coordinates
  if( directed ) {
    for( k in 1:K_net ){ # k<-1
      for( t in 1:T_net ){ # t<-1
        linpred_ijtk[,,t,k] <- baseline_tk[t,k] + coord_ith[[1]][,t,] %*% t(coord_ith[[2]][,t,])
      }
    }; rm(k,t)
  } else {
    for( k in 1:K_net ){
      for( t in 1:T_net ){
        linpred_ijtk[,,t,k] <- baseline_tk[t,k] + coord_ith[,t,] %*% t(coord_ith[,t,])
      }
    }; rm(k,t)
  }
  
  
  # Layer-specific latent coordinates
  if(!is.null(coord_ithk)) {
    if( directed ) {
      for( k in 1:K_net ){
        for( t in 1:T_net ){
          linpred_ijtk[,,t,k] <- linpred_ijtk[,,t,k] + coord_ithk[[1]][,t,,k] %*% t(coord_ithk[[2]][,t,,k])
        }
      }; rm(k,t)
    } else {
      for( k in 1:K_net ){
        for( t in 1:T_net ){
          linpred_ijtk[,,t,k] <- linpred_ijtk[,,t,k] + coord_ithk[,t,,k] %*% t(coord_ithk[,t,,k])
        }
      }; rm(k,t)
    }
  }
  
  return(linpred_ijtk)
}
