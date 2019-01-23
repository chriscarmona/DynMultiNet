#' @title
#'    Bayesian Learning of Dynamic Multilayer Networks with binary data
#'
#' @description
#'    \code{DynMultiNet_bin} Implements model from Durante and Dunson, 2018
#'
#' @param baseline_tk Matrix. Edge specific baseline factor of the linear predictor
#' @param add_eff_it_shared Array. Global additive effects of agent i at time t
#' @param add_eff_itk Array. Layer Specific additive effects of agent i at time t. Only for multilayer networks.
#' @param coord_ith_shared Array. Global latent coordinates of element i in dimension h
#' @param coord_ithk Array. Layer Specific Latent coordinates of element i in dimension h. Only for multilayer networks.
#' @param beta_edge_tp Array. Edge Specific coefficients associated with external covariates x_ijtkp.
#' @param x_ijtkp Array. Edge Specific external covariates.
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
                              add_eff_it_shared=NULL, add_eff_itk=NULL,
                              coord_ith_shared, coord_ithk=NULL,
                              beta_edge_tp=NULL, x_ijtkp=NULL,
                              directed=FALSE ) {
  
  T_net <- dim(baseline_tk)[1]
  K_net <- dim(baseline_tk)[2]
  if(!directed){
    V_net <- dim(coord_ith_shared)[1]
  } else {
    V_net <- dim(coord_ith_shared[[1]])[1]
  }
  
  if(!is.null(beta_edge_tp)){
    P_pred <- dim(beta_edge_tp)[2]
    
    if( is.null(x_ijtkp) | (!all.equal(dim(x_ijtkp),c(V_net,V_net,T_net,K_net,P_pred))) ){
      stop('There is a problem with "beta_edge_tp" and "x_ijtkp"')
    }
  }
  
  # Initialise linpred
  linpred_ijtk <- array( data=0,
                         dim=c(V_net,V_net,T_net,K_net) )
  
  # Layer baseline process
  for( k in 1:K_net ){ # k<-1
    for( t in 1:T_net ){ # t<-1
      linpred_ijtk[,,t,k] <- linpred_ijtk[,,t,k] + baseline_tk[t,k]
    }
  }
  
  # Global additive effects
  if(!is.null(add_eff_it_shared)) {
    if( !directed ) {
      for( k in 1:K_net ){
        for( t in 1:T_net ){
          linpred_ijtk[,,t,k] <- linpred_ijtk[,,t,k] + matrix(add_eff_it_shared[,t],V_net,V_net) + matrix(add_eff_it_shared[,t],V_net,V_net,byrow=T)
        }
      }
    } else {
      for( k in 1:K_net ){ # k<-1
        for( t in 1:T_net ){ # t<-1
          linpred_ijtk[,,t,k] <- linpred_ijtk[,,t,k] + matrix(add_eff_it_shared[,t,1],V_net,V_net) + matrix(add_eff_it_shared[,t,2],V_net,V_net,byrow=T)
        }
      }
    }
  }
  
  # Layer-specific additive effects
  if(!is.null(add_eff_itk)) {
    if( !directed ) {
      for( k in 1:K_net ){
        for( t in 1:T_net ){
          linpred_ijtk[,,t,k] <- linpred_ijtk[,,t,k] + matrix(add_eff_itk[,t,k],V_net,V_net) + matrix(add_eff_itk[,t,k],V_net,V_net,byrow=T)
        }
      }
    } else {
      for( k in 1:K_net ){ # k<-1
        for( t in 1:T_net ){ # t<-1
          linpred_ijtk[,,t,k] <- linpred_ijtk[,,t,k] + matrix(add_eff_itk[,t,k,1],V_net,V_net) + matrix(add_eff_itk[,t,k,2],V_net,V_net,byrow=T)
        }
      }
    }
  }
  
  # Global latent coordinates
  if( !directed ) {
    for( k in 1:K_net ){
      for( t in 1:T_net ){
        linpred_ijtk[,,t,k] <- linpred_ijtk[,,t,k] + coord_ith_shared[,t,] %*% t(coord_ith_shared[,t,])
      }
    }
  } else {
    for( k in 1:K_net ){ # k<-1
      for( t in 1:T_net ){ # t<-1
        linpred_ijtk[,,t,k] <- linpred_ijtk[,,t,k] + coord_ith_shared[[1]][,t,] %*% t(coord_ith_shared[[2]][,t,])
      }
    }
  }
  
  # Layer-specific latent coordinates
  if(!is.null(coord_ithk)) {
    if( !directed ) {
      for( k in 1:K_net ){
        for( t in 1:T_net ){
          linpred_ijtk[,,t,k] <- linpred_ijtk[,,t,k] + coord_ithk[,t,,k] %*% t(coord_ithk[,t,,k])
        }
      }
    } else {
      for( k in 1:K_net ){
        for( t in 1:T_net ){
          linpred_ijtk[,,t,k] <- linpred_ijtk[,,t,k] + coord_ithk[[1]][,t,,k] %*% t(coord_ithk[[2]][,t,,k])
        }
      }
    }
  }
  
  if( !is.null(beta_edge_tp) ){
    for( l in 1:P_pred ){
      for( t in 1:T_net ){
        linpred_ijtk[,,t,] <- linpred_ijtk[,,t,] + beta_edge_tp[t,l] * x_ijtkp[,,t,,l]
      }
    }
  }
  
  return(linpred_ijtk)
}
