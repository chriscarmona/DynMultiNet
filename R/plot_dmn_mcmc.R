#' @title
#'     Plotting posterior MCMC results for the function \code{dmn_sampling}
#'
#' @description
#'     Plotting method for objects inheriting from class "\code{dmn_mcmc}".
#'
#' @param x an object of class "\code{plot_dmn_mcmc}"
#' @param param Character. Specifies the parameter that will be plotted. Possible parameters are:
#'        \describe{
#'            \item{\code{"pi_ijtk"}}{Edge probabilities between node_i and node_j at time t in layer_k.}
#'            \item{\code{"mu_tk"}}{Baseline process for edge probability at layer_k.}
#'            \item{\code{"x_ith_shared"}}{Shared latent coordinates for node_i in dimension h, corresponding to edge probability.}
#'            \item{\code{"x_ithk"}}{Layer-specific latent coordinates for node_i at layer_k in dimension h, corresponding to edge probability.}
#'            \item{\code{"tau_h_shared"}}{MCMC chain for the shrinkage parameter of the global latent space in dimension h, corresponding to edge probability.}
#'            \item{\code{"tau_h_k"}}{MCMC chain for the shrinkage parameter of the layer-specific latent space in dimension h, corresponding to edge probability.}
#'            \item{\code{"r_ijtk"}}{Expected value of edge weight between node_i and node_j at time t in layer_k.}
#'            \item{\code{"sigma_w_k"}}{Variance of edge weight at layer_k.}
#'            \item{\code{"lambda_tk_mcmc"}}{Baseline process for edge weight at layer_k.}
#'            \item{\code{"u_ith_shared"}}{Shared latent coordinates for node_i in dimension h, corresponding to edge weight.}
#'            \item{\code{"x_ithk"}}{Layer-specific latent coordinates for node_i at layer_k in dimension h, corresponding to edge weight}
#'            \item{\code{"rho_h_shared"}}{MCMC chain for the shrinkage parameter of the global latent space in dimension h, corresponding to edge weight}
#'            \item{\code{"rho_h_k"}}{MCMC chain for the shrinkage parameter of the layer-specific latent space in dimension h, corresponding to edge weight}
#'            }
#' @param node_i character/numeric. Id of a node in the network.
#' @param node_j character/numeric. Id of a node in the network.
#' @param layer_k character/numeric. Id of a layer in the network.
#' @param h numeric. Dimension of the latent space.
#' @param pred_p character/numeric. Id of a predictor in the exogenous data.
#' @param lat_space character In a directed network, specifies if the latent coordinate is in the "send" or "receive" space.
#' @param cred_int_type character. Type of credible intervals: Empirical or normal approximation
#' @param n_burn Integer. Number of iterations discarded as part of the MCMC warming up period at the beginning of the chain.
#' @param n_thin Integer. Number of iterations discarded for thining the chain (reducing the autocorrelation). We keep 1 of every n_thin iterations.
#' 
#' @import coda
#' @import ggplot2
#' @importFrom stats qnorm
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' plot_dmn_mcmc( x=dmn_mcmc,
#'                param="pi_ijtk",
#'                node_i=dmn_mcmc$node_all[2],
#'                node_j=dmn_mcmc$node_all[1],
#'                layer_k=1 )
#'                
#' plot_dmn_mcmc( x=dmn_mcmc,
#'                param="mu_tk",
#'                layer_k=dmn_mcmc$layer_all[1] )
#'                
#' plot_dmn_mcmc( x=dmn_mcmc,
#'                param="x_ith_shared",
#'                node_i=dmn_mcmc$node_all[1],
#'                h=1 )
#' 
#' plot_dmn_mcmc( x=dmn_mcmc,
#'                param="x_ithk",
#'                node_i=dmn_mcmc$node_all[1],
#'                h=1, layer_k=1 )
#' }
#' @export
#' 

plot_dmn_mcmc <- function( x,
                           param = c( "pi_ijtk",
                                      "mu_tk", "x_ith_shared", "x_ithk",
                                      "tau_h_shared","tau_h_k",
                                      "r_ijtk","sigma_w_k",
                                      "lambda_tk","u_ith_shared","u_ithk",
                                      "rho_h_shared","rho_h_k" )[1],
                           node_i=NULL, node_j=NULL, layer_k=NULL, h=NULL, pred_p=NULL,
                           lat_space=NULL,
                           cred_int_type = c("empirical","norm")[1],
                           n_burn=0, n_thin=1 ) {
  
  if( !is.element( param , c( "pi_ijtk",
                              "mu_tk", "x_ith_shared", "x_ithk",
                              "tau_h_shared","tau_h_k",
                              "r_ijtk","sigma_w_k",
                              "lambda_tk","u_ith_shared","u_ithk",
                              "rho_h_shared","rho_h_k" ) ) ) {stop("param=",param," not supported.")}
  
  directed <- x$directed
  weighted <- x$weighted
  V_net <- dim(x$y_ijtk)[1]
  T_net <- dim(x$y_ijtk)[3]
  K_net <- dim(x$y_ijtk)[4]
  
  H_dim <- x$H_dim
  R_dim <- x$R_dim
  
  if( !weighted & is.element( param , c( "r_ijtk_mcmc","sigma_w_k",
                                         "lambda_tk_mcmc","u_ith_shared_mcmc","u_ithk_mcmc",
                                         "rho_h_shared_mcmc","rho_h_k_mcmc" ) ) ) {stop("param=",param," not applicable for unweighted networks.")}
  n_iter_mcmc <- dim(x$mu_tk_mcmc)[3]
  
  # cat("\n",
  #     "V_net=",V_net,", T_net=",T_net,", K_net=",K_net,"\n",
  #     "H_dim=",H_dim,", R_dim=",R_dim,"\n",
  #     "n_iter_mcmc=",n_iter_mcmc,
  #     "\n\n", sep="")
  
  iter_out_mcmc <- seq(from=n_burn+1,to=n_iter_mcmc,by=n_thin)
  
  if( is.element(param,"pi_ijtk") ){
    if(is.null(node_i)){ node_i=x$node_all[2]; warning("Plotting node_i=",node_i," as node_i was not specified") }
    if(!is.element(node_i,x$node_all)) {stop("node_i=",node_i," is not a valid node")}
    i <- match(node_i,x$node_all)
    
    if(is.null(node_j)){ node_j=x$node_all[1]; warning("Plotting node_j=",node_j," as node_j was not specified") }
    if(!is.element(node_j,x$node_all)) {stop("node_i=",node_j," is not a valid node")}
    j <- match(node_j,x$node_all)
    
    if(!directed){
      ij_aux <- c(i,j)
      i <- max(ij_aux); node_i <- x$node_all[i]
      j <- min(ij_aux); node_j <- x$node_all[j]
    }
    
    if(is.null(layer_k)){ layer_k=x$layer_all[1]; warning("Plotting layer_k=",layer_k," as layer_k was not specified") }
    if(!is.element(layer_k,x$layer_all)) {stop("layer_k=",layer_k," is not a valid layer")}
    k <- match(layer_k,x$layer_all)
    
    param_mcmc_chain <- coda::mcmc( t(x$pi_ijtk_mcmc[i,j,,k,iter_out_mcmc]) )
    
  } else if( is.element(param,"mu_tk") ){
    
    if(is.null(layer_k)){ layer_k=x$layer_all[1]; warning("Plotting layer_k=",layer_k," as layer_k was not specified") }
    if(!is.element(layer_k,x$layer_all)) {stop("layer_k=",layer_k," is not a valid layer")}
    k <- match(layer_k,x$layer_all)
    
    param_mcmc_chain <- coda::mcmc( t(x$mu_tk_mcmc[,k,iter_out_mcmc]) )
    
  } else if( is.element(param,"x_ith_shared") ) {
    
    if(is.null(node_i)){ node_i=x$node_all[2]; warning("Plotting node_i=",node_i," as node_i was not specified") }
    if(!is.element(node_i,x$node_all)) {stop("node_i=",node_i," is not a valid node")}
    i <- match(node_i,x$node_all)
    
    if(is.null(h)){ h=1; warning("Plotting h=",h," as h was not specified") }
    if((h<=0)|(h>H_dim)){ stop("h must be an integer between 0 and H_dim=",H_dim) }
    
    if(directed){
      if(is.null(lat_space)){ lat_space="send"; warning("Plotting lat_space=",lat_space," as lat_space was not specified") }
      if(!is.element(lat_space,c("send","receive"))) {stop("lat_space=",lat_space,' must be "send" or "receive".')}
      dir <- match(lat_space,c("send","receive"))
      param_mcmc_chain <- coda::mcmc( t(x$x_ith_shared_mcmc[[dir]][i,,h,iter_out_mcmc]) )
    } else {
      param_mcmc_chain <- coda::mcmc( t(x$x_ith_shared_mcmc[i,,h,iter_out_mcmc]) )
    }
    
  } else if( is.element(param,"x_ithk") ) {
    
    if(is.null(node_i)){ node_i=x$node_all[2]; warning("Plotting node_i=",node_i," as node_i was not specified") }
    if(!is.element(node_i,x$node_all)) {stop("node_i=",node_i," is not a valid node")}
    i <- match(node_i,x$node_all)
    
    if(is.null(h)){ h=1; warning("Plotting h=",h," as h was not specified") }
    if((h<=0)|(h>H_dim)){ stop("h must be an integer between 0 and H_dim=",H_dim) }
    
    if(is.null(layer_k)){ layer_k=x$layer_all[1]; warning("Plotting layer_k=",layer_k," as layer_k was not specified") }
    if(!is.element(layer_k,x$layer_all)) {stop("layer_k=",layer_k," is not a valid layer")}
    k <- match(layer_k,x$layer_all)
    
    if(directed){
      if(is.null(lat_space)){ lat_space="send"; warning("Plotting lat_space=",lat_space," as lat_space was not specified") }
      if(!is.element(lat_space,c("send","receive"))) {stop("lat_space=",lat_space,' must be "send" or "receive".')}
      dir <- match(lat_space,c("send","receive"))
      param_mcmc_chain <- coda::mcmc( t(x$x_ithk_mcmc[[dir]][i,,h,k,iter_out_mcmc]) )
    } else {
      param_mcmc_chain <- coda::mcmc( t(x$x_ithk_mcmc[i,,h,k,iter_out_mcmc]) )
    }
  } else if( is.element(param,"r_ijtk") ){
    
    if(is.null(node_i)){ node_i=x$node_all[2]; warning("Plotting node_i=",node_i," as node_i was not specified") }
    if(!is.element(node_i,x$node_all)) {stop("node_i=",node_i," is not a valid node")}
    i <- match(node_i,x$node_all)
    
    if(is.null(node_j)){ node_j=x$node_all[1]; warning("Plotting node_j=",node_j," as node_j was not specified") }
    if(!is.element(node_j,x$node_all)) {stop("node_i=",node_j," is not a valid node")}
    j <- match(node_j,x$node_all)
    
    if(is.null(layer_k)){ layer_k=x$layer_all[1]; warning("Plotting layer_k=",layer_k," as layer_k was not specified") }
    if(!is.element(layer_k,x$layer_all)) {stop("layer_k=",layer_k," is not a valid layer")}
    k <- match(layer_k,x$layer_all)
    
    param_mcmc_chain <- coda::mcmc( t(x$r_ijtk_mcmc[i,j,,k,iter_out_mcmc]) )
    
  } else if( is.element(param,"sigma_w_k") ){
    
    param_mcmc_chain <- coda::mcmc( t(x$sigma_w_k_mcmc[,iter_out_mcmc]) )
    colnames(param_mcmc_chain) <- paste("layer_",1:K_net,sep="")
    p <- bayesplot::mcmc_trace(param_mcmc_chain,pars=colnames(param_mcmc_chain))
    p <- p + labs(title="sigma_w_k",subtitle="MCMC trace")
    return(p)
    
  } else if( is.element(param,"lambda_tk") ){
    
    if(is.null(layer_k)){ layer_k=x$layer_all[1]; warning("Plotting layer_k=",layer_k," as layer_k was not specified") }
    if(!is.element(layer_k,x$layer_all)) {stop("layer_k=",layer_k," is not a valid layer")}
    k <- match(layer_k,x$layer_all)
    
    param_mcmc_chain <- coda::mcmc( t(x$lambda_tk_mcmc[,k,iter_out_mcmc]) )
    
  }
  
  cred_int_probs = c(0.05,0.25)
  cred_int_quantiles = sort(c(cred_int_probs/2,1-cred_int_probs/2))
  
  summary_mcmc <- summary(param_mcmc_chain,quantiles=cred_int_quantiles)
  colnames(summary_mcmc[[2]]) <- paste("qmcmc_",100*cred_int_quantiles,sep="")
  summary_mcmc <- data.frame( cbind( time=x$time_all,
                                     summary_mcmc[[1]],
                                     summary_mcmc[[2]] ))
  
  summary_mcmc_aux <- sapply( cred_int_quantiles, function(x) { qnorm( p=x,
                                                                       mean=summary_mcmc$Mean,
                                                                       sd=summary_mcmc$SD ) },
                              simplify=TRUE )
  colnames(summary_mcmc_aux) <- paste("qnorm_",100*cred_int_quantiles,sep="")
  summary_mcmc <- cbind(summary_mcmc,summary_mcmc_aux)
  rownames(summary_mcmc) <- NULL
  
  if( is.element(cred_int_type,"empirical") ){
    p <- ggplot() +
      # geom_ribbon( aes( ymin=Mean-SD,
      #                   ymax=Mean+SD,
      #                   x=time ),
      #              fill="grey80", data=summary_mcmc ) +
      geom_ribbon( aes( ymin=get(paste("qmcmc_",100*cred_int_probs[1]/2,sep="")),
                        ymax=get(paste("qmcmc_",100*(1-cred_int_probs[1]/2),sep="")),
                        x=time ),
                   fill="grey50", data=summary_mcmc, alpha=0.25 ) +
      geom_ribbon( aes( ymin=get(paste("qmcmc_",100*cred_int_probs[2]/2,sep="")),
                        ymax=get(paste("qmcmc_",100*(1-cred_int_probs[2]/2),sep="")),
                        x=time ),
                   fill="grey50", data=summary_mcmc, alpha=0.50 ) +
      geom_line( aes(y=Mean,x=time),col="red", data=summary_mcmc )
    
  } else if( is.element(cred_int_type,"norm") ){
    p <- ggplot() +
      geom_ribbon( aes( ymin=get(paste("qnorm_",100*cred_int_probs[1]/2,sep="")),
                        ymax=get(paste("qnorm_",100*(1-cred_int_probs[1]/2),sep="")),
                        x=time ),
                   fill="grey50", data=summary_mcmc, alpha=0.25 ) +
      geom_ribbon( aes( ymin=get(paste("qnorm_",100*cred_int_probs[2]/2,sep="")),
                        ymax=get(paste("qnorm_",100*(1-cred_int_probs[2]/2),sep="")),
                        x=time ),
                   fill="grey50", data=summary_mcmc, alpha=0.50 ) +
      geom_line( aes(y=Mean,x=time),col="red", data=summary_mcmc )
  }
  
  if( is.element(param,"pi_ijtk") ){
    p <- p + geom_point( aes( y=as.numeric(x$y_ijtk[i,j,,k]>0),
                              x=x$time_all[x$time_all_idx_net] ) ) +
      labs(x="time",y="pi_ijtk",title="pi_ijtk",subtitle=paste(node_i,"->",node_j,", layer_k=",layer_k,sep=""))+
      coord_cartesian(ylim=c(0,1))
  } else if( is.element(param,"mu_tk") ){
    p <- p + labs(x="time",y="mu_tk",title="mu_tk",subtitle=paste("layer_k=",layer_k,sep=""))
  } else if( is.element(param,"x_ith_shared") ) {
    p <- p + labs(x="time",y="x_ith_shared",title="x_ith_shared",subtitle=paste("node_i=",node_i,", h=",h,sep=""))
  } else if( is.element(param,"x_ithk") ) {
    p <- p + labs(x="time",y="x_ithk",title="x_ithk",subtitle=paste("node_i=",node_i,", h=",h,", layer_k=",layer_k,sep=""))
  } else if( is.element(param,"r_ijtk") ){
    p <- p + geom_line( aes( y=x$y_ijtk[i,j,,k],
                              x=x$time_all[x$time_all_idx_net]), col="blue" ) +
      labs(x="time",y="r_ijtk",title="Expected weight",subtitle=paste(node_i,"->",node_j,", layer_k=",layer_k,sep=""))
  } else if( is.element(param,"lambda_tk") ){
    p <- p + labs(x="time",y="lambda_tk",title="lambda_tk",subtitle=paste("layer_k=",layer_k,sep=""))
  }
  
  return(p)
  
}
