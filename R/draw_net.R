#' @title
#'     Drawing Inferred Network
#'
#' @description
#'     Drawing Inferred Network
#' 
#' @param x an object of class "\code{DynMultiNet}"
#' @param character, One of c("pi_ijtk","mu_ijtk", "E_y_ijtk")
#' 
#' @import igraph
#' 
#' @export
draw_net <- function( x,
                      param = c( "pi_ijtk",
                                 "mu_ijtk",
                                 "E_y_ijtk" )[1],
                      layer_k=NULL,
                      time=NULL,
                      FUN=median,
                      max_edge_width=1,
                      v_size=7,
                      no_id=FALSE,
                      ... ) {
  
  if(!is.element(param,c("pi_ijtk","mu_ijtk","E_y_ijtk")) | length(param)>1 ){ stop('"param" not supported') }
  if(is.null(time)){ time=1; warning('Plotting time=1 as "time" was not specified') }
  if(is.null(layer_k)){ layer_k=x$layer_all[1]; warning("Plotting layer_k=",layer_k," as layer_k was not specified") }
  k <- match(layer_k,x$layer_all)
  
  if( is.element(param,c("pi_ijtk","mu_ijtk")) ){
    net_data_obs <- reshape2::melt( apply(x[[paste(param,"_mcmc",sep="")]][,,time,k,],1:2,FUN) )
  } else if( is.element(param,c("E_y_ijtk")) ){
    net_data_obs <- reshape2::melt( apply(x[["pi_ijtk_mcmc"]][,,time,k,]*x[["mu_ijtk_mcmc"]][,,time,k,],1:2,FUN) )
  }
  colnames(net_data_obs) <- c("source","target","weight")
  net_data_obs <- net_data_obs[!is.na(net_data_obs[,3]),]
  
  net <- igraph::graph_from_data_frame( net_data_obs , directed=x$directed )
  if(no_id){ V(net)$name <- rank(-strength(net,mode="all")) }
  
  # Node strength
  # node_strength_i <- data.frame( name=V(net)$name,
  #                                in_degree_prob=strength(net,mode="in"),
  #                                out_degree_prob=strength(net,mode="out"),
  #                                sum_degree_prob=strength(net,mode="all") )
  # rownames(node_strength_i)<-NULL
  # node_strength_i[order(node_strength_i$sum_degree_prob,decreasing=T),]
  
  # vertex sizes
  vertex_size <- setNames( strength(net,mode="all"),
                           V(net)$name );
  vertex_size <- v_size*2*plogis(scale(vertex_size)[,1])
  
  igraph::plot.igraph( net,
                       # vertex.label.color="black",
                       # vertex.label.cex=1,
                       vertex.size=vertex_size,
                       # layout=layout_nicely,
                       layout=layout_with_fr,
                       # layout=layout_with_drl,
                       # layout=layout_with_graphopt,
                       edge.arrow.size=0,
                       edge.width=E(net)$weight/max(1,E(net)$weight)*max_edge_width,
                       ... )
}
