
#' @export
get_y_ijtk_from_edges <- function( edges_data,
                                   quiet=FALSE ) {
  
  if(!quiet){
    cat("Procesing Network data...\n")
  }
  
  node_all <- sort(unique(unlist(edges_data[,c("source","target")])))
  V_net <- length(node_all)
  
  time_all <- sort(unique(unlist(edges_data$time)))
  T_net <- length(time_all)
  
  layer_all <- sort(unique(unlist(edges_data$layer)))
  K_net <- length(layer_all)
  
  y_ijtk <- array( data=NA,
                   dim=c(V_net,V_net,T_net,K_net) )
  
  for( k in 1:K_net) {
    for( t in 1:T_net) {
      #t<-1;k<-1
      y_ijtk[,,t,k][lower.tri(y_ijtk[,,t,k])] <- 0
    }
  }
  for( row_i in 1:nrow(edges_data) ){
    # row_i <- 1
    aux_ij <- match(edges_data[row_i,c("source","target")],node_all)
    i <- max(aux_ij)
    j <- min(aux_ij)
    t <- match(edges_data[row_i,"time"],time_all)
    k <- match(edges_data[row_i,"layer"],layer_all)
    if(edges_data[row_i,"weight"]>0){
      y_ijtk[i,j,t,k] <- 1
    }
  }
  for( k in 1:K_net) {
    for( t in 1:T_net) {
      #t<-1;k<-1
      diag(y_ijtk[,,t,k]) <- NA
    }
  }
  if(!quiet){
    cat("done!\n")
  }
  return(y_ijtk)
  
}
