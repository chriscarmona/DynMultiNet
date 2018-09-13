
#' @importFrom stats var
#' @keywords internal
R_hat.mcmc <- function(x,m) {
  # Computes the Potential Scale Reduction Coefficient
  # Gelman et al. (2014) sec. 11.4 page 285
  n <- floor(length(x)/m)
  chain_split <- matrix(NA,nrow=n,ncol=m)
  for(i in 1:m) {
    chain_split[,i] <- x[(1+n*(i-1)):(n*i)]
  }
  B <- n * stats::var(apply(chain_split,2,mean))
  W <- mean(apply(chain_split,2,stats::var))
  var_hat <- ((n-1)/n)*W + (1/n)*B
  R_hat <- sqrt(var_hat/W)
  R_hat
}



#' @keywords internal
bounce_limit <- function(x,a,b){
  while( (x<a) || (x>b) ) {
    if(x < a) {
      x <- a + (a-x)
    }
    if(x > b) {
      x <- b - (x-b)
    }
  }
  return(x)
}



#' @keywords internal
get_y_ijtk_from_edges <- function( net_data,
                                   directed=TRUE,
                                   weighted=TRUE,
                                   self_edges=FALSE,
                                   quiet=FALSE ) {
  
  #### Start: Checking inputs ####
  colnames_net_data <- c("source","target","time","layer","weight")
  if( !all(is.element(colnames_net_data,colnames(net_data))) ) {
    stop('"net_data" must contain the following columns: ',paste(colnames_net_data,collapse=", ") )
  }
  net_data <- net_data[,colnames_net_data]
  rm(colnames_net_data)
  #### End: Checking inputs ####
  
  if(!quiet){
    cat("Procesing Network data...\n")
  }
  
  node_all <- sort(unique(unlist(net_data[,c("source","target")])))
  V_net <- length(node_all)
  
  time_all <- sort(unique(unlist(net_data$time)))
  T_net <- length(time_all)
  
  layer_all <- sort(unique(unlist(net_data$layer)))
  K_net <- length(layer_all)
  
  y_ijtk <- array( data=0,
                   dim=c(V_net,V_net,T_net,K_net),
                   dimnames=list(node_all,node_all,time_all,layer_all) )
  
  # for( k in 1:K_net) {
  #   for( t in 1:T_net) {
  #     #t<-1;k<-1
  #     y_ijtk[,,t,k][lower.tri(y_ijtk[,,t,k])] <- 0
  #   }
  # }
  
  for( row_i in 1:nrow(net_data) ){
    # row_i <- 1
    k <- match(net_data[row_i,"layer"],layer_all)
    t <- match(net_data[row_i,"time"],time_all)
    ij <- match(net_data[row_i,c("source","target")],node_all)
    y_ijtk[ ij[1], ij[2], t, k ] <- net_data$weight[row_i]
  }; rm(row_i,k,t,ij)
  
  if(!self_edges){
    for( k in 1:K_net) {
      for( t in 1:T_net) { #t<-1;k<-1
        diag(y_ijtk[,,t,k]) <- NA
      }
    }; rm(k,t)
  }
  
  if(!directed){
    for( k in 1:K_net) {
      for( t in 1:T_net) { #t<-1;k<-1
        y_ijtk[,,t,k][lower.tri(y_ijtk[,,t,k])] <- y_ijtk[,,t,k][lower.tri(y_ijtk[,,t,k])] + t(y_ijtk[,,t,k])[lower.tri(y_ijtk[,,t,k])]
        y_ijtk[,,t,k][upper.tri(y_ijtk[,,t,k])] <- NA
      }
    }; rm(k,t)
  }
  
  if(!weighted){
    y_ijtk[y_ijtk!=0] <- 1
    # y_ijtk[y_ijtk>0] <- 1
    # y_ijtk[y_ijtk<0] <- -1
  }
  
  if(!quiet){
    cat("done!\n")
  }
  return(y_ijtk)
  
}


#' @keywords internal
get_z_pred <- function( pred_data,
                        node_all, time_all, layer_all,
                        quiet=FALSE ) {
  
  V_net <- length(node_all)
  T_net <- length(time_all)
  K_net <- length(layer_all)
  
  #### Start: Checking inputs ####
  if(!is.null(pred_data)){
    colnames_pred_data <- c("source","target","time","layer","id","pred")
    if( !all(is.element(colnames_pred_data,colnames(pred_data))) ) {
      stop('"pred_data" must contain the following columns: ',paste(colnames_pred_data,collapse=", ") )
    }
    pred_data <- pred_data[,colnames_pred_data]
    if( any(is.na(pred_data[,c("time","id","pred")])) ) {
      stop('"pred_data" does NOT allow NAs in columns: ',paste(c("time","id","pred"),collapse=", ") )
    }
    
    rm(colnames_pred_data)
  }
  #### End: Checking inputs ####
  
  pred_all <- NULL
  pred_id_layer<-NULL; pred_id_edge<-NULL
  z_tp<-NULL; z_tkp<-NULL; z_ijtkp<-NULL
  
  if(!is.null(pred_data)) {
    #colnames_pred_data <- c("source","target","time","layer","id","pred")
    
    ## Predictors ##
    pred_all <- sort(unique(unlist(pred_data$id)))
    P_pred <- length(pred_all)
    
    cat("Procesing predictors data...\n")
    
    ## Layer specific ##
    cond_pred <- apply( matrix(!is.na(pred_data),nrow=nrow(pred_data),ncol=6), 1, identical, c(F,F,T,T,T,T) )
    if( any(cond_pred) ) {
      pred_id_layer <- pred_data[cond_pred,c("id","layer")]
      pred_id_layer <- pred_id_layer[!duplicated(pred_id_layer),]
      z_tkp <- array( NA, dim=c(T_net,K_net,P_pred) )
      for(row_i in which(cond_pred)){ # row_i <- which(cond_pred)[1]
        t <- match(pred_data[row_i,"time"],time_all)
        k <- match(pred_data[row_i,"layer"],layer_all)
        p <- match(pred_data[row_i,"id"],pred_all)
        z_tkp[t,k,p] <- pred_data[row_i,"pred"]
      }
      rm(row_i,t,k,p)
    }
    ## Edge specific ##
    # The edge specific predictors MUST specify also the asscociated layer.
    # there will be one coefficient for each predictor-layer combination
    cond_pred <- apply( matrix(!is.na(pred_data),nrow=nrow(pred_data),ncol=6), 1, identical, c(T,T,T,T,T,T) )
    if( any(cond_pred) ) {
      pred_id_edge <- pred_data[cond_pred,c("id","layer")]
      pred_id_edge <- pred_id_edge[!duplicated(pred_id_edge),]
      z_ijtkp <- array( NA, dim=c(V_net,V_net,T_net,K_net,P_pred) )
      for(row_i in which(cond_pred)){ # row_i <- which(cond_pred)[1]
        i <- match(pred_data[row_i,c("source")],node_all)
        j <- match(pred_data[row_i,c("target")],node_all)
        t <- match(pred_data[row_i,"time"],time_all)
        k <- match(pred_data[row_i,"layer"],time_all)
        p <- match(pred_data[row_i,"id"],pred_all)
        z_ijtkp[i,j,t,k,p] <- pred_data[row_i,"pred"]
      }
      rm(row_i,i,j,t,k,p)
    }
    rm(cond_pred)
    
    cat("done!\n")
  }
  
  pred_net <- list( pred_all=pred_all,
                    pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                    z_tp=z_tp, z_tkp=z_tkp, z_ijtkp=z_ijtkp )
  
  return(pred_net)
  
}
