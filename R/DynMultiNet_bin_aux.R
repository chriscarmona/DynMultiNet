
#' @export
get_y_ijtk_from_edges <- function( net_data,
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
  
  y_ijtk <- array( data=NA,
                   dim=c(V_net,V_net,T_net,K_net) )
  
  for( k in 1:K_net) {
    for( t in 1:T_net) {
      #t<-1;k<-1
      y_ijtk[,,t,k][lower.tri(y_ijtk[,,t,k])] <- 0
    }
  }
  for( row_i in 1:nrow(net_data) ){
    # row_i <- 1
    aux_ij <- match(net_data[row_i,c("source","target")],node_all)
    i <- max(aux_ij)
    j <- min(aux_ij)
    t <- match(net_data[row_i,"time"],time_all)
    k <- match(net_data[row_i,"layer"],layer_all)
    if(net_data[row_i,"weight"]>0){
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
  pred_id_global<-NULL; pred_id_layer<-NULL; pred_id_edge<-NULL
  z_tp<-NULL; z_tkp<-NULL; z_ijtkp<-NULL
  beta_z_global<-NULL; beta_z_layer<-NULL; beta_z_edge<-NULL
  
  if(!is.null(pred_data)) {
    #colnames_pred_data <- c("source","target","time","layer","id","pred")
    
    ## Predictors ##
    pred_all <- sort(unique(unlist(pred_data$id)))
    P_pred <- length(pred_all)
    
    cat("Procesing predictors data...\n")
    ## Global ##
    cond_pred <- apply( matrix(!is.na(pred_data),nrow=nrow(pred_data),ncol=6), 1, identical, c(F,F,T,F,T,T) )
    if( any(cond_pred) ) {
      pred_id_global <- pred_data[cond_pred,c("id","layer")]
      pred_id_global <- pred_id_global[!duplicated(pred_id_global),]
      z_tp <- matrix(NA,nrow=T_net,ncol=P_pred)
      for(row_i in which(cond_pred)){ # row_i <- which(cond_pred)[1]
        t <- match(pred_data[row_i,"time"],time_all)
        p <- match(pred_data[row_i,"id"],pred_all)
        z_tp[t,p] <- pred_data[cond_pred,][row_i,"pred"]
      }
      rm(row_i,t,p)
    }
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
                    pred_id_global=pred_id_global, pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                    z_tp=z_tp, z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                    beta_z_global=beta_z_global, beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge )
  
  return(pred_net)
  
}
