
.onUnload <- function (libpath) {
  library.dynam.unload("DynMultiNet", libpath)
}

#' @keywords internal
nato0 <- function(x){x[is.na(x)]<-0;x}

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



#' @export
get_y_ijtk_from_edges <- function( net_data,
                                   directed=FALSE,
                                   weighted=FALSE,
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



#' @keywords internal
get_x_ijtkp_mat <- function( x_ijtkp,
                             directed=FALSE,
                             weighted=FALSE ) {
  # Transforms the 5-dim x_ijtkp into a matrix x_ijtkp_mat
  # Let beta=(beta_1(t_1),...,beta_P(t_1),...,beta_1(t_T),,beta_P(t_T))
  # then
  # Y = x_ijtkp_mat %*% beta
  
  V_net <- dim(x_ijtkp)[1]
  T_net <- dim(x_ijtkp)[3]
  K_net <- dim(x_ijtkp)[4]
  P_pred <- dim(x_ijtkp)[5]
  
  if(!directed) {
    valid_idx <- which(c(lower.tri(diag(V_net))))
  } else {
    valid_idx <- which(c(diag(V_net))==0)
  }
  x_ijtkp_mat <- matrix(0,length(valid_idx)*T_net*K_net,P_pred*T_net)
  for(k in 1:K_net) { #k<-1
    for(t in 1:T_net) { #t<-1
      idx_row <- (k-1)*T_net*length(valid_idx)+(t-1)*length(valid_idx)+1:length(valid_idx)
      idx_col <- t+(1:P_pred-1)*T_net
      x_ijtkp_mat[idx_row,idx_col] <- matrix( c(x_ijtkp[,,t,k,]), nrow=V_net*V_net, ncol=P_pred )[valid_idx,]
    }
  }
  
  return(x_ijtkp_mat)
}



#' @keywords internal
get_GHW_t <- function( sigma_U,
                       sigma_A,
                       time_all,
                       nGP_mat_approx=FALSE ) {
  # From Zhu2013
  # state_beta(t) = c(U(t),U'(t),A(t))
  # state_beta(t+1) = G(t)*state_beta(t)+w(t)
  if( !all.equal(time_all,sort(time_all,decreasing=F)) ){stop("time_all is not ordered increasingly")}
  T_net <- length(time_all)
  diff_time_all <- c(diff(time_all),1)
  
  G <- array(diag(3),dim=c(3,3,T_net))
  G[1,2,] <- G[2,3,] <- diff_time_all
  if(!nGP_mat_approx){
    G[1,3,] <- (1.0/2)*diff_time_all^2
  }
  
  H <- array(diag(3),dim=c(3,3,T_net))
  if(nGP_mat_approx){
    H <- H[,2:3,]
  }
  
  if(!nGP_mat_approx){
    W <- Wchol <- array(diag(3),dim=c(3,3,T_net))
    # diag
    W[1,1,] <- (1.0/3)*sigma_U^2*diff_time_all^3 + (1.0/20)*sigma_A^2*diff_time_all^5
    W[2,2,] <- sigma_U^2*diff_time_all + (1.0/3)*sigma_A^2*diff_time_all^3
    W[3,3,] <- sigma_A^2*diff_time_all
    # out-diag
    W[2,1,] <- W[1,2,] <- (1.0/2)*sigma_U^2*diff_time_all^2 + (1.0/8)*sigma_A^2*diff_time_all^4
    W[3,1,] <- W[1,3,] <- (1.0/6)*sigma_A^2*diff_time_all^3
    W[3,2,] <- W[2,3,] <- (1.0/2)*sigma_A^2*diff_time_all^2
  } else {
    W <- Wchol <- array(diag(2),dim=c(2,2,T_net))
    W[1,1,] <- sigma_U^2*diff_time_all
    W[2,2,] <- sigma_A^2*diff_time_all
  }
  for(t in 1:T_net){
    Wchol[,,t] = chol( W[,,t] )
  }
  
  
  return( list( G=G,H=H,W=W,Wchol=Wchol,
                nGP_mat_approx=nGP_mat_approx ) )
}



#' @keywords internal
get_nGP_mat_net <- function( nGP_sigma_net,
                             nGP_mat_approx=FALSE,
                             directed=FALSE ) {
  
  V_net <- nGP_sigma_net$V_net
  T_net <- nGP_sigma_net$T_net
  K_net <- nGP_sigma_net$K_net
  
  time_all <- nGP_sigma_net$time_all
  # Output #
  nGP_mat_net <- list( baseline_k = NULL,
                       coord_i = NULL,
                       coord_ik = NULL,
                       add_eff_i = NULL,
                       add_eff_ik = NULL,
                       time_all=time_all,
                       nGP_mat_approx=nGP_mat_approx,
                       directed=directed )
  if(!nGP_mat_approx){
    # exact
    nGP_mat_net$baseline_k = list( G=array(NA,dim=c(K_net,3,3,T_net)),
                                   H=array(NA,dim=c(K_net,3,3,T_net)),
                                   Wchol=array(NA,dim=c(K_net,3,3,T_net)) )
    nGP_mat_net$coord_i = list( G=array(NA,dim=c(V_net,3,3,T_net)),
                                H=array(NA,dim=c(V_net,3,3,T_net)),
                                Wchol=array(NA,dim=c(V_net,3,3,T_net)) )
    nGP_mat_net$coord_ik = list( G=array(NA,dim=c(V_net,K_net,3,3,T_net)),
                                 H=array(NA,dim=c(V_net,K_net,3,3,T_net)),
                                 Wchol=array(NA,dim=c(V_net,K_net,3,3,T_net)) )
  } else {
    # approx
    nGP_mat_net$baseline_k = list( G=array(NA,dim=c(K_net,3,3,T_net)),
                                   H=array(NA,dim=c(K_net,3,2,T_net)),
                                   Wchol=array(NA,dim=c(K_net,2,2,T_net)) )
    nGP_mat_net$coord_i = list( G=array(NA,dim=c(V_net,3,3,T_net)),
                                H=array(NA,dim=c(V_net,3,2,T_net)),
                                Wchol=array(NA,dim=c(V_net,2,2,T_net)) )
    nGP_mat_net$coord_ik = list( G=array(NA,dim=c(V_net,K_net,3,3,T_net)),
                                 H=array(NA,dim=c(V_net,K_net,3,2,T_net)),
                                 Wchol=array(NA,dim=c(V_net,K_net,2,2,T_net)) )
  }
  
  # Computation #
  for(k in 1:K_net) {
    # variance of baseline process #
    aux <- get_GHW_t( sigma_U=nGP_sigma_net$baseline_k[k,1],
                      sigma_A=nGP_sigma_net$baseline_k[k,2],
                      time_all=time_all,
                      nGP_mat_approx=nGP_mat_approx )
    nGP_mat_net$baseline_k$G[k,,,] <- aux$G
    nGP_mat_net$baseline_k$H[k,,,] <- aux$H
    nGP_mat_net$baseline_k$Wchol[k,,,] <- aux$Wchol
  }
  
  # variance of latent coordinates #
  for(i in 1:V_net) {
    # variance of global coordinates #
    aux <- get_GHW_t( sigma_U=nGP_sigma_net$coord_i[k,1],
                      sigma_A=nGP_sigma_net$coord_i[k,2],
                      time_all=time_all,
                      nGP_mat_approx=nGP_mat_approx )
    nGP_mat_net$coord_i$G[i,,,] <- aux$G
    nGP_mat_net$coord_i$H[i,,,] <- aux$H
    nGP_mat_net$coord_i$Wchol[i,,,] <- aux$Wchol
    # variance of layer-specific coordinates #
    if(K_net>1){
      for(k in 1:K_net) {
        aux <- get_GHW_t( sigma_U=nGP_sigma_net$coord_ik[i,k,1],
                          sigma_A=nGP_sigma_net$coord_ik[i,k,2],
                          time_all=time_all,
                          nGP_mat_approx=nGP_mat_approx )
        nGP_mat_net$coord_ik$G[i,k,,,] <- aux$G
        nGP_mat_net$coord_ik$H[i,k,,,] <- aux$H
        nGP_mat_net$coord_ik$Wchol[i,k,,,] <- aux$Wchol
      }
    }
  }
  # variance of additive effects #
  if(!is.null(nGP_mat_net$add_eff_i)) {
    if(!nGP_mat_approx){
      # exact
      nGP_mat_net$add_eff_i = list( G=array(NA,dim=c(V_net,3,3,T_net)),
                                    H=array(NA,dim=c(V_net,3,3,T_net)),
                                    Wchol=array(NA,dim=c(V_net,3,3,T_net)) )
    } else {
      # approx
      nGP_mat_net$add_eff_i = list( G=array(NA,dim=c(V_net,3,3,T_net)),
                                    H=array(NA,dim=c(V_net,3,2,T_net)),
                                    Wchol=array(NA,dim=c(V_net,2,2,T_net)) )
    }
    
    for(i in 1:V_net) {
      aux <- get_GHW_t( sigma_U=nGP_sigma_net$add_eff_i[i,1],
                        sigma_A=nGP_sigma_net$add_eff_i[i,2],
                        time_all=time_all,
                        nGP_mat_approx=nGP_mat_approx )
      nGP_mat_net$add_eff_i$G[i,,,] <- aux$G
      nGP_mat_net$add_eff_i$H[i,,,] <- aux$H
      nGP_mat_net$add_eff_i$Wchol[i,,,] <- aux$Wchol
    }
  }
  # variance of layer-specific additive effects #
  if(K_net>1){
    if(!is.null(nGP_mat_net$add_eff_ik)) {
      if(!nGP_mat_approx){
        # exact
        nGP_mat_net$add_eff_ik = list( G=array(NA,dim=c(V_net,K_net,3,3,T_net)),
                                       H=array(NA,dim=c(V_net,K_net,3,3,T_net)),
                                       Wchol=array(NA,dim=c(V_net,K_net,3,3,T_net)) )
      } else {
        # approx
        nGP_mat_net$add_eff_ik = list( G=array(NA,dim=c(V_net,K_net,3,3,T_net)),
                                       H=array(NA,dim=c(V_net,K_net,3,2,T_net)),
                                       Wchol=array(NA,dim=c(V_net,K_net,2,2,T_net)) )
      }
      
      for(k in 1:K_net) {
        aux <- get_GHW_t( sigma_U=nGP_sigma_net$add_eff_ik[i,k,1],
                          sigma_A=nGP_sigma_net$add_eff_ik[i,k,2],
                          time_all=time_all,
                          nGP_mat_approx=nGP_mat_approx )
        nGP_mat_net$add_eff_ik$G[i,k,,,] <- aux$G
        nGP_mat_net$add_eff_ik$H[i,k,,,] <- aux$H
        nGP_mat_net$add_eff_ik$Wchol[i,k,,,] <- aux$Wchol
      }
    }
  }
  
  return(nGP_mat_net)
}
