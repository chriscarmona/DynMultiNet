#' @title
#'    Bayesian Learning of Dynamic Multilayer Networks with binary data
#'
#' @description
#'    \code{DynMultiNet_bin} Implements model from Durante and Dunson, 2018
#'
#' @param net_data Data frame with network information
#' @param pred_data Data frame with linked predictors information
#' @param H_dim Latent space dimension
#' @param n_iter_mcmc number of iterations for the MCMC
#' @param out_file Indicates a file (.RData) where the output should be saved
#' @param log_file Indicates a file (.txt) where the log of the process should be saved
#' @param quiet_mcmc Silent mode, no progress update
#'
#'
#' @details
#'    The model assumes a latent variable approach
#'    
#'    \code{net_data} must be a data frame with the following columns
#'        \describe{
#'        \item{\code{source}}{Start Node}
#'        \item{\code{target}}{End node}
#'        \item{\code{weight}}{Edge weight}
#'        \item{\code{time}}{time associated with the edge}
#'        \item{\code{layer}}{Layer associated with the edge}
#'        }
#'
#' @return
#'    A list with the following components:
#' \describe{
#'     \item{\code{theta_mcmc}}{Matrix with the chain of the parameters in the model.}
#' }
#'
#'
#' @examples
#' 
#'    DynMultiNet_bin( net_data,
#'                     pred_data=NULL,
#'                     H_dim=10, k_x=0.10, k_mu=0.10, a_1=2, a_2=2.5,
#'                     n_iter_mcmc=100000,
#'                     out_file=NULL, log_file=NULL,
#'                     quiet_mcmc=FALSE )
#' 
#' @useDynLib DynMultiNet
#' 
#' @import BayesLogit
#' @importFrom abind abind
#' 
#' @export

DynMultiNet_bin <- function( net_data,
                             pred_data=NULL,
                             H_dim=10, k_x=0.10, k_mu=0.10, a_1=2, a_2=2.5,
                             n_iter_mcmc=100000,
                             out_file=NULL, log_file=NULL,
                             quiet_mcmc=FALSE ) {
  
  #### Start: Checking inputs ####
  colnames_net_data <- c("source","target","weight","time","layer")
  if( !all(is.element(colnames_net_data,colnames(net_data))) ) {
    stop('"net_data" must contain the following columns: ',paste(colnames_net_data,collapse=", ") )
  }
  net_data <- net_data[,colnames_net_data]
  #### End: Checking inputs ####
  
  
  #### Start: Global parameters ####
  node_all <- sort(unique(unlist(net_data[,c("source","target")])))
  V_net <- length(node_all)
  
  time_all <- sort(unique(unlist(net_data$time)))
  T_net <- length(time_all)
  
  layer_all <- sort(unique(unlist(net_data$layer)))
  K_net <- length(layer_all)
  #### End: Global parameters ####
  
  
  #### Start: MCMC initialization ####
  # Edge between actors i and j at time t in layer k
  cat("Procesing data...\n")
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
  cat("done!\n")
  
  # Augmented Polya-gamma data
  w_ijtk <- y_ijtk
  w_ijtk[!is.na(w_ijtk)] <- 0
  
  # Baseline parameter #
  # at time t for layer k
  mu_tk <- matrix( data=0,
                   nrow=T_net,
                   ncol=K_net )
  
  # Covariance matrix for baseline mu_t
  mu_tk_sigma <- array( data=0,
                        dim=c(T_net,T_net,K_net) )
  for(k in 1:K_net) { mu_tk_sigma[,,k]<-diag(T_net) }
  rm(k)
  
  # Covariance matrix prior for baseline mu_t
  mu_t_sigma_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_mu){ exp(-k*(x-y)^2) } )
  
  # Latent coordinates #
  # shared: hth coordinate of actor v at time t shared across the different layers
  x_iht <- array( data=0.01,
                  #data=runif(V_net*H_dim*T_net,-1,1),
                  dim=c(V_net,H_dim,T_net) )
  x_iht_mat <- aperm(a=x_iht,perm=c(1,3,2))
  dim(x_iht_mat) <- c(V_net,T_net*H_dim)
  if( !all(x_iht_mat[1,1:T_net]==x_iht[1,1,]) ){stop("there is a problem arranging x_iht into x_iht_mat")}
  x_iht_mat_mcmc <- array(NA,c(V_net,T_net*H_dim,n_iter_mcmc))
  
  if( K_net>1 ){
    # by layer: hth coordinate of actor v at time t specific to layer k
    x_ihtk <- array( data=runif(V_net*H_dim*T_net*K_net),
                     dim=c(V_net,H_dim,T_net,K_net) )
  }
  # Covariance matrix prior for coordinates x_t
  x_t_sigma_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_x){ exp(-k*(x-y)^2) } )
  
  # Linear predictor for the probability of an edge between actors i and j at time t in layer k
  s_ijtk <- array( data=NA,
                   dim=c(V_net,V_net,T_net,K_net) )
  for( k in 1:K_net ){
    for( t in 1:T_net ){
      #t<-1
      s_ijtk[,,t,k] <- mu_tk[t,k] + x_iht[,,t] %*% t(x_iht[,,t])
    }
  }
  s_ijtk[is.na(y_ijtk)]<-NA
  rm(k,t)
  
  # Probability of an edge between actors i and j at time t in layer k
  pi_ijtk <- plogis(s_ijtk)
  #pi_ijt[,,1]
  # all( abs(qlogis( pi_ijt ) - s_ijt)<1e-6,na.rm=T ) # TRUE
  
  # Shrinkage Parameters
  v_dim <- matrix( NA, nrow=H_dim, ncol=1 )
  v_dim[1,1] <- rgamma(n=1,shape=a_1,rate=1); v_dim[-1,1] <- rgamma(n=H_dim-1,shape=a_2,rate=1)
  tau_h <- matrix(cumprod(v_dim), nrow=H_dim, ncol=1 )
  # 1/tau_h
  #### End: MCMC initialization ####
  
  
  
  #### Start: MCMC Sampling ####
  if(!quiet_mcmc){ cat("Sampling MCMC ...\n") }
  for ( iter_i in 1:n_iter_mcmc) {
    
    ### Step 1. Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
    #browser()
    for( i in 2:V_net) {
      for( j in 1:(i-1)) {
        for( t in 1:T_net) {
          for( k in 1:K_net) {
            # i<-1;j<-1;t<-1;k<-1
            w_ijtk[i,j,t,k] <- rpg( num=1, h=1, z=s_ijtk[i,j,t,k] )
          }
        }
      }
    }
    rm(i,j,t,k)
    
    
    ### Step 2. Sample mu_t_k from its conditional N-variate Gaussian posterior ###
    aux_sum_w_tk <- apply(w_ijtk,c(3,4),sum,na.rm=T)
    for( k in 1:K_net) {
      #k <- 1
      aux_mat <- solve( diag(aux_sum_w_tk[,k]) + solve(mu_t_sigma_prior) )
      if(!isSymmetric(aux_mat)) {aux_mat[upper.tri(aux_mat)] <- t(aux_mat)[upper.tri(aux_mat)]}
      mu_tk_sigma[,,k] <- aux_mat
      
      aux_vec1 <- matrix(NA,T_net,1)
      for( t in 1:T_net) {
        aux_vec1[t,] <- sum( y_ijtk[,,t,k] - 0.5 - w_ijtk[,,t,k] * (x_iht[,,t] %*% t(x_iht[,,t])), na.rm=TRUE )
      }
      mu_old <- mu_tk[,k]
      
      
      mu_t_new <- mvtnorm::rmvnorm( n=1,
                                     mean=mu_tk_sigma[,,k] %*% aux_vec1,
                                     sigma=mu_tk_sigma[,,k] )
      mu_t_new <- matrix(mu_t_new)
      if(T){
        k<-1
        plot(y=mu_tk[,k], x=1:T_net, col="blue", type="l", ylim=range(mu_tk[,k],mu_t_new),main="mu_tk[,k]")
        lines(y=mu_t_new, x=1:T_net, col="red")
        legend("topright",col=c("blue","red"),legend=paste("iter ",c(iter_i-1,iter_i),sep=""),lty=1)
        #if(iter_i==6){browser()}
      }
      mu_tk[,k] <- mu_t_new
    }
    
    rm(aux_sum_w_tk,k,aux_mat,aux_vec1)
    #write.csv(mu_tk[,1],paste(out_dir,"mu_t_iter_",iter_i,".csv",sep=""),row.names=F)
    
    ### Step 3. For each unit, block-sample the set of time-varying latent coordinates x_iht ###
    for( k in 1:K_net) {
      for(i in 1:V_net) {
        #k<-1; i<-3
        y_i <- c(t(rbind( matrix(y_ijtk[i,1:i,,k],i,T_net)[-i,],
                          matrix(y_ijtk[i:V_net,i,,k],V_net-i+1,T_net)[-1,] ) ))
        y_i <- matrix(y_i)
        w_i <- c(t(rbind( matrix(w_ijtk[i,1:i,,k],i,T_net)[-i,],
                          matrix(w_ijtk[i:V_net,i,,k],V_net-i+1,T_net)[-1,] ) ))
        w_diag_i <- diag(w_i)
        
        x_tilde_i <- x_iht_mat[rep((1:V_net)[-i],each=T_net),]
        
        x_tilde_i_rm <- matrix(TRUE,T_net,T_net*H_dim)
        for(t in 1:T_net){
          x_tilde_i_rm[t,seq(t,T_net*H_dim,by=T_net)] <- F
        }
        x_tilde_i_rm <- do.call(rbind, replicate(V_net-1, x_tilde_i_rm, simplify=FALSE))
        
        x_tilde_i[x_tilde_i_rm] <- 0
        
        if(F){
          s_i <- c(t(rbind( matrix(s_ijtk[i,1:i,,k],i,T_net)[-i,],
                            matrix(s_ijtk[i:V_net,i,,k],V_net-i+1,T_net)[-1,] ) ))
          s_i <- matrix(s_i)
          x_i <- matrix(x_iht_mat[i,])
          all(s_i == x_tilde_i %*% x_i)
          rm(s_i,x_i)
        }
        x_i_cov <- solve( t(x_tilde_i) %*% w_diag_i %*% x_tilde_i + kronecker( diag(as.numeric(tau_h)), solve(x_t_sigma_prior) ) )
        if(!isSymmetric(x_i_cov)) {x_i_cov[upper.tri(x_i_cov)] <- t(x_i_cov)[upper.tri(x_i_cov)]}
        x_i_mean <- x_i_cov %*% ( t(x_tilde_i) %*% ( y_i - kronecker(matrix(1,V_net-1,1),matrix(1,T_net,1)*0.5) - w_diag_i %*% kronecker(matrix(1,V_net-1,1),mu_tk[,k]) ) )
        x_i <- mvtnorm::rmvnorm( n=1,
                                 mean=x_i_mean,
                                 sigma=x_i_cov )
        x_iht_mat[i,] <- x_i
      }
    }
    rm(k,i,y_i,w_i,w_diag_i,x_tilde_i,x_tilde_i_rm,x_i_cov,x_i_mean)
    # redefine x_iht with the new sampled values
    x_iht <- x_iht_mat
    dim(x_iht) <- c(V_net,T_net,H_dim)
    x_iht <- aperm(a=x_iht,perm=c(1,3,2))
    if( !all(x_iht_mat[1,1:T_net]==x_iht[1,1,]) ){stop("there is a problem arranging x_iht into x_iht_mat")}
    
    #write.csv(x_iht_mat,paste(out_dir,"x_iht_mat_iter_",iter_i,".csv",sep=""),row.names=F)
    
    # Linear predictor and link probabilities
    for( k in 1:K_net ){
      for( t in 1:T_net ){
        s_ijtk[,,t,k] <- mu_tk[t,k] + x_iht[,,t] %*% t(x_iht[,,t])
      }
    }
    s_ijtk[is.na(y_ijtk)]<-NA
    rm(k,t)
    pi_ijtk <- plogis(s_ijtk)
    
    
    
    ### Step 4. Sample the global shrinkage hyperparameters from conditional gamma distributions ###
    for(h in 2:H_dim) {
      tau_h <- matrix(cumprod(v_dim), nrow=H_dim, ncol=1 )
      tau_minush_l <- matrix(tau_h, nrow=H_dim, ncol=H_dim )
      tau_minush_l[upper.tri(tau_minush_l)] <- NA ; tau_minush_l[1,1] <- NA
      tau_minush_l <- tau_minush_l / matrix(v_dim, nrow=H_dim, ncol=H_dim, byrow=T)
      aux_1 <- apply(tau_minush_l,2,sum,na.rm=TRUE)
      aux_2 <- vector(mode="numeric",length=V_net)
      for( i in 1:V_net ){aux_2[i]<-matrix(x_iht[i,h,],nrow=1) %*% solve(x_t_sigma_prior) %*% matrix(x_iht[i,h,],ncol=1)}
      if(h==1){
        v_dim[h,1] <- rgamma( n=1,
                              shape = a_1+0.5*(V_net*T_net*H_dim),
                              rate = 1+0.5*aux_1[h]*sum(aux_2) )
      }
      if(h>1){
        v_dim[h,1] <- rgamma( n=1,
                              shape = a_1+0.5*(V_net*T_net*(H_dim-h+1)),
                              rate = 1+0.5*aux_1[h]*sum(aux_2) )
      }
    }
    tau_h <- matrix(cumprod(v_dim), nrow=H_dim, ncol=1 )
    rm(h,i,tau_minush_l,aux_1,aux_2)
    
    
    if(!quiet_mcmc){
      if( is.element(iter_i, floor(n_iter_mcmc*seq(0.05,1,0.05)) ) ) {
        cat(round(100*iter_i/n_iter_mcmc),"% ",sep="")
      }
    }
  }
  if(!quiet_mcmc){cat("\nMCMC finished!\n")}
  #### End: MCMC Sampling ####
  
  
  
  out <- list(NULL)
  return( out )
  
}
