#' @title
#'    MCMC algorithm for Dynamic Multilayer directed unweighted graphs
#'
#' @description
#'    \code{mcmc_d_1_w_1} Implements a Gibbs sampler MCMC algorithm for Dynamic Multilayer directed unweighted graphs.
#'
#' @param y_ijtk Array. Network data, with entry \code{y_ijtk[i,j,t,k]} representing the link from node i to node j at time t in layer k. 0 indicates no link.
#' @param node_all Character vector. Id's of nodes in the network.
#' @param time_all Numeric vector. Timestamps of all relevant epochs for the MCMC, those observed and those for forecast
#' @param layer_all Character vector. Id's of layers in the network.
#' @param H_dim Integer. Latent space dimension.
#' @param R_dim Integer. Latent space dimension, for layer specific latent vectors.
#' @param delta Positive scalar. Hyperparameter controlling for the smoothness in the dynamic of latent coordinates. Larger=smoother.
#' @param shrink_lat_space Boolean. Indicates if the space should be shrinked probabilistically.
#' @param a_1 Positive scalar. Hyperparameter controlling for number of effective dimensions in the latent space.
#' @param a_2 Positive scalar. Hyperparameter controlling for number of effective dimensions in the latent space.
#' @param procrustes_lat Boolean. Indicates if the latent coordinates should be stabilised using the procustres transformation.
#' @param n_iter_mcmc Integer. Number of iterations for the MCMC.
#' @param n_burn Integer. Number of iterations discarded as part of the MCMC warming up period at the beginning of the chain.
#' @param n_thin Integer. Number of iterations discarded for thining the chain (reducing the autocorrelation). We keep 1 of every n_thin iterations.
#' @param rds_file String. Indicates a file (.rds) where the output will be saved.
#' @param log_file String. Indicates a file (.txt) where the log of the process will be saved.
#' @param quiet_mcmc Boolean. Indicates if silent mode is preferes, if \code{FALSE} progress update is displayed.
#' @param parallel_mcmc Boolean. Indicates if the mcmc should be processed in parallel.
#' 
#' @import foreach
#' @import BayesLogit
#' @importFrom MCMCpack procrustes
#' 
#' @keywords internal
#' 


mcmc_d_1_w_1 <- function( y_ijtk,
                          node_all, time_all, layer_all,
                          
                          H_dim=10, R_dim=10,
                          
                          delta=36,
                          
                          shrink_lat_space=TRUE,
                          a_1=2, a_2=2.5,
                          
                          procrustes_lat=TRUE,
                          
                          n_iter_mcmc=10000, n_burn=floor(n_iter_mcmc/3), n_thin=3,
                          
                          rds_file=NULL, log_file=NULL,
                          quiet_mcmc=FALSE,
                          parallel_mcmc=FALSE ) {
  
  # This software deal with continuos (weighted) edges
  # y_ijtk[y_ijtk<0] <- 0 # Shall we restrict to being positive?
  
  V_net <- length(node_all)
  T_net <- length(time_all)
  K_net <- length(layer_all)
  
  # assume no self-edges
  diag_y_idx <- matrix(FALSE,V_net,V_net); diag(diag_y_idx)<-TRUE
  diag_y_idx <- array(diag_y_idx,dim=dim(y_ijtk))
  y_ijtk[diag_y_idx] <- 0
  
  ### iterations that will be reported ###
  # after burn-in period and thinning
  iter_out_mcmc <- seq(from=n_burn+1,to=n_iter_mcmc,by=n_thin)
  n_iter_mcmc_out <- length(iter_out_mcmc)
  
  
  
  #### Start: MCMC initialization ####
  # Edge between actors i and j at time t in layer k
  
  
  ### GAUSSIAN MODEL FOR WEIGHTS ###
  
  # y_ijtk ~ norm( mu_ijtk[t,k] , sigma_k[k]^2  )
  # mu_ijtk[t,k] = theta_tk[t,k] + t(u_ith[i,t,]) * v_ith[j,t,] + t(u_ithk[i,t,,k]) * v_ithk[j,t,,k]
  # u_ith[i,,h] ~ Norm( 0 , C(t,t) )
  # theta_tk[,k] ~ Norm( 0 , C(t,t) )
  
  
  
  # Baseline parameter for weights #
  # at time t for layer k
  theta_tk <- matrix( data=runif(T_net*K_net),
                      nrow=T_net,
                      ncol=K_net )
  theta_tk_mcmc <- array( NA, dim=c(T_net,K_net,n_iter_mcmc_out) )
  
  # Latent coordinates #
  # shared: hth coordinate of actor v at time t shared across the different layers
  # One latent space for each direction #
  uv_ith_shared <- list( send=array( data=runif(V_net*T_net*H_dim,-1,1),
                                    dim=c(V_net,T_net,H_dim) ),
                        receive=array( data=runif(V_net*T_net*H_dim,-1,1),
                                       dim=c(V_net,T_net,H_dim) ) )
  uv_ith_shared_mcmc <- list( send=array(NA,c(V_net,T_net,H_dim,n_iter_mcmc_out)),
                             receive=array(NA,c(V_net,T_net,H_dim,n_iter_mcmc_out)) )
  
  # by layer: hth coordinate of actor v at time t specific to layer k
  if( K_net>1 ){
    # One latent space for each direction #
    uv_ithk <- list( send=array( data=runif(V_net*T_net*R_dim*K_net),
                                dim=c(V_net,T_net,R_dim,K_net) ),
                    receive=array( data=runif(V_net*T_net*R_dim*K_net),
                                   dim=c(V_net,T_net,R_dim,K_net) ) )
    uv_ithk_mcmc <- list( send=array(NA,c(V_net,T_net,R_dim,K_net,n_iter_mcmc_out)),
                         receive=array(NA,c(V_net,T_net,R_dim,K_net,n_iter_mcmc_out)) )
  } else {
    uv_ithk <- NULL
    uv_ithk_mcmc <- NULL
  }
  
  # Weight variance for layer k
  sigma_k <- apply( y_ijtk, MARGIN=4, FUN=sd, na.rm=T)
  
  # Covariance matrix prior for latent coordinates
  cov_gp_prior <- outer( time_all, time_all, FUN=function(x,y,k=delta){ exp(-((x-y)/delta)^2) } )
  diag(cov_gp_prior) <- diag(cov_gp_prior) + 1e-3 # numerical stability
  cov_gp_prior_inv <- solve(cov_gp_prior)
  
  # Linear predictor for the mean of the weight between actors i and j at time t in layer k
  mu_ijtk <- get_linpred_ijtk( baseline_tk=theta_tk,
                               coord_ith=uv_ith_shared,
                               coord_ithk=uv_ithk,
                               directed=TRUE )
  
  
  
  ### LOGISTIC MODEL FOR LINKS ###
  
  # Associated binary adjacency matrix
  z_ijtk <- y_ijtk
  z_ijtk[y_ijtk!=0] <- 1
  
  # Augmented Polya-gamma data
  w_ijtk <- y_ijtk
  w_ijtk[!is.na(w_ijtk)] <- 0
  
  # Baseline parameter for link #
  # at time t for layer k
  eta_tk <- matrix( #data=0,
    data=runif(T_net*K_net),
    nrow=T_net,
    ncol=K_net )
  eta_tk_mcmc <- array( NA, dim=c(T_net,K_net,n_iter_mcmc_out) )
  
  
  # Latent coordinates #
  # shared: hth coordinate of actor v at time t shared across the different layers
  # One latent space for each direction #
  ab_ith_shared <- list( send=array( data=runif(V_net*T_net*H_dim,-1,1),
                                    dim=c(V_net,T_net,H_dim) ),
                        receive=array( data=runif(V_net*T_net*H_dim,-1,1),
                                       dim=c(V_net,T_net,H_dim) ) )
  ab_ith_shared_mcmc <- list( send=array(NA,c(V_net,T_net,H_dim,n_iter_mcmc_out)),
                             receive=array(NA,c(V_net,T_net,H_dim,n_iter_mcmc_out)) )
  
  # by layer: hth coordinate of actor v at time t specific to layer k
  if( K_net>1 ){
    # One latent space for each direction #
    ab_ithk <- list( send=array( data=runif(V_net*T_net*R_dim*K_net),
                                dim=c(V_net,T_net,R_dim,K_net) ),
                    receive=array( data=runif(V_net*T_net*R_dim*K_net),
                                   dim=c(V_net,T_net,R_dim,K_net) ) )
    ab_ithk_mcmc <- list( send=array(NA,c(V_net,T_net,R_dim,K_net,n_iter_mcmc_out)),
                         receive=array(NA,c(V_net,T_net,R_dim,K_net,n_iter_mcmc_out)) )
  } else {
    ab_ithk <- NULL
    ab_ithk_mcmc <- NULL
  }
  
  # Linear predictor for the probability of an edge between actors i and j at time t in layer k
  gamma_ijtk <- get_linpred_ijtk( baseline_tk=eta_tk,
                                  coord_ith=ab_ith_shared,
                                  coord_ithk=ab_ithk,
                                  directed=TRUE )
  
  # Probability of an edge between actors i and j at time t in layer k
  pi_ijtk_mcmc <- array(NA, dim=c(V_net,V_net,T_net,K_net,n_iter_mcmc_out))
  
  #pi_ijt[,,1]
  # all( abs(qlogis( pi_ijt ) - s_ijt)<1e-6,na.rm=T ) # TRUE
  
  ## Shrinkage Parameters ##
  
  # shared #
  v_shrink_shared <- matrix(1, nrow=H_dim, ncol=1 )
  if(shrink_lat_space){
    v_shrink_shared[1,1] <- rgamma(n=1,shape=a_1,rate=1); v_shrink_shared[-1,1] <- rgamma(n=H_dim-1,shape=a_2,rate=1)
  }
  tau_h_shared <- matrix(cumprod(v_shrink_shared), nrow=H_dim, ncol=1 )
  if(shrink_lat_space){
    tau_h_shared_mcmc <- matrix(NA, nrow=H_dim, ncol=n_iter_mcmc_out )
  } else {
    tau_h_shared_mcmc <- NULL
  }
  # duplicate, one latent space for each direction #
  v_shrink_shared <- list( send=v_shrink_shared,
                           receive=v_shrink_shared )
  tau_h_shared <- list( send=tau_h_shared,
                        receive=tau_h_shared )
  tau_h_shared_mcmc <- list( send=tau_h_shared_mcmc,
                             receive=tau_h_shared_mcmc )
  
  # Creates one for the weights and one for the link
  v_shrink_shared_weight <- v_shrink_shared_link <- v_shrink_shared
  rm(v_shrink_shared)
  tau_h_shared_weight <- tau_h_shared_link <- tau_h_shared
  rm(tau_h_shared)
  tau_h_shared_weight_mcmc <- tau_h_shared_link_mcmc <- tau_h_shared_mcmc
  rm(tau_h_shared_mcmc)
  
  # layer-specific #
  if( K_net>1 ){
    v_shrink_k <- matrix(1, nrow=R_dim, ncol=K_net )
    if(shrink_lat_space){
      v_shrink_k[1,] <- rgamma(n=K_net,shape=a_1,rate=1); v_shrink_k[-1,] <- rgamma(n=K_net*(R_dim-1),shape=a_2,rate=1)
    }
    tau_h_k <- matrix(apply(v_shrink_k,2,cumprod), nrow=R_dim, ncol=K_net )
    if(shrink_lat_space){
      tau_h_k_mcmc <- array( NA, dim=c(R_dim,K_net,n_iter_mcmc_out) )
    } else {
      tau_h_k_mcmc <- NULL
    }
    
    # duplicate, one latent space for each direction #
    v_shrink_k <- list( send=v_shrink_k,
                        receive=v_shrink_k )
    tau_h_k <- list( send=tau_h_k,
                     receive=tau_h_k )
    tau_h_k_mcmc <- list( send=tau_h_k_mcmc,
                          receive=tau_h_k_mcmc )
  } else {
    v_shrink_k <- NULL
    tau_h_k <- NULL
    tau_h_k_mcmc <- NULL
  }
  
  # Creates one for the weights and one for the link
  v_shrink_k_weight <- v_shrink_k_link <- v_shrink_k
  rm(v_shrink_k)
  tau_h_k_weight <- tau_h_k_link <- tau_h_k
  rm(tau_h_k)
  tau_h_k_weight_mcmc <- tau_h_k_link_mcmc <- tau_h_k_mcmc
  rm(tau_h_k_mcmc)
  
  
  ### Missing links ###
  y_ijtk_orig <- y_ijtk
  
  y_ijtk_miss <- FALSE
  y_ijtk_miss_idx <- NULL
  y_ijtk_imp_mcmc <- NULL
  
  if( any(is.na(y_ijtk)) ) {
    y_ijtk_miss <- TRUE
    
    # identifies missing data indices
    y_ijtk_miss_idx <- which( is.na(y_ijtk), arr.ind=TRUE )
    colnames(y_ijtk_miss_idx) <- c("i","j","t","k")
    
    # Only consider missing data outside the diagonal
    y_ijtk_miss_idx <- y_ijtk_miss_idx[y_ijtk_miss_idx[,1]!=y_ijtk_miss_idx[,2],]
    
    # MCMC chain for missing values #
    y_ijtk_imp_mcmc <- matrix( NA, nrow = n_iter_mcmc_out, ncol=nrow(y_ijtk_miss_idx) )
    
    # Imputing y_ijtk and z_ijtk #
    if(F){ # WE WILL NOT IMPUTE HERE, as it make harder to converge for the inference in the other parameters
      # First, impute links in z_ijtk #
      # We will initialise the imputation using the average number of links for every pair in a given layer
      for(k in 1:K_net){ # k <-1
        idx_aux <- y_ijtk_miss_idx[,4]==k
        if(sum(idx_aux)>0){
          probs_aux <- apply(z_ijtk[,,,k],1:2,mean,na.rm=T)
          weight_mean_aux <- apply(y_ijtk[,,,k],1:2,mean,na.rm=T)
          weight_sd_aux <- apply(y_ijtk[,,,k],1:2,sd,na.rm=T)
          # Sample link #
          z_ijtk[y_ijtk_miss_idx[idx_aux,]] <- rbinom( n=sum(idx_aux),
                                                            size = 1,
                                                            prob=probs_aux[y_ijtk_miss_idx[idx_aux,1:2]] )
          # Sample weight #
          y_ijtk[y_ijtk_miss_idx[idx_aux,]] <- rnorm( n=sum(idx_aux),
                                                      mean=weight_mean_aux[y_ijtk_miss_idx[idx_aux,1:2]],
                                                      sd=weight_sd_aux[y_ijtk_miss_idx[idx_aux,1:2]] )
          # Mixture #
          y_ijtk[y_ijtk_miss_idx[idx_aux,]] <- apply(cbind(0,y_ijtk[y_ijtk_miss_idx[idx_aux,]] * z_ijtk[y_ijtk_miss_idx[idx_aux,]]),1,max)
          
        }
      }
    }
    
  }
  
  #### End: MCMC initialization ####
  
  
  # Keep execution time #
  mcmc_clock <- Sys.time()
  if( !is.null(log_file) ) {
    cat("MCMC Starting time:\n",as.character(mcmc_clock),"\n\n",
        "---------------------------\n\n\n",
        "iter_i , mcmc_acum_minutes , Sys.time \n",
        file=log_file, append=T )
  }
  
  #### Start: MCMC Sampling ####
  if(!quiet_mcmc){ cat("Sampling MCMC ...\n") }
  
  for ( iter_i in 1:n_iter_mcmc) { # iter_i <- 1
    #cat(iter_i,",")
    
    ##### SAMPLING WEIGHTS #####
    
    ### Step W1. Sample theta_tk from its conditional N-variate Gaussian posterior ###
    out_aux <- sample_baseline_tk_weight( theta_tk=theta_tk,
                                          y_ijtk=y_ijtk, mu_ijtk=mu_ijtk,
                                          sigma_k=sigma_k,
                                          theta_t_cov_prior_inv=cov_gp_prior_inv,
                                          directed=TRUE )
    theta_tk <- out_aux$theta_tk
    mu_ijtk <- out_aux$mu_ijtk # This updates mu_ijtk except for the diagonal
    mu_ijtk[diag_y_idx] <- NA
    
    # MCMC chain #
    if(is.element(iter_i,iter_out_mcmc)){
      theta_tk_mcmc[,,match(iter_i,iter_out_mcmc)] <- theta_tk
    }
    
    # mu_ijtk_new <- get_linpred_ijtk( baseline_tk=theta_tk,
    #                                  coord_ith=uv_ith_shared,
    #                                  coord_ithk=uv_ithk,
    #                                  directed=TRUE )
    # mu_ijtk_new[diag_y_idx] <- NA
    # all.equal(mu_ijtk,mu_ijtk_new)
    
    ### Step W2. For each unit, block-sample the set of time-varying latent coordinates x_ith ###
    ### SHARED Latent Coordinates ###
    out_aux <- sample_coord_ith_shared_weight( uv_ith_shared=uv_ith_shared,
                                               uv_t_sigma_prior_inv=cov_gp_prior_inv,
                                               tau_h=tau_h_shared_weight,
                                               y_ijtk=y_ijtk, mu_ijtk=mu_ijtk,
                                               sigma_k=sigma_k,
                                               directed=TRUE )
    uv_ith_shared <- out_aux$uv_ith_shared
    mu_ijtk <- out_aux$mu_ijtk
    
    # MCMC chain #
    if(is.element(iter_i,iter_out_mcmc)){
      uv_ith_shared_mcmc[[1]][,,,match(iter_i,iter_out_mcmc)] <- uv_ith_shared[[1]]
      uv_ith_shared_mcmc[[2]][,,,match(iter_i,iter_out_mcmc)] <- uv_ith_shared[[2]]
    }
    
    
    
    ### LAYER SPECIFIC Latent Coordinates ###
    if( K_net>1 ) {
      ### Step 3A. For each unit, block-sample the EDGE SPECIFIC set of time-varying latent coordinates x_ithk ###
      out_aux <- sample_coord_ithk_weight( uv_ithk=uv_ithk,
                                           uv_t_sigma_prior_inv=cov_gp_prior_inv,
                                           tau_h=tau_h_k,
                                           y_ijtk=y_ijtk, mu_ijtk=mu_ijtk,
                                           sigma_k=sigma_k,
                                           directed=TRUE )
      uv_ithk <- out_aux$uv_ithk
      mu_ijtk <- out_aux$mu_ijtk
      # mu_ijtk[diag_y_idx] <- NA
      
      # MCMC chain for x_ithk #
      if(is.element(iter_i,iter_out_mcmc)){
        uv_ithk_mcmc[[1]][,,,,match(iter_i,iter_out_mcmc)] <- uv_ithk[[1]]
        uv_ithk_mcmc[[2]][,,,,match(iter_i,iter_out_mcmc)] <- uv_ithk[[2]]
      }
      
    }
    
    ### Step W3. Sample weight variance ###
    # TO BE IMPLEMENTED #
    
    
    
    ##### SAMPLING LINKS #####
    
    ### Step L1. Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
    w_ijtk <- sample_w_ijtk_DynMultiNet_bin( w_ijtk=w_ijtk,
                                             s_ijtk=gamma_ijtk,
                                             directed=TRUE )
    
    
    
    ### Step L2_mu. Sample eta_tk from its conditional N-variate Gaussian posterior ###
    out_aux <- sample_mu_tk_DynMultiNet_bin( mu_tk=eta_tk,
                                             y_ijtk=y_ijtk, w_ijtk=w_ijtk, s_ijtk=gamma_ijtk,
                                             mu_t_cov_prior_inv=cov_gp_prior_inv,
                                             directed=TRUE )
    eta_tk <- out_aux$mu_tk
    gamma_ijtk <- out_aux$s_ijtk # This updates gamma_ijtk except for the diagonal
    gamma_ijtk[diag_y_idx] <- NA
    
    # MCMC chain #
    if(is.element(iter_i,iter_out_mcmc)){
      eta_tk_mcmc[,,match(iter_i,iter_out_mcmc)] <- eta_tk
    }
    
    
    ### Step L2_beta. Sample beta_z_layer and beta_z_edge from its conditional N-variate Gaussian posterior ###
    # TO BE IMPLEMENTED #
    
    ### Step L3. For each unit, block-sample the set of time-varying latent coordinates x_ith ###
    ### SHARED Latent Coordinates ###
    out_aux <- sample_x_ith_shared_DynMultiNet_bin( x_ith_shared=x_ith_shared,
                                                    x_t_sigma_prior_inv=x_t_sigma_prior_inv,
                                                    tau_h=tau_h_shared,
                                                    y_ijtk=y_ijtk,
                                                    w_ijtk=w_ijtk,
                                                    s_ijtk=gamma_ijtk,
                                                    directed=TRUE )
    x_ith_shared <- out_aux$x_ith_shared
    gamma_ijtk <- out_aux$s_ijtk
    # gamma_ijtk[diag_y_idx] <- NA
    
    # MCMC chain #
    if(is.element(iter_i,iter_out_mcmc)){
      x_ith_shared_mcmc[[1]][,,,match(iter_i,iter_out_mcmc)] <- x_ith_shared[[1]]
      x_ith_shared_mcmc[[2]][,,,match(iter_i,iter_out_mcmc)] <- x_ith_shared[[2]]
    }
    
    
    
    ### LAYER SPECIFIC Latent Coordinates ###
    if( K_net>1 ) {
      ### Step 3A. For each unit, block-sample the EDGE SPECIFIC set of time-varying latent coordinates x_ithk ###
      out_aux <- sample_x_ithk_DynMultiNet_bin( x_ithk=x_ithk,
                                                x_t_sigma_prior_inv=x_t_sigma_prior_inv,
                                                tau_h=tau_h_k,
                                                y_ijtk=y_ijtk, w_ijtk=w_ijtk, s_ijtk=gamma_ijtk,
                                                directed=TRUE )
      x_ithk <- out_aux$x_ithk
      gamma_ijtk <- out_aux$s_ijtk
      # gamma_ijtk[diag_y_idx] <- NA
      
      # MCMC chain for x_ithk #
      if(is.element(iter_i,iter_out_mcmc)){
        x_ithk_mcmc[[1]][,,,,match(iter_i,iter_out_mcmc)] <- x_ithk[[1]]
        x_ithk_mcmc[[2]][,,,,match(iter_i,iter_out_mcmc)] <- x_ithk[[2]]
      }
      
      
    }
    
    
    
    ### Edge probabilities ###
    if(is.element(iter_i,iter_out_mcmc)){
      gamma_ijtk[diag_y_idx] <- NA
      pi_ijtk_mcmc[,,,,match(iter_i,iter_out_mcmc)] <- plogis(gamma_ijtk)
    }
    
    
    
    ### Impute missing links ###
    if( y_ijtk_miss ) {
      y_ijtk[y_ijtk_miss_idx] <- rbinom( n=nrow(y_ijtk_miss_idx), size=1, prob=plogis(gamma_ijtk[y_ijtk_miss_idx]) )
      
      # MCMC chain #
      if(is.element(iter_i,iter_out_mcmc)){
        y_ijtk_imp_mcmc[match(iter_i,iter_out_mcmc),] <- y_ijtk[y_ijtk_miss_idx]
      }
    }
    
    
    
    ### Step L4. Sample the global shrinkage hyperparameters from conditional gamma distributions ###
    if(shrink_lat_space){
      for( dir_i in 1:2 ) {
        v_shrink_shared[[dir_i]] <- sample_v_shrink_DynMultiNet_bin( v_shrink_shared[[dir_i]],
                                                                     a_1, a_2,
                                                                     x_ith_shared[[dir_i]],
                                                                     x_t_sigma_prior_inv )
        tau_h_shared[[dir_i]] <- matrix(cumprod(v_shrink_shared[[dir_i]]), nrow=H_dim, ncol=1 )
        if(is.element(iter_i,iter_out_mcmc)){
          tau_h_shared_mcmc[[dir_i]][,match(iter_i,iter_out_mcmc)] <- tau_h_shared[[dir_i]]
        }
      }
      
      if(K_net>1){
        for( dir_i in 1:2 ) {
          for(k in 1:K_net) {
            v_shrink_k[[dir_i]][,k] <- sample_v_shrink_DynMultiNet_bin( v_shrink_k[[dir_i]][,k,drop=F],
                                                                        a_1, a_2,
                                                                        x_ithk[[dir_i]][,,,k],
                                                                        x_t_sigma_prior_inv )
          }
          tau_h_k[[dir_i]] <- matrix(apply(v_shrink_k[[dir_i]],2,cumprod), nrow=R_dim, ncol=K_net )
          if(is.element(iter_i,iter_out_mcmc)){
            tau_h_k_mcmc[[dir_i]][,,match(iter_i,iter_out_mcmc)] <- tau_h_k[[dir_i]]
          }
        }
      }
    }
    
    # display MCMC progress #
    if( is.element(iter_i, floor(n_iter_mcmc*seq(0,1,0.05)[-1]) ) ) {
      if(!quiet_mcmc){
        cat(round(100*iter_i/n_iter_mcmc),"% ",sep="")
      }
      if( !is.null(log_file) ) {
        cat(iter_i," , ", as.numeric(difftime(Sys.time(),mcmc_clock,units="mins"))," , ", as.character(Sys.time()),"\n",
            file=log_file,append=TRUE )
      }
    }
    
    # save MCMC progress #
    if( is.element(iter_i, floor(n_iter_mcmc*seq(0,1,0.25)[-1]) ) & iter_i>min(iter_out_mcmc,na.rm=T) ) {
      if(!is.null(rds_file)){
        dmn_mcmc <- list( y_ijtk=y_ijtk_orig,
                          
                          directed=TRUE,
                          weighted=FALSE,
                          
                          n_chains_mcmc=NULL,
                          n_iter_mcmc=n_iter_mcmc, n_burn=n_burn, n_thin=n_thin,
                          
                          H_dim=H_dim, R_dim=R_dim,
                          k_x=k_x, k_mu=k_mu, k_p=k_p,
                          
                          shrink_lat_space=shrink_lat_space,
                          a_1=a_1, a_2=a_2,
                          
                          procrustes_lat=procrustes_lat,
                          
                          node_all=node_all, time_all=time_all, layer_all=layer_all,
                          
                          # For link probabilities #
                          pi_ijtk_mcmc=pi_ijtk_mcmc,
                          eta_tk_mcmc=eta_tk_mcmc,
                          x_ith_shared_mcmc=x_ith_shared_mcmc,
                          x_ithk_mcmc=x_ithk_mcmc,
                          tau_h_shared_mcmc=tau_h_shared_mcmc,
                          tau_h_k_mcmc=tau_h_k_mcmc,
                          
                          # imputed links
                          y_ijtk_miss_idx=y_ijtk_miss_idx,
                          y_ijtk_imp_mcmc=y_ijtk_imp_mcmc,
                          
                          # For link weights #
                          r_ijtk_mcmc = NULL,
                          sigma_w_k_mcmc = NULL,
                          lambda_tk_mcmc = NULL,
                          u_ith_shared_mcmc = NULL,
                          u_ithk_mcmc = NULL,
                          rho_h_shared_mcmc = NULL,
                          rho_h_k_mcmc = NULL,
                          
                          pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                          beta_z_layer_mcmc=beta_z_layer_mcmc,
                          beta_z_edge_mcmc=beta_z_edge_mcmc )
        dmn_mcmc <- structure( dmn_mcmc, class="dmn_mcmc" )
        saveRDS(dmn_mcmc,file=rds_file)
      }
    }
  }
  if(!quiet_mcmc){cat("\nMCMC finished!\n")}
  #### End: MCMC Sampling ####
  if( !is.null(log_file) ) {
    cat("\n\n---------------------------\n\n",
        "Finishing time:\n",as.character(Sys.time()),"\n\n",
        file=log_file, append=TRUE)
  }
  
  dmn_mcmc <- list( y_ijtk=y_ijtk_orig,
                    
                    directed=TRUE,
                    weighted=FALSE,
                    
                    n_chains_mcmc=NULL,
                    n_iter_mcmc=n_iter_mcmc, n_burn=n_burn, n_thin=n_thin,
                    
                    H_dim=H_dim, R_dim=R_dim,
                    k_x=k_x, k_mu=k_mu, k_p=k_p,
                    
                    shrink_lat_space=shrink_lat_space,
                    a_1=a_1, a_2=a_2,
                    
                    procrustes_lat=procrustes_lat,
                    
                    node_all=node_all, time_all=time_all, layer_all=layer_all,
                    
                    # For link probabilities #
                    pi_ijtk_mcmc=pi_ijtk_mcmc,
                    eta_tk_mcmc=eta_tk_mcmc,
                    x_ith_shared_mcmc=x_ith_shared_mcmc,
                    x_ithk_mcmc=x_ithk_mcmc,
                    tau_h_shared_mcmc=tau_h_shared_mcmc,
                    tau_h_k_mcmc=tau_h_k_mcmc,
                    
                    # imputed links
                    y_ijtk_miss_idx=y_ijtk_miss_idx,
                    y_ijtk_imp_mcmc=y_ijtk_imp_mcmc,
                    
                    # For link weights #
                    r_ijtk_mcmc = NULL,
                    sigma_w_k_mcmc = NULL,
                    lambda_tk_mcmc = NULL,
                    u_ith_shared_mcmc = NULL,
                    u_ithk_mcmc = NULL,
                    rho_h_shared_mcmc = NULL,
                    rho_h_k_mcmc = NULL,
                    
                    pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                    beta_z_layer_mcmc=beta_z_layer_mcmc,
                    beta_z_edge_mcmc=beta_z_edge_mcmc )
  
  dmn_mcmc <- structure( dmn_mcmc, class="dmn_mcmc" )
  
  return( dmn_mcmc )
  
}
