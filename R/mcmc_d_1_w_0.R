#' @title
#'    MCMC algorithm for Dynamic Multilayer directed unweighted graphs
#'
#' @description
#'    \code{mcmc_d_1_w_0} Implements a Gibbs sampler MCMC algorithm for Dynamic Multilayer directed unweighted graphs.
#'
#' @param y_ijtk Array. Network data, with entry \code{y_ijtk[i,j,t,k]} representing the link from node i to node j at time t in layer k. 0 indicates no link.
#' @param node_all Character vector. Id's of nodes in the network.
#' @param time_all Numeric vector. Timestamps of all relevant epochs for the MCMC, those observed and those for forecast
#' @param layer_all Character vector. Id's of layers in the network.
#' @param x_ijtkp Array. Edge Specific external covariates.
#' @param H_dim Integer. Latent space dimension.
#' @param R_dim Integer. Latent space dimension, for layer specific latent vectors.
#' @param add_eff_link Boolean. Indicates if dynamic additive effects by node should be considered.
#' @param class_dyn String. Specifies prior for latent elements. "GP" for Gaussian Processes,"nGP" for nested Gaussian Processes.
#' @param delta Positive scalar. Hyperparameter controlling for the smoothness in the dynamic of latent coordinates. Larger=smoother.
#' @param n_chains_mcmc Integer. Number of parallel MCMC chains.
#' @param n_iter_mcmc Integer. Number of iterations for the MCMC.
#' @param n_burn Integer. Number of iterations discarded as part of the MCMC warming up period at the beginning of the chain.
#' @param n_thin Integer. Number of iterations discarded for thining the chain (reducing the autocorrelation). We keep 1 of every n_thin iterations.
#' @param keep_y_ijtk_imp Boolean. Indicates wheter the chain with imputed missing links will be saved (FALSE by default)
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


mcmc_d_1_w_0 <- function( y_ijtk,
                          node_all, time_all, layer_all,
                          
                          x_ijtkp=NULL,
                          
                          H_dim=10, R_dim=10,
                          
                          add_eff_link=FALSE,
                          
                          class_dyn=c("GP","nGP")[1],
                          
                          delta=1,
                          
                          n_chains_mcmc=1,
                          n_iter_mcmc=10000, n_burn=floor(n_iter_mcmc/4), n_thin=3,
                          
                          keep_y_ijtk_imp=FALSE,
                          
                          rds_file=NULL, log_file=NULL,
                          quiet_mcmc=FALSE,
                          parallel_mcmc=FALSE ) {
  
  
  shrink_lat_space=FALSE
  a_1=2; a_2=2.5
  procrustes_lat=FALSE
  
  
  # This function only deal with binary edges (non-weighted)
  y_ijtk[y_ijtk>0] <- 1
  
  V_net <- length(node_all)
  T_net <- length(time_all)
  K_net <- length(layer_all)
  
  # nGP matrices approximated or exact (nGP_mat_approx=F means exact)
  nGP_mat_approx = F
  
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
  
  # Augmented Polya-gamma data
  w_ijtk <- y_ijtk
  w_ijtk[!is.na(w_ijtk)] <- 0
  
  # Baseline parameter #
  # at time t for layer k
  eta_tk <- matrix( runif(T_net*K_net), T_net, K_net )
  eta_tk_mcmc <- array( NA, dim=c(T_net,K_net,n_iter_mcmc_out) )
  if( class_dyn=="nGP" ){
    alpha_eta_tk <- array( runif(T_net*K_net*3), dim=c(T_net,K_net,3) )
    alpha_eta_tk[,,1] <- eta_tk
  }
  
  # Latent coordinates #
  # shared: hth coordinate of actor v at time t shared across the different layers
  # One latent space for each direction #
  ab_ith_shared <- list( send=array( runif(V_net*T_net*H_dim,-1,1), dim=c(V_net,T_net,H_dim) ),
                         receive=array( runif(V_net*T_net*H_dim,-1,1), dim=c(V_net,T_net,H_dim) ) )
  if( class_dyn=="nGP" ){
    alpha_ab_ith_shared <- list( send=array( runif(V_net*T_net*H_dim,-1,1), dim=c(V_net,T_net,H_dim,3) ),
                                 receive=array( runif(V_net*T_net*H_dim,-1,1), dim=c(V_net,T_net,H_dim,3) ) )
    alpha_ab_ith_shared[[1]][,,,1] <- ab_ith_shared[[1]]
    alpha_ab_ith_shared[[2]][,,,1] <- ab_ith_shared[[2]]
  }
  ab_ith_shared_mcmc <- list( send=array(NA,c(V_net,T_net,H_dim,n_iter_mcmc_out)),
                              receive=array(NA,c(V_net,T_net,H_dim,n_iter_mcmc_out)) )
  
  # by layer: hth coordinate of actor v at time t specific to layer k
  ab_ithk <- NULL
  ab_ithk_mcmc <- NULL
  alpha_ab_ithk <- NULL
  if( K_net>1 ){
    # One latent space for each direction #
    ab_ithk <- list( send=array( runif(V_net*T_net*R_dim*K_net), dim=c(V_net,T_net,R_dim,K_net) ),
                     receive=array( runif(V_net*T_net*R_dim*K_net), dim=c(V_net,T_net,R_dim,K_net) ) )
    ab_ithk_mcmc <- list( send=array(NA,c(V_net,T_net,R_dim,K_net,n_iter_mcmc_out)),
                         receive=array(NA,c(V_net,T_net,R_dim,K_net,n_iter_mcmc_out)) )
    if( class_dyn=="nGP" ){
      alpha_ab_ithk <- list( send=array( runif(V_net*T_net*R_dim*K_net*3,-0.1,0.1), dim=c(V_net,T_net,R_dim,K_net,3) ),
                             receive=array( runif(V_net*T_net*R_dim*K_net*3,-0.1,0.1), dim=c(V_net,T_net,R_dim,K_net,3) ) )
      alpha_ab_ithk[[1]][,,,,1] <- ab_ithk[[1]]
      alpha_ab_ithk[[2]][,,,,1] <- ab_ithk[[2]]
    }
  }
  
  ### Dynamic additive effects for each node ###
  sp_link_it_shared <- NULL
  sp_link_it_shared_mcmc <- NULL
  alpha_sp_link_it_shared<-NULL
  sp_link_itk <- NULL
  sp_link_itk_mcmc <- NULL
  alpha_sp_link_itk <- NULL
  if(add_eff_link){
    sp_link_it_shared <- array(0,dim=c(V_net,T_net,2))
    sp_link_it_shared_mcmc <- array(NA,dim=c(V_net,T_net,2,n_iter_mcmc_out))
    if( class_dyn=="nGP" ){
      alpha_sp_link_it_shared <- array( 0, dim=c(V_net,T_net,2,3) )
      alpha_sp_link_it_shared[,,,1] <- sp_link_it_shared
    }
    if(FALSE&(K_net>1)){ # Only consider global additive effects
      sp_link_itk <- array(0,dim=c(V_net,T_net,K_net,2))
      sp_link_itk_mcmc <- array(NA,dim=c(V_net,T_net,K_net,2,n_iter_mcmc_out))
      if( class_dyn=="nGP" ){
        alpha_sp_link_itk <- array( 0, dim=c(V_net,T_net,K_net,2,3) )
        alpha_sp_link_itk[,,,,1] <- sp_link_itk
      }
    }
  }
  
  # Coefficients for external covariates #
  beta_lambda_tp <- NULL
  beta_lambda_tp_mcmc <- NULL
  if(!is.null(x_ijtkp)) {
    P_pred <- dim(x_ijtkp)[5]
    beta_lambda_tp <- matrix( runif(T_net*P_pred), nrow=T_net, ncol=P_pred )
    beta_lambda_tp_mcmc <- array( NA, dim=c(T_net,P_pred,n_iter_mcmc_out) )
    x_ijtkp_mat <- get_x_ijtkp_mat( x_ijtkp=x_ijtkp,
                                    directed=FALSE,
                                    weighted=FALSE )
    if( class_dyn=="nGP" ){
      alpha_beta_lambda_tp <- array( 0, dim=c(V_net,T_net,K_net,3) )
      alpha_beta_lambda_tp[,,,1] <- beta_lambda_tp
    }
  }
  
  # Linear predictor for the probability of an edge between actors i and j at time t in layer k
  gamma_ijtk <- get_linpred_ijtk( baseline_tk=eta_tk,
                                  coord_ith_shared=ab_ith_shared,
                                  coord_ithk=ab_ithk,
                                  add_eff_it_shared=sp_link_it_shared,
                                  add_eff_itk=sp_link_itk,
                                  directed=TRUE )
  gamma_ijtk[diag_y_idx] <- NA
  
  # Probability of an edge between actors i and j at time t in layer k
  pi_ijtk_mcmc <- array(NA, dim=c(V_net,V_net,T_net,K_net,n_iter_mcmc_out))
  
  
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
  
  # Covariance in Latent dynamic elements #
  cov_gp_prior <- cov_gp_prior_inv <- NULL
  nGP_sigma_net <- nGP_mat <- NULL
  if( class_dyn=="GP" ){
    # Covariance matrix prior for latent coordinates
    cov_gp_prior <- outer( time_all, time_all, FUN=function(x,y,k=delta){ exp(-((x-y)/delta)^2) } )
    diag(cov_gp_prior) <- diag(cov_gp_prior) + 1e-3 # numerical stability
    cov_gp_prior_inv <- solve(cov_gp_prior)
  } else if( class_dyn=="nGP" ){
    # Hyperparameters for variance of nGPs
    a <- b <- 0.01
    nGP_sigma_net <- get_nGP_sigma_net( a=a, b=b,
                                        time_all=time_all,
                                        alpha_baseline_tk=alpha_eta_tk,
                                        alpha_coord_ith_shared=alpha_ab_ith_shared,
                                        alpha_coord_ithk=alpha_ab_ithk,
                                        alpha_add_eff_it_shared=alpha_sp_link_it_shared,
                                        alpha_add_eff_itk=alpha_sp_link_itk,
                                        directed=TRUE )
    nGP_mat <- get_nGP_mat_net( nGP_sigma_net,
                                nGP_mat_approx=nGP_mat_approx,
                                directed=TRUE )
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
  
  ### Missing links ###
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
    if(keep_y_ijtk_imp){
      # CAUTION: may take too much disk space
      y_ijtk_imp_mcmc <- matrix( NA, nrow = n_iter_mcmc_out, ncol=nrow(y_ijtk_miss_idx) )
    } else {
      y_ijtk_imp_mcmc <- NULL 
    }
  }
  
  #### Start: MCMC Sampling ####
  if(!quiet_mcmc){ cat("Sampling MCMC ...\n") }
  
  for ( iter_i in 1:n_iter_mcmc) { # iter_i <- 1
    #cat(iter_i,",")
    
    
    
    ### Step 1. Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
    w_ijtk <- sample_pg_w_ijtk_link( w_ijtk=w_ijtk,
                                     gamma_ijtk=gamma_ijtk,
                                     directed=TRUE )
    
    
    ### Step 2. Sample eta_tk from its conditional N-variate Gaussian posterior ###
    out_aux <- sample_baseline_tk_link( eta_tk=eta_tk,
                                        y_ijtk=y_ijtk, w_ijtk=w_ijtk, gamma_ijtk=gamma_ijtk,
                                        class_dyn=class_dyn,
                                        eta_t_cov_prior_inv=cov_gp_prior_inv,
                                        nGP_mat=nGP_mat$baseline_k,
                                        alpha_eta_tk=alpha_eta_tk,
                                        directed=TRUE )
    eta_tk <- out_aux$eta_tk
    gamma_ijtk <- out_aux$gamma_ijtk # This updates gamma_ijtk except for the diagonal
    
    
    eta_tk <- out_aux$eta_tk
    gamma_ijtk <- out_aux$gamma_ijtk # This updates ONLY the lower triangular matrices in gamma_ijtk
    gamma_ijtk[diag_y_idx] <- NA
    if( class_dyn=="nGP" ){
      alpha_eta_tk <- out_aux$alpha_eta_tk
    }
    rm(out_aux)
    
    # # Checking consistency of linear predictor gamma_ijtk
    # gamma_ijtk_old <- gamma_ijtk
    # gamma_ijtk_old[diag_y_idx] <- NA
    # gamma_ijtk <- get_linpred_ijtk( baseline_tk=eta_tk,
    #                                 add_eff_it_shared=sp_link_it_shared,
    #                                 add_eff_itk=sp_link_itk,
    #                                 coord_ith_shared=ab_ith_shared,
    #                                 coord_ithk=ab_ithk,
    #                                 directed=TRUE )
    # gamma_ijtk[diag_y_idx] <- NA
    # all.equal(gamma_ijtk,gamma_ijtk_old)
    
    # MCMC chain #
    if(is.element(iter_i,iter_out_mcmc)){
      eta_tk_mcmc[,,match(iter_i,iter_out_mcmc)] <- eta_tk
    }
    
    ### Step 3. For each unit, block-sample the set of time-varying latent coordinates x_ith ###
    ### SHARED Latent Coordinates ###
    out_aux <- sample_coord_ith_shared_link( ab_ith=ab_ith_shared,
                                             y_ijtk=y_ijtk, w_ijtk=w_ijtk, gamma_ijtk=gamma_ijtk,
                                             class_dyn=class_dyn,
                                             ab_t_sigma_prior_inv=cov_gp_prior_inv,
                                             tau_h=tau_h_shared,
                                             alpha_ab_ith=alpha_ab_ith_shared,
                                             nGP_mat=nGP_mat$coord_i,
                                             directed=TRUE )
    ab_ith_shared <- out_aux$ab_ith
    gamma_ijtk <- out_aux$gamma_ijtk
    # gamma_ijtk[diag_y_idx] <- NA
    if( class_dyn=="nGP" ){
      alpha_ab_ith_shared <- out_aux$alpha_ab_ith
    }
    rm(out_aux)
    
    # MCMC chain #
    if(is.element(iter_i,iter_out_mcmc)){
      ab_ith_shared_mcmc[[1]][,,,match(iter_i,iter_out_mcmc)] <- ab_ith_shared[[1]]
      ab_ith_shared_mcmc[[2]][,,,match(iter_i,iter_out_mcmc)] <- ab_ith_shared[[2]]
    }
    
    ### LAYER SPECIFIC Latent Coordinates ###
    if( K_net>1 ) {
      ### Step 3A. For each unit, block-sample the EDGE SPECIFIC set of time-varying latent coordinates ab_ithk ###
      out_aux <- sample_coord_ithk_link( ab_ithk=ab_ithk,
                                         y_ijtk=y_ijtk, w_ijtk=w_ijtk, gamma_ijtk=gamma_ijtk,
                                         class_dyn=class_dyn,
                                         ab_t_sigma_prior_inv=cov_gp_prior_inv,
                                         tau_h=tau_h_k,
                                         alpha_ab_ithk=alpha_ab_ithk,
                                         nGP_mat=nGP_mat$coord_ik,
                                         directed=TRUE )
      ab_ithk <- out_aux$ab_ithk
      gamma_ijtk <- out_aux$gamma_ijtk
      # gamma_ijtk[diag_y_idx] <- NA
      if( class_dyn=="nGP" ){
        alpha_ab_ithk <- out_aux$alpha_ab_ithk
      }
      
      # MCMC chain for ab_ithk #
      if(is.element(iter_i,iter_out_mcmc)){
        ab_ithk_mcmc[[1]][,,,,match(iter_i,iter_out_mcmc)] <- ab_ithk[[1]]
        ab_ithk_mcmc[[2]][,,,,match(iter_i,iter_out_mcmc)] <- ab_ithk[[2]]
      }
    }
    
    
    ### Step 4 Sample global additive effects ###
    if(add_eff_link){
      out_aux <- sample_add_eff_it_shared_link( sp_it_shared=sp_link_it_shared,
                                                sp_t_cov_prior_inv=cov_gp_prior_inv,
                                                y_ijtk=y_ijtk, w_ijtk=w_ijtk, gamma_ijtk=gamma_ijtk,
                                                directed=TRUE )
      sp_link_it_shared <- out_aux$sp_it_shared
      gamma_ijtk <- out_aux$gamma_ijtk
      
      # MCMC chain #
      if(is.element(iter_i,iter_out_mcmc)){
        sp_link_it_shared_mcmc[,,,match(iter_i,iter_out_mcmc)] <- sp_link_it_shared
      }
      
      
      ### Step L2_add_eff_itk. Sample layer-specific additive effects ###
      if(FALSE&!is.null(sp_link_itk)){ # Only consider global additive effects
        out_aux <- sample_add_eff_itk_link( sp_itk=sp_link_itk,
                                            sp_t_cov_prior_inv=cov_gp_prior_inv,
                                            y_ijtk=y_ijtk, w_ijtk=w_ijtk, gamma_ijtk=gamma_ijtk,
                                            directed=TRUE )
        sp_link_itk <- out_aux$sp_itk
        gamma_ijtk <- out_aux$gamma_ijtk
        
        # MCMC chain #
        if(is.element(iter_i,iter_out_mcmc)){
          sp_link_itk_mcmc[,,,,match(iter_i,iter_out_mcmc)] <- sp_link_itk
        }
      }
    }
    
    
    ### Step 6A. Compute Edge probabilities ###
    if(is.element(iter_i,iter_out_mcmc)){
      gamma_ijtk[diag_y_idx] <- NA
      pi_ijtk_mcmc[,,,,match(iter_i,iter_out_mcmc)] <- plogis(gamma_ijtk)
    }
    
    ### Step 6B. Impute missing links ###
    if( keep_y_ijtk_imp & y_ijtk_miss ) { # requires too much disk memory
      # MCMC chain #
      if(is.element(iter_i,iter_out_mcmc)){
        Y_imp <- rbinom( n=nrow(y_ijtk_miss_idx), size=1, prob=plogis(gamma_ijtk[y_ijtk_miss_idx]) )
        y_ijtk_imp_mcmc[match(iter_i,iter_out_mcmc),] <- Y_imp
      }
    }
    
    ### Step 7. Sample variance terms ###
    
    if( class_dyn=="GP" ){
      # Sample the global shrinkage hyperparameters from conditional gamma distributions #
      if(shrink_lat_space){
        for( dir_i in 1:2 ) {
          v_shrink_shared[[dir_i]] <- sample_v_shrink_DynMultiNet_bin( v_shrink_shared[[dir_i]],
                                                                       a_1, a_2,
                                                                       ab_ith_shared[[dir_i]],
                                                                       cov_gp_prior_inv )
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
                                                                          ab_ithk[[dir_i]][,,,k],
                                                                          cov_gp_prior_inv )
            }
            tau_h_k[[dir_i]] <- matrix(apply(v_shrink_k[[dir_i]],2,cumprod), nrow=R_dim, ncol=K_net )
            if(is.element(iter_i,iter_out_mcmc)){
              tau_h_k_mcmc[[dir_i]][,,match(iter_i,iter_out_mcmc)] <- tau_h_k[[dir_i]]
            }
          }
        }
      }
    } else if( class_dyn=="nGP" ){
      nGP_sigma_net <- get_nGP_sigma_net( a=a, b=b,
                                          time_all=time_all,
                                          alpha_baseline_tk=alpha_eta_tk,
                                          alpha_coord_ith_shared=alpha_ab_ith_shared,
                                          alpha_coord_ithk=alpha_ab_ithk,
                                          alpha_add_eff_it_shared=alpha_sp_link_it_shared,
                                          alpha_add_eff_itk=alpha_sp_link_itk,
                                          directed=TRUE )
      nGP_mat <- get_nGP_mat_net( nGP_sigma_net,
                                  nGP_mat_approx=nGP_mat_approx,
                                  directed=TRUE )
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
        dmn_mcmc <- list( y_ijtk=y_ijtk,
                          
                          directed=TRUE,
                          weighted=FALSE,
                          
                          n_chains_mcmc=n_chains_mcmc,
                          n_iter_mcmc=n_iter_mcmc, n_burn=n_burn, n_thin=n_thin,
                          
                          H_dim=H_dim, R_dim=R_dim,
                          delta=delta,
                          
                          shrink_lat_space=shrink_lat_space,
                          a_1=a_1, a_2=a_2,
                          
                          procrustes_lat=procrustes_lat,
                          
                          node_all=node_all, time_all=time_all, layer_all=layer_all,
                          
                          # For link probabilities #
                          pi_ijtk_mcmc=pi_ijtk_mcmc,
                          eta_tk_mcmc=eta_tk_mcmc,
                          sp_link_it_shared_mcmc=sp_link_it_shared_mcmc,
                          sp_link_itk_mcmc=sp_link_itk_mcmc,
                          ab_ith_shared_mcmc=ab_ith_shared_mcmc,
                          ab_ithk_mcmc=ab_ithk_mcmc,
                          tau_h_shared_mcmc=tau_h_shared_mcmc,
                          tau_h_k_mcmc=tau_h_k_mcmc,
                          
                          # imputed links
                          y_ijtk_miss_idx=y_ijtk_miss_idx,
                          y_ijtk_imp_mcmc=y_ijtk_imp_mcmc,
                          
                          # For link weights #
                          mu_ijtk_mcmc = NULL,
                          sp_weight_it_shared_mcmc=NULL,
                          sp_weight_itk_mcmc=NULL,
                          sigma_k_mcmc = NULL,
                          theta_tk_mcmc = NULL,
                          uv_ith_shared_mcmc = NULL,
                          uv_ithk_mcmc = NULL,
                          rho_h_shared_mcmc = NULL,
                          rho_h_k_mcmc = NULL
                          
                          )
        
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
  
  dmn_mcmc <- list( y_ijtk=y_ijtk,
                    
                    directed=TRUE,
                    weighted=FALSE,
                    
                    n_chains_mcmc=n_chains_mcmc,
                    n_iter_mcmc=n_iter_mcmc, n_burn=n_burn, n_thin=n_thin,
                    
                    H_dim=H_dim, R_dim=R_dim,
                    delta=delta,
                    
                    shrink_lat_space=shrink_lat_space,
                    a_1=a_1, a_2=a_2,
                    
                    procrustes_lat=procrustes_lat,
                    
                    node_all=node_all, time_all=time_all, layer_all=layer_all,
                    
                    # For link probabilities #
                    pi_ijtk_mcmc=pi_ijtk_mcmc,
                    eta_tk_mcmc=eta_tk_mcmc,
                    sp_link_it_shared_mcmc=sp_link_it_shared_mcmc,
                    sp_link_itk_mcmc=sp_link_itk_mcmc,
                    ab_ith_shared_mcmc=ab_ith_shared_mcmc,
                    ab_ithk_mcmc=ab_ithk_mcmc,
                    tau_h_shared_mcmc=tau_h_shared_mcmc,
                    tau_h_k_mcmc=tau_h_k_mcmc,
                    
                    # imputed links
                    y_ijtk_miss_idx=y_ijtk_miss_idx,
                    y_ijtk_imp_mcmc=y_ijtk_imp_mcmc,
                    
                    # For link weights #
                    mu_ijtk_mcmc = NULL,
                    sp_weight_it_shared_mcmc=NULL,
                    sp_weight_itk_mcmc=NULL,
                    sigma_k_mcmc = NULL,
                    theta_tk_mcmc = NULL,
                    uv_ith_shared_mcmc = NULL,
                    uv_ithk_mcmc = NULL,
                    rho_h_shared_mcmc = NULL,
                    rho_h_k_mcmc = NULL )
  
  dmn_mcmc <- structure( dmn_mcmc, class="dmn_mcmc" )
  
  return( dmn_mcmc )
  
}
