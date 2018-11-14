
#' @keywords internal
mcmc_d_1_w_0 <- function( y_ijtk,
                          node_all, time_all, layer_all,
                          time_all_idx_net,
                          
                          pred_all,
                          pred_id_layer, pred_id_edge,
                          z_tkp, z_ijtkp,
                          
                          H_dim=10, R_dim=10,
                          k_x=0.10, k_mu=0.10, k_p=0.10,
                          a_1=2, a_2=2.5,
                          
                          n_iter_mcmc=10000, n_burn=n_iter_mcmc/2, n_thin=3,
                          
                          out_file=NULL, log_file=NULL,
                          quiet_mcmc=FALSE,
                          parallel_mcmc=FALSE ) {
  
  time_net <- time_all[time_all_idx_net]
  
  V_net <- length(node_all)
  T_net <- length(time_all)
  K_net <- length(layer_all)
  
  
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
  mu_tk <- matrix( #data=0,
    data=runif(T_net*K_net),
    nrow=T_net,
    ncol=K_net )
  mu_tk_mcmc <- array( NA, dim=c(T_net,K_net,n_iter_mcmc_out) )
  
  # Covariance matrix prior for baseline mu_tk
  mu_t_cov_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_mu){ exp(-k*(x-y)^2) } )
  diag(mu_t_cov_prior) <- diag(mu_t_cov_prior) + 1e-3 # numerical stability
  mu_t_cov_prior_inv <- solve(mu_t_cov_prior)
  
  # Covariance matrix prior for parameters beta
  beta_t_cov_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_p){ exp(-k*(x-y)^2) } )
  diag(beta_t_cov_prior) <- diag(beta_t_cov_prior) + 1e-3 # numerical stability
  beta_t_cov_prior_inv <- solve(beta_t_cov_prior)
  
  # Latent coordinates #
  # shared: hth coordinate of actor v at time t shared across the different layers
  x_ith_shared <- array( #data=0,
    data=runif(V_net*T_net*H_dim,-1,1),
    dim=c(V_net,T_net,H_dim) )
  x_ith_shared_mcmc <- array(NA,c(V_net,T_net,H_dim,n_iter_mcmc_out))
  
  # One latent space for each direction #
  x_ith_shared <- list( send=x_ith_shared,
                        receive=x_ith_shared )
  x_ith_shared_mcmc <- list( send=x_ith_shared_mcmc,
                             receive=x_ith_shared_mcmc )
  
  # by layer: hth coordinate of actor v at time t specific to layer k
  if( K_net>1 ){
    x_ithk <- array( #data=0,
      data=runif(V_net*T_net*R_dim*K_net),
      dim=c(V_net,T_net,R_dim,K_net) )
    x_ithk_mcmc <- array(NA,c(V_net,T_net,R_dim,K_net,n_iter_mcmc_out))
    
    # One latent space for each direction #
    x_ithk <- list( send=x_ithk,
                    receive=x_ithk )
    x_ithk_mcmc <- list( send=x_ithk_mcmc,
                         receive=x_ithk_mcmc )
  } else {
    x_ithk <- NULL
    x_ithk_mcmc <- NULL
  }
  
  
  # Covariance matrix prior for coordinates x_t
  x_t_sigma_prior <- outer( time_all, time_all, FUN=function(x,y,k=k_x){ exp(-k*(x-y)^2) } )
  diag(x_t_sigma_prior) <- diag(x_t_sigma_prior) + 1e-3 # numerical stability
  x_t_sigma_prior_inv <- solve(x_t_sigma_prior)
  
  # Predictor coefficients
  beta_z_layer <- beta_z_layer_mcmc <- NULL
  beta_z_edge <- beta_z_edge_mcmc <- NULL
  if(!is.null(pred_all)) {
    if(!is.null(pred_id_layer)){
      beta_z_layer <- matrix(0,nrow=T_net,ncol=nrow(pred_id_layer))
      beta_z_layer_mcmc <- array(NA, dim=c(T_net,nrow(pred_id_layer),n_iter_mcmc_out) )
    }
    if(!is.null(pred_id_edge)){
      beta_z_edge <- matrix(0,nrow=T_net,ncol=nrow(pred_id_edge))
      beta_z_edge_mcmc <- array(NA, dim=c(T_net,nrow(pred_id_edge),n_iter_mcmc_out) )
    }
  }
  
  # Linear predictor for the probability of an edge between actors i and j at time t in layer k
  s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                x_ith_shared=x_ith_shared, x_ithk=x_ithk,
                                pred_all=pred_all, layer_all=layer_all,
                                z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                directed=TRUE )
  
  # Probability of an edge between actors i and j at time t in layer k
  pi_ijtk_mcmc <- array(NA, dim=c(V_net,V_net,T_net,K_net,n_iter_mcmc_out))
  
  #pi_ijt[,,1]
  # all( abs(qlogis( pi_ijt ) - s_ijt)<1e-6,na.rm=T ) # TRUE
  
  # Shrinkage Parameters
  v_shrink_shared <- matrix(NA, nrow=H_dim, ncol=1 )
  v_shrink_shared[1,1] <- rgamma(n=1,shape=a_1,rate=1); v_shrink_shared[-1,1] <- rgamma(n=H_dim-1,shape=a_2,rate=1)
  tau_h_shared <- matrix(cumprod(v_shrink_shared), nrow=H_dim, ncol=1 )
  tau_h_shared_mcmc <- matrix(NA, nrow=H_dim, ncol=n_iter_mcmc_out )
  
  # One latent space for each direction #
  v_shrink_shared <- list( send=v_shrink_shared,
                           receive=v_shrink_shared )
  tau_h_shared <- list( send=tau_h_shared,
                        receive=tau_h_shared )
  tau_h_shared_mcmc <- list( send=tau_h_shared_mcmc,
                             receive=tau_h_shared_mcmc )
  
  # 1/tau_h
  if( K_net>1 ){
    v_shrink_k <- matrix(NA, nrow=R_dim, ncol=K_net )
    v_shrink_k[1,] <- rgamma(n=K_net,shape=a_1,rate=1); v_shrink_k[-1,] <- rgamma(n=K_net*(R_dim-1),shape=a_2,rate=1)
    tau_h_k <- matrix(apply(v_shrink_k,2,cumprod), nrow=R_dim, ncol=K_net )
    tau_h_k_mcmc <- array( NA, dim=c(R_dim,K_net,n_iter_mcmc_out) )
    # One latent space for each direction #
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
    
    # Check singularity of implicit logistic regression for mu and x_ith, in the first iteration
    check_Y<-FALSE
    if(iter_i==2) {check_Y<-FALSE}
    
    #cat(iter_i,",")
    
    
    
    ### Step 1. Update each augmented data w_ijtk from the full conditional Polya-gamma posterior ###
    w_ijtk <- sample_w_ijtk_DynMultiNet_bin( w_ijtk=w_ijtk,
                                             s_ijtk=s_ijtk,
                                             directed=TRUE )
    
    
    
    ### Step 2_mu. Sample mu_tk from its conditional N-variate Gaussian posterior ###
    mu_tk <- sample_mu_tk_DynMultiNet_bin( mu_tk=mu_tk,
                                           y_ijtk=y_ijtk, w_ijtk=w_ijtk, s_ijtk=s_ijtk,
                                           mu_t_cov_prior_inv=mu_t_cov_prior_inv,
                                           directed=TRUE,
                                           use_cpp=TRUE,
                                           parallel_mcmc=parallel_mcmc,
                                           check_Y=check_Y )
    # MCMC chain #
    if(is.element(iter_i,iter_out_mcmc)){
      mu_tk_mcmc[,,match(iter_i,iter_out_mcmc)] <- mu_tk
    }
    
    # update linear predictor
    s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                  x_ith_shared=x_ith_shared, x_ithk=x_ithk,
                                  pred_all=pred_all, layer_all=layer_all,
                                  z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                  beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                  pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                  directed=TRUE )
    
    
    
    ### Step 2_beta. Sample beta_z_layer and beta_z_edge from its conditional N-variate Gaussian posterior ###
    if(!is.null(pred_all)){
      # Layer specific
      if(!is.null(beta_z_layer)&!is.null(pred_id_layer)){
        beta_z_layer <- sample_beta_z_layer_DynMultiNet_bin( beta_z_layer,
                                                             z_tkp, pred_id_layer, pred_all, layer_all,
                                                             y_ijtk, w_ijtk, s_ijtk,
                                                             beta_t_cov_prior_inv,
                                                             use_cpp=TRUE,
                                                             check_Y=check_Y )
        if(is.element(iter_i,iter_out_mcmc)){
          beta_z_layer_mcmc[,,match(iter_i,iter_out_mcmc)] <- beta_z_layer
        }
        
        # update linear predictor
        s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                      x_ith_shared=x_ith_shared, x_ithk=x_ithk,
                                      pred_all=pred_all, layer_all=layer_all,
                                      z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                      beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                      pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                      directed=TRUE )
      }
      # Edge specific
      if(!is.null(beta_z_edge)&!is.null(pred_id_edge)){
        beta_z_edge <- sample_beta_z_edge_DynMultiNet_bin( beta_z_edge,
                                                           z_ijtkp, pred_id_edge, pred_all, layer_all,
                                                           y_ijtk, w_ijtk, s_ijtk,
                                                           beta_t_cov_prior_inv )
        if(is.element(iter_i,iter_out_mcmc)){
          beta_z_edge_mcmc[,,match(iter_i,iter_out_mcmc)] <- beta_z_edge
        }
        # update linear predictor
        s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                      x_ith_shared=x_ith_shared, x_ithk=x_ithk,
                                      pred_all=pred_all, layer_all=layer_all,
                                      z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                      beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                      pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                      directed=TRUE )
      }
    }
    
    ### Step 3. For each unit, block-sample the set of time-varying latent coordinates x_ith ###
    ### SHARED Latent Coordinates ###
    x_ith_shared <- sample_x_ith_shared_DynMultiNet_bin( x_ith_shared=x_ith_shared,
                                                         x_t_sigma_prior_inv=x_t_sigma_prior_inv,
                                                         tau_h=tau_h_shared,
                                                         y_ijtk=y_ijtk, w_ijtk=w_ijtk, s_ijtk=s_ijtk,
                                                         directed=TRUE,
                                                         check_Y=check_Y )
    # MCMC chain #
    if(is.element(iter_i,iter_out_mcmc)){
      x_ith_shared_mcmc[[1]][,,,match(iter_i,iter_out_mcmc)] <- x_ith_shared[[1]]
      x_ith_shared_mcmc[[2]][,,,match(iter_i,iter_out_mcmc)] <- x_ith_shared[[2]]
    }
    
    # update linear predictor
    s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                  x_ith_shared=x_ith_shared, x_ithk=x_ithk,
                                  pred_all=pred_all, layer_all=layer_all,
                                  z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                  beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                  pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                  directed=TRUE )
    
    
    ### LAYER SPECIFIC Latent Coordinates ###
    if( K_net>1 ) {
      ### Step 3A. For each unit, block-sample the EDGE SPECIFIC set of time-varying latent coordinates x_ithk ###
      x_ithk <- sample_x_ithk_DynMultiNet_bin( x_ithk=x_ithk,
                                               x_t_sigma_prior_inv=x_t_sigma_prior_inv,
                                               tau_h=tau_h_k,
                                               y_ijtk=y_ijtk, w_ijtk=w_ijtk, s_ijtk=s_ijtk,
                                               directed=TRUE,
                                               parallel_mcmc=parallel_mcmc,
                                               check_Y=check_Y )
      
      # MCMC chain for x_ithk #
      if(is.element(iter_i,iter_out_mcmc)){
        x_ithk_mcmc[[1]][,,,,match(iter_i,iter_out_mcmc)] <- x_ithk[[1]]
        x_ithk_mcmc[[2]][,,,,match(iter_i,iter_out_mcmc)] <- x_ithk[[2]]
      }
      
      # update linear predictor
      s_ijtk <- get_linpred_s_ijtk( y_ijtk=y_ijtk, mu_tk=mu_tk,
                                    x_ith_shared=x_ith_shared, x_ithk=x_ithk,
                                    pred_all=pred_all, layer_all=layer_all,
                                    z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                                    beta_z_layer=beta_z_layer, beta_z_edge=beta_z_edge,
                                    pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                                    directed=TRUE )
    }
    
    
    
    # Edge probabilities #
    if(is.element(iter_i,iter_out_mcmc)){
      pi_ijtk_mcmc[,,,,match(iter_i,iter_out_mcmc)] <- plogis(s_ijtk)
    }
    
    
    
    ### Step 4. Sample the global shrinkage hyperparameters from conditional gamma distributions ###
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
      if(!is.null(out_file)){
        dmn_mcmc <- list( y_ijtk=y_ijtk,
                          
                          directed=TRUE,
                          weighted=FALSE,
                          
                          n_chains_mcmc=NULL,
                          n_iter_mcmc=n_iter_mcmc, n_burn=n_burn, n_thin=n_thin,
                          time_all_idx_net=time_all_idx_net,
                          
                          a_1=a_1, a_2=a_2,
                          k_x=k_x, k_mu=k_mu, k_p=k_p,
                          
                          node_all=node_all, time_all=time_all, layer_all=layer_all,
                          
                          # For link probabilities #
                          pi_ijtk_mcmc=pi_ijtk_mcmc,
                          mu_tk_mcmc=mu_tk_mcmc,
                          x_ith_shared_mcmc=x_ith_shared_mcmc,
                          x_ithk_mcmc=x_ithk_mcmc,
                          tau_h_shared_mcmc=tau_h_shared_mcmc,
                          tau_h_k_mcmc=tau_h_k_mcmc,
                          
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
        save(dmn_mcmc,file=out_file)
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
                    
                    n_chains_mcmc=NULL,
                    n_iter_mcmc=n_iter_mcmc, n_burn=n_burn, n_thin=n_thin,
                    time_all_idx_net=time_all_idx_net,
                    
                    a_1=a_1, a_2=a_2,
                    k_x=k_x, k_mu=k_mu, k_p=k_p,
                    
                    node_all=node_all, time_all=time_all, layer_all=layer_all,
                    
                    # For link probabilities #
                    pi_ijtk_mcmc=pi_ijtk_mcmc,
                    mu_tk_mcmc=mu_tk_mcmc,
                    x_ith_shared_mcmc=x_ith_shared_mcmc,
                    x_ithk_mcmc=x_ithk_mcmc,
                    tau_h_shared_mcmc=tau_h_shared_mcmc,
                    tau_h_k_mcmc=tau_h_k_mcmc,
                    
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