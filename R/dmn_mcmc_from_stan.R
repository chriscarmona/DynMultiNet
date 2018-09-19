
dmn_mcmc_from_stan <- function( stan_fit,
                                directed=FALSE,
                                weighted=FALSE ) {
  
  posterior <- rstan::extract(stan_fit)
  
  if( !directed & !weighted ) {
    dmn_mcmc <- list( y_ijtk=NULL,
                      
                      directed = directed,
                      weighted = weighted,
                      
                      n_chains_mcmc=NULL,
                      n_iter_mcmc = NULL, n_burn = NULL, n_thin = NULL,
                      
                      a_1=NULL, a_2=NULL,
                      k_mu=NULL, k_x=NULL,
                      
                      lp=posterior[["lp__"]],
                      
                      pi_ijtk_mcmc = aperm(posterior[["pi_ij_tk"]],c(4,5,2,3,1)),
                      
                      node_all=NULL, time_all=NULL, layer_all=NULL,
                      time_all_idx_net=NULL,
                      
                      mu_tk_mcmc = posterior[["mu_tk"]],
                      
                      x_ith_shared_mcmc = aperm(posterior[["x_ti_h_shared"]],c(4,3,2,1)),
                      x_ithk_mcmc = aperm(posterior[["x_ti_hk"]],c(5,4,2,3,1)),
                      tau_h_shared_mcmc = aperm(posterior[["tau_h_shared"]],c(2,1)),
                      tau_h_k_mcmc = aperm(posterior[["tau_hk"]],c(2,3,1)),
                      
                      pred_id_layer=NULL, pred_id_edge=NULL,
                      beta_z_layer_mcmc=NULL,
                      beta_z_edge_mcmc=NULL )
  } else if( directed & !weighted ) {
    dmn_mcmc <- list( y_ijtk=NULL,
                      
                      directed = directed,
                      weighted = weighted,
                      
                      n_chains_mcmc=NULL,
                      n_iter_mcmc = NULL, n_burn = NULL, n_thin = NULL,
                      time_all_idx_net=NULL,
                      
                      a_1=NULL, a_2=NULL,
                      k_mu=NULL, k_x=NULL,
                      
                      pi_ijtk_mcmc = aperm(posterior[["pi_ij_tk"]],c(4,5,2,3,1)),
                      
                      node_all=NULL, time_all=NULL, layer_all=NULL,
                      
                      mu_tk_mcmc = aperm(posterior[["mu_tk"]],c(2,3,1)),
                      
                      x_ith_shared_mcmc = list( aperm(posterior[["x_ti_h_shared"]][,1,,,],c(4,3,2,1)),
                                                aperm(posterior[["x_ti_h_shared"]][,2,,,],c(4,3,2,1)) ),
                      x_ithk_mcmc = list( aperm(posterior[["x_ti_hk"]][,1,,,,],c(5,4,2,3,1)),
                                          aperm(posterior[["x_ti_hk"]][,2,,,,],c(5,4,2,3,1)) ),
                      tau_h_shared_mcmc = list( aperm(posterior[["tau_h_shared"]][,1,],c(2,1) ),
                                                aperm(posterior[["tau_h_shared"]][,1,],c(2,1)) ),
                      tau_h_k_mcmc = list( aperm(posterior[["tau_hk"]][,1,,],c(2,3,1)),
                                           aperm(posterior[["tau_hk"]][,1,,],c(2,3,1)) ),
                      
                      pred_id_layer=NULL, pred_id_edge=NULL,
                      beta_z_layer_mcmc=NULL,
                      beta_z_edge_mcmc=NULL )
  } else if( directed & weighted ) {
    dmn_mcmc <- list( y_ijtk=NULL,
                      
                      directed = directed,
                      weighted = weighted,
                      
                      n_chains_mcmc=NULL,
                      n_iter_mcmc = NULL, n_burn = NULL, n_thin = NULL,
                      time_all_idx_net=NULL,
                      
                      a_1=NULL, a_2=NULL,
                      k_mu=NULL, k_x=NULL,
                      
                      node_all=NULL, time_all=NULL, layer_all=NULL,
                      
                      lp=posterior[["lp__"]],
                      
                      # For link probabilities #
                      pi_ijtk_mcmc = aperm(posterior[["pi_ij_tk"]],c(4,5,2,3,1)),
                      
                      mu_tk_mcmc = aperm(posterior[["mu_kt"]],c(3,2,1)),
                      x_ith_shared_mcmc = list( aperm(posterior[["x_it_h_shared"]][,1,,,],c(3,4,2,1)),
                                                aperm(posterior[["x_it_h_shared"]][,2,,,],c(3,4,2,1)) ),
                      x_ithk_mcmc = list( aperm(posterior[["x_it_hk"]][,1,,,,],c(4,5,2,3,1)),
                                          aperm(posterior[["x_it_hk"]][,2,,,,],c(4,5,2,3,1)) ),
                      tau_h_shared_mcmc = list( aperm(posterior[["tau_h_shared"]][,1,],c(2,1) ),
                                                aperm(posterior[["tau_h_shared"]][,1,],c(2,1)) ),
                      tau_h_k_mcmc = list( aperm(posterior[["tau_hk"]][,1,,],c(2,3,1)),
                                           aperm(posterior[["tau_hk"]][,1,,],c(2,3,1)) ),
                      
                      # For link weights #
                      r_ijtk_mcmc = aperm(posterior[["r_ij_tk"]],c(4,5,2,3,1)),
                      sigma_w_k_mcmc = aperm(posterior[["sigma_w_k"]],c(2,1)),
                      
                      lambda_tk_mcmc = aperm(posterior[["lambda_kt"]],c(3,2,1)),
                      u_ith_shared_mcmc = list( aperm(posterior[["u_it_h_shared"]][,1,,,],c(3,4,2,1)),
                                                aperm(posterior[["u_it_h_shared"]][,2,,,],c(3,4,2,1)) ),
                      u_ithk_mcmc = list( aperm(posterior[["u_it_hk"]][,1,,,,],c(4,5,2,3,1)),
                                          aperm(posterior[["u_it_hk"]][,2,,,,],c(4,5,2,3,1)) ),
                      rho_h_shared_mcmc = list( aperm(posterior[["rho_h_shared"]][,1,],c(2,1) ),
                                                aperm(posterior[["rho_h_shared"]][,1,],c(2,1)) ),
                      rho_h_k_mcmc = list( aperm(posterior[["rho_hk"]][,1,,],c(2,3,1)),
                                           aperm(posterior[["rho_hk"]][,1,,],c(2,3,1)) ),
                      
                      pred_id_layer=NULL, pred_id_edge=NULL,
                      beta_z_layer_mcmc=NULL,
                      beta_z_edge_mcmc=NULL )
  } else {
    stop('"dmn_mcmc_from_stan" does not supported your specification.')
  }
  dmn_mcmc <- structure( dmn_mcmc, class="dmn_mcmc" )
  
  return( dmn_mcmc )
}

dmn_mcmc_from_stan_failed <- function( sample_file,
                                       n_chains_mcmc,
                                       n_iter_mcmc,n_burn,
                                       V_net,T_all,K_net,H_dim,R_dim,
                                       directed=FALSE,
                                       weighted=FALSE ) {
  
  if( directed & weighted ) {
    
    pars=c( "lp__","pi_ij_tk",
            "mu_tk",
            "x_it_h_shared","x_it_hk",
            "tau_h_shared","tau_hk",
            "r_ij_tk",
            "sigma_w_k",
            "lambda_tk",
            "u_it_h_shared","u_it_hk",
            "rho_h_shared","rho_hk" )
    
    
    chain_files <- paste(sample_file,"_",1:n_chains_mcmc,".csv",sep="")
    posterior <- list(NULL)
    
    cat("\nLoading chains from output csv files...\n\n")
    
    for(file_i in 1:n_chains_mcmc) { # file_i<-1
      cat("chain =",file_i,"\n")
      # finding starting position of chains
      tmp = readLines( chain_files[file_i] )
      skip.rows <- which(substr(tmp,1,1)=="#")
      tmp = tmp[-(skip.rows)]
      tmpFile = tempfile()
      on.exit(unlink(tmpFile))
      writeLines(tmp,tmpFile)
    
      chain_stan_names <- data.table::fread( tmpFile,
                                             nrows=0,
                                             data.table=FALSE )
      chain_stan_names <- colnames(chain_stan_names)
      col_list <- list(NULL)
      cat("par = ")
      for(par_i in 1:length(pars)) { # par_i<-3
        
        col_list[[par_i]] <- which(substr(chain_stan_names,1,nchar(pars[par_i]))==pars[par_i])
        chain_aux <- data.table::fread( tmpFile,
                                        select=col_list[[par_i]],
                                        data.table=FALSE )
        if(par_i==1){
          if( nrow(chain_aux)<n_iter_mcmc ) {
            warning('The number of mcmc iterations is less than n_iter_mcmc=',n_iter_mcmc)
          }
          if( (nrow(chain_aux)>n_burn)&(n_burn>0) ) {
            cat("Eliminating n_burn=",n_burn,"iterations\n")
            chain_aux <- chain_aux[-(1:n_burn),]
          } else {
            warning('The number of mcmc iterations is less than n_burn=',n_burn,", no iterations will be discarded!")
          }
        }
        if(file_i==1){
          posterior[[par_i]] <- as.matrix(chain_aux)
        } else {
          posterior[[par_i]] <- rbind(posterior[[par_i]],as.matrix(chain_aux))
        }
        cat(pars[par_i],", ")
      }
      cat("\n")
    }
    names(posterior) <- pars
    
    n_iter_mcmc_eff <- nrow(posterior[["mu_tk"]])
    
    ### For link probabilities ###
    
    if( all(dim(posterior[["pi_ij_tk"]])==c(n_iter_mcmc_eff,T_all*K_net*V_net*V_net)) ){
      posterior[["pi_ij_tk"]] <- array(data=c(posterior[["pi_ij_tk"]]),dim=c(n_iter_mcmc_eff,T_all,K_net,V_net,V_net))
    } else { stop('Error loading samples from "pi_ij_tk"') }
    
    if( all(dim(posterior[["mu_kt"]])==c(n_iter_mcmc_eff,T_all*K_net)) ){
      posterior[["mu_kt"]] <- array(data=c(posterior[["mu_kt"]]),dim=c(n_iter_mcmc_eff,K_net,T_all))
    } else { stop('Error loading samples from "mu_kt"') }
    
    if( all(dim(posterior[["x_it_h_shared"]])==c(n_iter_mcmc_eff,2*H_dim*T_all*V_net)) ){
      posterior[["x_it_h_shared"]] <- array(data=c(posterior[["x_it_h_shared"]]),dim=c(n_iter_mcmc_eff,2,H_dim,V_net,T_all))
    } else { stop('Error loading samples from "x_it_h_shared"') }
    
    if( all(dim(posterior[["x_it_hk"]])==c(n_iter_mcmc_eff,2*R_dim*K_net*T_all*V_net)) ){
      posterior[["x_it_hk"]] <- array(data=c(posterior[["x_it_hk"]]),dim=c(n_iter_mcmc_eff,2,R_dim,K_net,V_net,T_all))
    } else { stop('Error loading samples from "x_it_hk"') }
    
    if( all(dim(posterior[["tau_h_shared"]])==c(n_iter_mcmc_eff,2*H_dim)) ){
      posterior[["tau_h_shared"]] <- array(data=c(posterior[["tau_h_shared"]]),dim=c(n_iter_mcmc_eff,2,H_dim))
    } else { stop('Error loading samples from "tau_h_shared"') }
    
    if( all(dim(posterior[["tau_hk"]])==c(n_iter_mcmc_eff,2*R_dim*K_net)) ){
      posterior[["tau_hk"]] <- array(data=c(posterior[["tau_hk"]]),dim=c(n_iter_mcmc_eff,2,R_dim,K_net))
    } else { stop('Error loading samples from "tau_hk"') }
    
    ### For link weights ###
    
    if( all(dim(posterior[["r_ij_tk"]])==c(n_iter_mcmc_eff,T_all*K_net*V_net*V_net)) ){
      posterior[["r_ij_tk"]] <- array(data=c(posterior[["r_ij_tk"]]),dim=c(n_iter_mcmc_eff,T_all,K_net,V_net,V_net))
    } else { stop('Error loading samples from "r_ij_tk"') }
    
    if( all(dim(posterior[["sigma_w_k"]])==c(n_iter_mcmc_eff,K_net)) ){
      posterior[["sigma_w_k"]] <- array(data=c(posterior[["sigma_w_k"]]),dim=c(n_iter_mcmc_eff,K_net))
    } else { stop('Error loading samples from "sigma_w_k"') }
    
    if( all(dim(posterior[["lambda_kt"]])==c(n_iter_mcmc_eff,T_all*K_net)) ){
      posterior[["lambda_kt"]] <- array(data=c(posterior[["lambda_kt"]]),dim=c(n_iter_mcmc_eff,K_net,T_all))
    } else { stop('Error loading samples from "lambda_kt"') }
    
    if( all(dim(posterior[["u_it_h_shared"]])==c(n_iter_mcmc_eff,2*H_dim*T_all*V_net)) ){
      posterior[["u_it_h_shared"]] <- array(data=c(posterior[["u_it_h_shared"]]),dim=c(n_iter_mcmc_eff,2,H_dim,V_net,T_all))
    } else { stop('Error loading samples from "u_it_h_shared"') }
    
    if( all(dim(posterior[["u_it_hk"]])==c(n_iter_mcmc_eff,2*R_dim*K_net*T_all*V_net)) ){
      posterior[["u_it_hk"]] <- array(data=c(posterior[["u_it_hk"]]),dim=c(n_iter_mcmc_eff,2,R_dim,K_net,V_net,T_all))
    } else { stop('Error loading samples from "u_it_hk"') }
    
    if( all(dim(posterior[["rho_h_shared"]])==c(n_iter_mcmc_eff,2*H_dim)) ){
      posterior[["rho_h_shared"]] <- array(data=c(posterior[["rho_h_shared"]]),dim=c(n_iter_mcmc_eff,2,H_dim))
    } else { stop('Error loading samples from "rho_h_shared"') }
    
    if( all(dim(posterior[["rho_hk"]])==c(n_iter_mcmc_eff,2*R_dim*K_net)) ){
      posterior[["rho_hk"]] <- array(data=c(posterior[["rho_hk"]]),dim=c(n_iter_mcmc_eff,2,R_dim,K_net))
    } else { stop('Error loading samples from "rho_hk"') }
    
    dmn_mcmc <- list( y_ijtk=NULL,
                      
                      directed = directed,
                      weighted = weighted,
                      
                      n_chains_mcmc=NULL,
                      n_iter_mcmc = NULL, n_burn = NULL, n_thin = NULL,
                      time_all_idx_net=NULL,
                      
                      a_1=NULL, a_2=NULL,
                      k_mu=NULL, k_x=NULL,
                      
                      node_all=NULL, time_all=NULL, layer_all=NULL,
                      
                      # For link probabilities #
                      pi_ijtk_mcmc = aperm(posterior[["pi_ij_tk"]],c(4,5,2,3,1)),
                      
                      mu_tk_mcmc = aperm(posterior[["mu_kt"]],c(3,2,1)),
                      x_ith_shared_mcmc = list( aperm(posterior[["x_it_h_shared"]][,1,,,],c(3,4,2,1)),
                                                aperm(posterior[["x_it_h_shared"]][,2,,,],c(3,4,2,1)) ),
                      x_ithk_mcmc = list( aperm(posterior[["x_it_hk"]][,1,,,,],c(4,5,2,3,1)),
                                          aperm(posterior[["x_it_hk"]][,2,,,,],c(4,5,2,3,1)) ),
                      tau_h_shared_mcmc = list( aperm(posterior[["tau_h_shared"]][,1,],c(2,1) ),
                                                aperm(posterior[["tau_h_shared"]][,1,],c(2,1)) ),
                      tau_h_k_mcmc = list( aperm(posterior[["tau_hk"]][,1,,],c(2,3,1)),
                                           aperm(posterior[["tau_hk"]][,1,,],c(2,3,1)) ),
                      
                      # For link weights #
                      r_ijtk_mcmc = aperm(posterior[["r_ij_tk"]],c(4,5,2,3,1)),
                      sigma_w_k_mcmc = aperm(posterior[["sigma_w_k"]],c(2,1)),
                      
                      lambda_tk_mcmc = aperm(posterior[["lambda_kt"]],c(3,2,1)),
                      u_ith_shared_mcmc = list( aperm(posterior[["u_it_h_shared"]][,1,,,],c(3,4,2,1)),
                                                aperm(posterior[["u_it_h_shared"]][,2,,,],c(3,4,2,1)) ),
                      u_ithk_mcmc = list( aperm(posterior[["u_it_hk"]][,1,,,,],c(4,5,2,3,1)),
                                          aperm(posterior[["u_it_hk"]][,2,,,,],c(4,5,2,3,1)) ),
                      rho_h_shared_mcmc = list( aperm(posterior[["rho_h_shared"]][,1,],c(2,1) ),
                                                aperm(posterior[["rho_h_shared"]][,1,],c(2,1)) ),
                      rho_h_k_mcmc = list( aperm(posterior[["rho_hk"]][,1,,],c(2,3,1)),
                                           aperm(posterior[["rho_hk"]][,1,,],c(2,3,1)) ),
                      
                      pred_id_layer=NULL, pred_id_edge=NULL,
                      beta_z_layer_mcmc=NULL,
                      beta_z_edge_mcmc=NULL )
  } else {
    stop('"dmn_mcmc_from_stan_failed" does not support your specification.')
  }
  dmn_mcmc <- structure( dmn_mcmc, class="dmn_mcmc" )
  
  saveRDS( dmn_mcmc , file=paste(out_file,".rds",sep="") )
  
  return( dmn_mcmc )
  
}
