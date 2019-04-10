
#' @import foreach
#' @importFrom rstan stan sampling
#' @importFrom matrixcalc is.positive.definite
#' @keywords internal
VB_stan <- function( y_ijtk,
                     node_all, time_all, layer_all,
                     time_all_idx_net,
                     
                     pred_all,
                     pred_id_layer, pred_id_edge,
                     z_tkp, z_ijtkp,
                     
                     H_dim=10, R_dim=10,
                     
                     delta_mu=1,
                     delta_x=1,
                     delta_p=1,
                     
                     delta_lambda=1,
                     delta_u=1,
                     delta_q=1,
                     
                     a_1=2, a_2=2.5,
                     
                     directed=directed,
                     weighted=weighted,
                     
                     n_chains_mcmc=4,
                     n_iter_mcmc=10000, n_burn=n_iter_mcmc/2, n_thin=1,
                     
                     time_fc=NULL,
                     
                     out_file="DynMultiNet_mcmc_result",
                     quiet_mcmc=FALSE,
                     only_read_csv_stan_fit=FALSE ) {
  
  V_net <- length(node_all)
  T_all <- length(time_all)
  T_net <- length(time_all_idx_net)
  K_net <- length(layer_all)
  
  y_tkij <- aperm(y_ijtk,c(3,4,1,2))
  y_tkij[is.na(y_tkij)]<-0 # Stan does not support NA (in y_ijtk) in data
  
  y_ijtk_aux <- y_ijtk
  y_ijtk_aux[y_ijtk_aux>0] <- 1
  mu_tk_mean <- apply(y_ijtk_aux,c(4),mean,na.rm=T)
  mu_tk_mean <- log(mu_tk_mean/(1-mu_tk_mean))
  mu_tk_mean <- matrix(mu_tk_mean, nrow=T_all, ncol=K_net, byrow=T)
  rm(y_ijtk_aux)
  
  y_ijtk_aux <- y_ijtk
  y_ijtk_aux[y_ijtk_aux==0] <- NA
  lambda_tk_mean <- apply(y_ijtk_aux,c(4),mean,na.rm=T)
  lambda_tk_mean <- matrix(lambda_tk_mean, nrow=T_all, ncol=K_net, byrow=T)
  rm(y_ijtk_aux)
  
  stan_data_input <- list( V_net=V_net,T_net=T_net,K_net=K_net,
                           H_dim=H_dim, R_dim=R_dim,
                           y_tkij=y_tkij,
                           
                           T_all=T_all,
                           time_all=time_all,
                           time_all_idx_net=time_all_idx_net,
                           
                           a_1=a_1, a_2=a_2,
                           
                           delta_mu=delta_mu,
                           delta_x=delta_x,
                           delta_p=delta_p,
                           mu_tk_mean = mu_tk_mean,
                           
                           delta_lambda=delta_lambda,
                           delta_u=delta_u,
                           delta_q=delta_q,
                           lambda_tk_mean = lambda_tk_mean )
  
  cat("MCMC will be saved in ",n_chains_mcmc," csv files:\n\n",paste(out_file,"_",1:n_chains_mcmc,".csv\n",sep=""),"\n\n",sep="")
  if( K_net>=1 & is.null(pred_all) & !directed & !weighted ) {
    
    stan_fit <- rstan::sampling( stanmodels$net_m_1_p_0_d_0_w_0,
                                 data = stan_data_input, 
                                 iter = n_iter_mcmc,warmup=n_burn,thin=n_thin,
                                 chains = n_chains_mcmc,
                                 verbose=!quiet_mcmc,
                                 sample_file=out_file,
                                 pars=c("pi_ij_tk",
                                        "mu_tk",
                                        "x_ti_h_shared","x_ti_hk",
                                        "tau_h_shared","tau_hk") )
    
    saveRDS( stan_fit , file=paste(out_file,"_stan.rds",sep="") )
    dmn_mcmc <- dmn_mcmc_from_stan( stan_fit=stan_fit,
                                    directed=directed,
                                    weighted=weighted )
    
  } else if( K_net>=1 & is.null(pred_all) & directed & !weighted ) {
    
    stan_fit <- rstan::sampling( stanmodels$net_m_1_p_0_d_1_w_0,
                                 data = stan_data_input, 
                                 iter = n_iter_mcmc,warmup=n_burn,thin=n_thin,
                                 chains = n_chains_mcmc,
                                 verbose=!quiet_mcmc,
                                 sample_file=out_file,
                                 pars=c("pi_ij_tk",
                                        "mu_tk",
                                        "x_ti_h_shared","x_ti_hk",
                                        "tau_h_shared","tau_hk") )
    
    saveRDS( stan_fit , file=paste(out_file,"_stan.rds",sep="") )
    dmn_mcmc <- dmn_mcmc_from_stan( stan_fit=stan_fit,
                                    directed=directed,
                                    weighted=weighted )
  } else if( K_net>=1 & is.null(pred_all) & directed & weighted ) {
    if( only_read_csv_stan_fit ){
      dmn_mcmc <- dmn_mcmc_from_stan_failed( sample_file=out_file,
                                             n_chains_mcmc=n_chains_mcmc,
                                             n_iter_mcmc=n_iter_mcmc,n_burn=n_burn,
                                             V_net=V_net,T_all=T_all,K_net=K_net,
                                             H_dim=H_dim,R_dim=R_dim,
                                             directed=directed,
                                             weighted=weighted )
    } else {
      tryCatch(
        {
          stan_fit <- rstan::sampling( stanmodels$net_m_1_p_0_d_1_w_1,
                                       data = stan_data_input, 
                                       iter = n_iter_mcmc,warmup=n_burn,thin=n_thin,
                                       chains = n_chains_mcmc,
                                       verbose=!quiet_mcmc,
                                       sample_file=out_file,
                                       pars=c("pi_ij_tk",
                                              "mu_t_k",
                                              "x_t_hi_shared","x_t_hki",
                                              "tau_h_shared","tau_hk",
                                              "r_ij_tk",
                                              "sigma_w_k",
                                              "lambda_t_k",
                                              "u_t_hi_shared","u_t_hki",
                                              "rho_h_shared","rho_hk"
                                       ) )
          
          saveRDS( stan_fit , file=paste(out_file,"_stan.rds",sep="") )
          dmn_mcmc <- dmn_mcmc_from_stan( stan_fit=stan_fit,
                                          directed=directed,
                                          weighted=weighted )
        },
        error=function(cond){
          message("Error message:")
          message(cond)
          dmn_mcmc <- dmn_mcmc_from_stan_failed( sample_file=out_file,
                                                 n_chains_mcmc=n_chains_mcmc,
                                                 n_iter_mcmc=n_iter_mcmc,n_burn=n_burn,
                                                 V_net=V_net,T_all=T_all,K_net=K_net,
                                                 H_dim=H_dim,R_dim=R_dim,
                                                 directed=directed,
                                                 weighted=weighted )
          return(dmn_mcmc)
        }
      )
    }
    
  } else {
    stop("Apologies, network not supported by DynMultiNet.")
  }
  
  dmn_mcmc$y_ijtk = y_ijtk
  dimnames(dmn_mcmc$pi_ijtk) = list( node_all,
                                     node_all,
                                     time_all,
                                     layer_all )
  
  if(weighted){
    dimnames(dmn_mcmc$r_ijtk) = dimnames(dmn_mcmc$pi_ijtk)
  }
  
  dmn_mcmc$n_chains_mcmc = n_chains_mcmc
  dmn_mcmc$n_iter_mcmc = n_iter_mcmc
  dmn_mcmc$n_burn = n_burn
  dmn_mcmc$n_thin = n_thin
  
  dmn_mcmc$node_all = node_all
  dmn_mcmc$time_all = time_all
  dmn_mcmc$layer_all = layer_all
  
  dmn_mcmc$a_1=a_1; dmn_mcmc$a_2=a_2
  
  dmn_mcmc$time_all_idx_net=time_all_idx_net
  
  saveRDS( dmn_mcmc , file=paste(out_file,".rds",sep="") )
  
  return( dmn_mcmc )
  
}
