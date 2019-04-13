
#' @import rstan
#' 
#' @keywords internal
#' 

VB_stan <- function( y_ijtk,
                     directed=FALSE, weighted=FALSE,
                     
                     node_all, time_all, layer_all,
                     
                     x_ijtkp=NULL,
                     
                     H_dim=10, R_dim=10,
                     
                     add_eff_weight=FALSE,
                     add_eff_link=FALSE,
                     
                     delta=36,
                     lat_mean=TRUE,
                     sigma_lat_mean=5,
                     
                     vb_algorithm=c("meanfield","fullrank")[2],
                     rds_file=NULL, log_file=NULL ) {
  
  V_net <- length(node_all)
  T_net <- length(time_all)
  K_net <- length(layer_all)
  
  y_tkij <- aperm(y_ijtk,c(3,4,1,2),drop=F)
  y_tkij[is.na(y_tkij)]<-1234567890 # Stan does not support NA (in y_ijtk) in data
  
  stan_data_input <- list( V_net=V_net, T_net=T_net, K_net=K_net,
                           H_dim=H_dim, R_dim=R_dim,
                           
                           y_tkij=y_tkij,
                           
                           T_net=T_net,
                           time_all=time_all,
                           
                           delta=delta )
  
  if( K_net==1 & is.null(x_ijtkp) & directed & weighted ) {
    
    stan_fit <- rstan::vb( stanmodels$net_m_0_p_0_d_1_w_1,
                           data = stan_data_input,
                           sample_file=out_file,
                           init=0,
                           pars=c("pi_ij_tk",
                                  
                                  "eta_t_k",
                                  "ab_t_hi_shared","ab_t_hki",
                                  "eta_bar_k",
                                  "ab_bar_hi_shared", "ab_bar_hki",
                                  
                                  "mu_ij_tk",
                                  
                                  "theta_t_k",
                                  "theta_bar_k",
                                  "uv_t_hi_shared","uv_t_hki",
                                  "uv_bar_hi_shared", "uv_bar_hki",
                                  
                                  "sigma_w_k"),
                           algorithm=vb_algorithm )
    saveRDS( stan_fit , file=paste(out_file,"_stan.rds",sep="") )
    return(stan_fit)
    dmn_mcmc <- dmn_mcmc_from_stan( stan_fit=stan_fit,
                                    directed=directed,
                                    weighted=weighted )
    
  } else if ( K_net>1 & is.null(x_ijtkp) & directed & weighted ) {
    
    stan_fit <- rstan::vb( stanmodels$net_m_1_p_0_d_1_w_1,
                           data = stan_data_input,
                           sample_file=out_file,
                           init=0,
                           pars=c("pi_ij_tk",
                                  
                                  "eta_t_k",
                                  "ab_t_hi_shared","ab_t_hki",
                                  "eta_bar_k",
                                  "ab_bar_hi_shared", "ab_bar_hki",
                                  
                                  "mu_ij_tk",
                                  
                                  "theta_t_k",
                                  "theta_bar_k",
                                  "uv_t_hi_shared","uv_t_hki",
                                  "uv_bar_hi_shared", "uv_bar_hki",
                                  
                                  "sigma_w_k"),
                           algorithm=vb_algorithm )
    saveRDS( stan_fit , file=paste(out_file,"_stan.rds",sep="") )
    return(stan_fit)
    dmn_mcmc <- dmn_mcmc_from_stan( stan_fit=stan_fit,
                                    directed=directed,
                                    weighted=weighted )
    
  } else {
    stop("Apologies, variational inference not supported by DynMultiNet at the moment")
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
  
  saveRDS( dmn_mcmc , file=paste(out_file,".rds",sep="") )
  
  return( dmn_mcmc )
  
}
