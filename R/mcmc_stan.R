
#' @import foreach
#' @importFrom rstan stan
#' @keywords internal
mcmc_stan <- function( y_ijtk,
                       node_all, time_all, layer_all,
                       
                       pred_all,
                       pred_id_layer, pred_id_edge,
                       z_tkp, z_ijtkp,
                       
                       H_dim=10, R_dim=10,
                       k_x=0.10, k_mu=0.10, k_p=0.10,
                       a_1=2, a_2=2.5,
                       
                       directed=directed,
                       weighted=weighted,
                       
                       n_iter_mcmc=10000, n_burn=n_iter_mcmc/2, n_thin=1,
                       n_chains_mcmc=4,
                       
                       out_file=NULL,
                       quiet_mcmc=FALSE ) {
  
  V_net <- length(node_all)
  T_net <- length(time_all)
  K_net <- length(layer_all)
  
  y_ijtk[is.na(y_ijtk)]<-0 # Stan does not support NA (in y_ijtk) in data
  stan_data_input <- list( V_net=V_net,T_net=T_net,K_net=K_net,
                           H_dim=H_dim, R_dim=R_dim,
                           y_ijtk=y_ijtk,
                           mu_t_cov_prior = outer( time_all, time_all, FUN=function(x,y,k=k_mu){ exp(-k*(x-y)^2) } ),
                           x_t_cov_prior = outer( time_all, time_all, FUN=function(x,y,k=k_x){ exp(-k*(x-y)^2) } ),
                           a_1=a_1, a_2=a_2,
                           k_mu=k_mu, k_x=k_x )
  
  if( K_net>=1 & is.null(pred_all) & !directed & !weighted ) {
    stan_fit <- rstan::sampling( stanmodels$net_m_1_p_0_d_0_w_0,
                                 data = stan_data_input, 
                                 iter = n_iter_mcmc,warmup=n_burn,thin=n_thin,
                                 chains = n_chains_mcmc,
                                 verbose=!quiet_mcmc,
                                 pars=c("pi_ij_tk",
                                        "mu_tk",
                                        "x_ti_h_shared","x_ti_hk",
                                        "tau_h_shared","tau_hk") )
    dmn_mcmc <- dmn_mcmc_from_stan( stan_fit,
                                    stan_data_input=stan_data_input,
                                    directed=directed,
                                    weighted=weighted )
    
  } else if( K_net>=1 & is.null(pred_all) & directed & !weighted ) {
    stan_fit <- rstan::sampling( stanmodels$net_m_1_p_0_d_1_w_0,
                                 data = stan_data_input, 
                                 iter = n_iter_mcmc,warmup=n_burn,thin=n_thin,
                                 chains = n_chains_mcmc,
                                 verbose=!quiet_mcmc,
                                 pars=c("pi_ij_tk",
                                        "mu_tk",
                                        "x_ti_h_shared","x_ti_hk",
                                        "tau_h_shared","tau_hk") )
    dmn_mcmc <- dmn_mcmc_from_stan( stan_fit,
                                    stan_data_input=stan_data_input,
                                    directed=directed,
                                    weighted=weighted )
  } else {
    stop("Not supported")
  }
  
  if(!is.null(out_file)){
    save( dmn_mcmc , file=out_file )
  }
  
  return( dmn_mcmc )
  
}
