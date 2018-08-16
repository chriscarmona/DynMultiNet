
#' @import foreach
#' @importFrom rstan stan
#' @keywords internal
mcmc_d_0_w_0_stan <- function( y_ijtk,
                               node_all, time_all, layer_all,
                               
                               pred_all,
                               pred_id_layer, pred_id_edge,
                               z_tkp, z_ijtkp,
                               
                               H_dim=10, R_dim=10,
                               k_x=0.10, k_mu=0.10, k_p=0.10,
                               a_1=2, a_2=2.5,
                               
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
                           # k_x=k_x, k_mu=k_mu,# k_p=k_p,
                           a_1=a_1, a_2=a_2 )
  
  DynMultiNet_stan <- rstan::sampling( stanmodels$net_d_0_w_0,
                                       data = stan_data_input, 
                                       iter = n_iter_mcmc,warmup=n_burn,thin=n_thin,
                                       chains = n_chains_mcmc,
                                       verbose=!quiet_mcmc,
                                       pars=c("pi_ij_tk",
                                              "mu_tk",
                                              "x_ti_h_shared","x_ti_hk",
                                              "tau_h_shared","tau_hk") )
  
  if(!is.null(out_file)){
    save(DynMultiNet_stan,file=out_file)
  }
  # DynMultiNet_mcmc <- structure( DynMultiNet_mcmc, class="DynMultiNet_mcmc" )
  
  return( DynMultiNet_mcmc )
  
}
