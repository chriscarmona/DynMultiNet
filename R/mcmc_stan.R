
#' @import foreach
#' @importFrom rstan stan
#' @importFrom matrixcalc is.positive.definite
#' @keywords internal
mcmc_stan <- function( y_ijtk,
                       node_all, time_all, layer_all,
                       time_all_idx_net,
                       
                       pred_all,
                       pred_id_layer, pred_id_edge,
                       z_tkp, z_ijtkp,
                       
                       H_dim=10, R_dim=10,
                       
                       k_mu=0.2, lambda_mu=1,
                       k_x=0.2, lambda_x=1,
                       k_p=0.2, lambda_p=1,
                       
                       a_1=2, a_2=2.5,
                       
                       directed=directed,
                       weighted=weighted,
                       
                       n_chains_mcmc=4,
                       n_iter_mcmc=10000, n_burn=n_iter_mcmc/2, n_thin=1,
                       
                       time_fc=NULL,
                       
                       out_file=NULL,
                       quiet_mcmc=FALSE ) {
  
  V_net <- length(node_all)
  T_all <- length(time_all)
  T_net <- length(time_all_idx_net)
  K_net <- length(layer_all)
  
  y_tkij <- aperm(y_ijtk,c(3,4,1,2))
  y_tkij[is.na(y_tkij)]<-0 # Stan does not support NA (in y_ijtk) in data
  
  mu_t_cov_prior = outer( time_all, time_all, FUN=function(x,y,k=k_mu,lambda=lambda_mu){ exp(-k*((x-y)/lambda)^2) } )
  if(!is.positive.definite(mu_t_cov_prior)){ stop('"mu_t_cov_prior" is not positive definite, increase the value of "k_mu"') }
  x_t_cov_prior = outer( time_all, time_all, FUN=function(x,y,k=k_x,lambda=lambda_x){ exp(-k*((x-y)/lambda)^2) } )
  if(!is.positive.definite(x_t_cov_prior)){ stop('"x_t_cov_prior" is not positive definite, increase the value of "k_x"') }
  
  mu_tk_mean <- apply(y_ijtk,c(4),mean,na.rm=T)
  mu_tk_mean <- log(mu_tk_mean/(1-mu_tk_mean))
  mu_tk_mean <- matrix(mu_tk_mean, nrow=T_all, ncol=K_net, byrow=T)
  
  stan_data_input <- list( V_net=V_net,T_net=T_net,K_net=K_net,
                           H_dim=H_dim, R_dim=R_dim,
                           y_tkij=y_tkij,
                           
                           k_mu=k_mu, lambda_mu=lambda_mu,
                           k_x=k_x, lambda_x=lambda_x,
                           k_p=k_p, lambda_p=lambda_p,
                           
                           a_1=a_1, a_2=a_2,
                           
                           T_all=T_all,
                           mu_tk_mean = mu_tk_mean,
                           # mu_t_cov_prior=mu_t_cov_prior,
                           # x_t_cov_prior=x_t_cov_prior,
                           time_all=time_all,
                           time_all_idx_net=time_all_idx_net )
  
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
    save( stan_fit , file=out_file )
    dmn_mcmc <- dmn_mcmc_from_stan( stan_fit=stan_fit,
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
    save( stan_fit , file=out_file )
    dmn_mcmc <- stan_fit
    dmn_mcmc <- dmn_mcmc_from_stan( stan_fit=stan_fit,
                                    directed=directed,
                                    weighted=weighted )
  } else {
    stop("Not supported")
  }
  
  dmn_mcmc$y_ijtk = y_ijtk
  dmn_mcmc$n_chains_mcmc = n_chains_mcmc
  dmn_mcmc$n_iter_mcmc =n_iter_mcmc
  dmn_mcmc$n_burn = n_burn
  dmn_mcmc$n_thin = n_thin
  
  dmn_mcmc$node_all = node_all
  dmn_mcmc$time_all = time_all
  dmn_mcmc$layer_all = layer_all
  
  dmn_mcmc$a_1=a_1; dmn_mcmc$a_2=a_2
  dmn_mcmc$k_mu=k_mu; dmn_mcmc$k_x=k_x
  
  dmn_mcmc$time_all_idx_net=time_all_idx_net
  
  if(!is.null(out_file)){
    save( dmn_mcmc , file=out_file )
  }
  
  return( dmn_mcmc )
  
}
