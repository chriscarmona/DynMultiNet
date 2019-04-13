
// saved as net_m_1_p_0_d_1_w_1.stan

// Supports Dynamic Networks with the following characteristics:
//   multi-layer: TRUE
//   predictors: FALSE
//   directed: TRUE
//   weigthed: TRUE

data {
  int<lower=1> V_net; // number of nodes
  int<lower=1> T_net; // number of time-steps in the observed network
  int<lower=1> K_net; // number of layers
  int<lower=1> H_dim; // number of dimensions in shared latent space
  int<lower=1> R_dim; // number of dimensions in layer-specific latent space
  
  real<lower=0> y_tkij[T_net,K_net,V_net,V_net];
  
  vector[T_net] time_all; // Time stamps of all time-steps
  
  real<lower=0> delta;
  real<lower=0> sigma_lat_mean;
}

transformed data {
  
  cov_matrix[T_net] eta_t_cov;
  cholesky_factor_cov[T_net] eta_t_cov_sqrt;
  cov_matrix[T_net] ab_t_cov;
  cholesky_factor_cov[T_net] ab_t_cov_sqrt;
  
  cov_matrix[T_net] theta_t_cov;
  cholesky_factor_cov[T_net] theta_t_cov_sqrt;
  cov_matrix[T_net] uv_t_cov;
  cholesky_factor_cov[T_net] uv_t_cov_sqrt;
  
  for(t1 in 1:T_net) {
    for(t2 in 1:t1) {
      // for link probability
      eta_t_cov[t1,t2] = exp(-((time_all[t1]-time_all[t2])/delta)^2);
      eta_t_cov[t2,t1] = exp(-((time_all[t1]-time_all[t2])/delta)^2);
      ab_t_cov[t1,t2] = exp(-((time_all[t1]-time_all[t2])/delta)^2);
      ab_t_cov[t2,t1] = exp(-((time_all[t1]-time_all[t2])/delta)^2);
      // for link weight
      theta_t_cov[t1,t2] = exp(-((time_all[t1]-time_all[t2])/delta)^2);
      theta_t_cov[t2,t1] = exp(-((time_all[t1]-time_all[t2])/delta)^2);
      uv_t_cov[t1,t2] = exp(-((time_all[t1]-time_all[t2])/delta)^2);
      uv_t_cov[t2,t1] = exp(-((time_all[t1]-time_all[t2])/delta)^2);
    }
  }
  // for link probability
  eta_t_cov = eta_t_cov + diag_matrix(rep_vector(1e-4, T_net));
  eta_t_cov_sqrt = cholesky_decompose( eta_t_cov );
  ab_t_cov = ab_t_cov + diag_matrix(rep_vector(1e-4, T_net));
  ab_t_cov_sqrt = cholesky_decompose( ab_t_cov );
  // for link weight
  theta_t_cov = theta_t_cov + diag_matrix(rep_vector(1e-4, T_net));
  theta_t_cov_sqrt = cholesky_decompose( theta_t_cov );
  uv_t_cov = uv_t_cov + diag_matrix(rep_vector(1e-4, T_net));
  uv_t_cov_sqrt = cholesky_decompose( uv_t_cov );
  
}

parameters {
  // for link probability
  vector[T_net] eta_t_k[K_net];
  vector[T_net] ab_t_hi_shared[2,H_dim,V_net];
  vector[T_net] ab_t_hki[2,R_dim,K_net,V_net];
  
  real eta_bar_k[K_net];
  real ab_bar_hi_shared[2,H_dim,V_net];
  real ab_bar_hki[2,R_dim,K_net,V_net];
  
  // for link weight
  vector[T_net] theta_t_k[K_net];
  vector[T_net] uv_t_hi_shared[2,H_dim,V_net];
  vector[T_net] uv_t_hki[2,R_dim,K_net,V_net];
  
  real theta_bar_k[K_net];
  real uv_bar_hi_shared[2,H_dim,V_net];
  real uv_bar_hki[2,R_dim,K_net,V_net];
  
  real log_sigma_w_k[K_net];
}

transformed parameters {
  // Linear predictors
  matrix[V_net,V_net] gamma_ij_tk[T_net,K_net]; // for link probability
  matrix[V_net,V_net] mu_ij_tk[T_net,K_net]; // for link weight
  
  // Link probabilities
  matrix<lower=0,upper=1>[V_net,V_net] pi_ij_tk[T_net,K_net];
  
  matrix[V_net,H_dim] ab_ih_t_shared[2,T_net];
  matrix[V_net,R_dim] ab_ih_tk[2,T_net,K_net];
  
  vector[T_net] eta_bar_t_k[K_net];
  vector[T_net] ab_bar_t_hi_shared[2,H_dim,V_net];
  vector[T_net] ab_bar_t_hki[2,R_dim,K_net,V_net];
  
  // Link weight
  matrix[V_net,H_dim] uv_ih_t_shared[2,T_net];
  matrix[V_net,R_dim] uv_ih_tk[2,T_net,K_net];
  
  vector[T_net] theta_bar_t_k[K_net];
  vector[T_net] uv_bar_t_hi_shared[2,H_dim,V_net];
  vector[T_net] uv_bar_t_hki[2,R_dim,K_net,V_net];
  
  for(dir in 1:2){
    
    // Shared latent coordinates //
    // rearrange 
    for(t in 1:T_net) {
      for(h in 1:H_dim) {
        for(i in 1:V_net) {
          ab_ih_t_shared[dir,t][i,h] = ab_t_hi_shared[dir,h,i][t];
          uv_ih_t_shared[dir,t][i,h] = uv_t_hi_shared[dir,h,i][t];
        }
      }
    }
  
    // Layer-specific latent coordinates //
    // rearrange
    for(t in 1:T_net) {
      for(k in 1:K_net) {
        for(h in 1:R_dim) {
          for(i in 1:V_net) {
            ab_ih_tk[dir,t,k][i,h] = ab_t_hki[dir,h,k,i][t];
            uv_ih_tk[dir,t,k][i,h] = uv_t_hki[dir,h,k,i][t];
          }
        }
      }
    }
    
    for(k in 1:K_net) {
      for(h in 1:R_dim) {
        for(i in 1:V_net) {
          eta_bar_t_k[k] = rep_vector(eta_bar_k[k],T_net);
          ab_bar_t_hi_shared[dir,h,i] = rep_vector(ab_bar_hi_shared[dir,h,i],T_net);
          ab_bar_t_hki[dir,h,k,i] = rep_vector(ab_bar_hki[dir,h,k,i],T_net);
          
          theta_bar_t_k[k] = rep_vector(theta_bar_k[k],T_net);
          uv_bar_t_hi_shared[dir,h,i] = rep_vector(uv_bar_hi_shared[dir,h,i],T_net);
          uv_bar_t_hki[dir,h,k,i] = rep_vector(uv_bar_hki[dir,h,k,i],T_net);
        }
      }
    }
    
  }
  
  // Linear predictors //
  for(t in 1:T_net) {
    for(k in 1:K_net) {
      // Baseline and Global coordinates
      gamma_ij_tk[t,k] = eta_t_k[k][t] + ab_ih_t_shared[1,t]*ab_ih_t_shared[2,t]';
      mu_ij_tk[t,k] = theta_t_k[k][t] + uv_ih_t_shared[1,t]*uv_ih_t_shared[2,t]';
      
      // Layer-specific coordinates
      if(K_net>1){
        gamma_ij_tk[t,k] += ab_ih_tk[1,t,k]*ab_ih_tk[2,t,k]';
        mu_ij_tk[t,k] += uv_ih_tk[1,t,k]*uv_ih_tk[2,t,k]';
      }
    }
  }
  
  // Link probabilities //
  for (t in 1:T_net) {
    for (k in 1:K_net) {
      pi_ij_tk[t,k] = inv_logit( gamma_ij_tk[t,k] );
    }
  }
}

model {
  // The likelihood
  for (t in 1:T_net) {
    for (k in 1:K_net) {
      for (i in 1:V_net) {
        for (j in 1:V_net) {
          if(i!=j){
            if( y_tkij[t,k,i,j]==0 ){
              target += log(1-pi_ij_tk[t,k][i,j]);
            } else if( y_tkij[t,k,i,j]!=1234567890 ){
              target += log(pi_ij_tk[t,k][i,j]) + normal_lpdf( y_tkij[t,k,i,j] | mu_ij_tk[t,k][i,j] , exp(log_sigma_w_k[k]) );
            }
          }
        }
      }
    }
  }
  
  // baseline processes //
  
  eta_t_k ~ multi_normal_cholesky( eta_bar_t_k , eta_t_cov_sqrt );
  theta_t_k ~ multi_normal_cholesky( theta_bar_t_k , theta_t_cov_sqrt );
  
  
  for(dir in 1:2){
    // Shared latent coordinates //
    for (h in 1:H_dim) {
      ab_t_hi_shared[dir,h] ~ multi_normal_cholesky( ab_bar_t_hi_shared[dir,h] , ab_t_cov_sqrt );
      uv_t_hi_shared[dir,h] ~ multi_normal_cholesky( uv_bar_t_hi_shared[dir,h] , uv_t_cov_sqrt );
      
    }
    
    // Layer-specific latent coordinates //
    for (h in 1:R_dim) {
      for (k in 1:K_net) {
        ab_t_hki[dir,h,k] ~ multi_normal_cholesky( ab_bar_t_hki[dir,h,k] , ab_t_cov_sqrt );
        uv_t_hki[dir,h,k] ~ multi_normal_cholesky( uv_bar_t_hki[dir,h,k] , uv_t_cov_sqrt );
      }
    }
    
  }
  
  // eta_bar_k ~ normal( 0.0 , sigma_lat_mean );
  // ab_bar_hi_shared ~ normal( 0.0 , sigma_lat_mean );
  // ab_bar_hki ~ normal( 0.0 , sigma_lat_mean );
  // theta_bar_k ~ normal( 0.0 , sigma_lat_mean );
  // uv_bar_hi_shared ~ normal( 0.0 , sigma_lat_mean );
  // uv_bar_hki ~ normal( 0.0 , sigma_lat_mean );
  
}

generated quantities {
  real<lower=0> sigma_w_k[K_net];
  sigma_w_k = exp(log_sigma_w_k);
}
