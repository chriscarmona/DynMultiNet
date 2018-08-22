
// saved as net_m_1_p_0_d_1_w_0.stan

// Supports Dynamic Networks with the following characteristics:
//   multi-layer: TRUE
//   predictors: FALSE
//   directed: TRUE
//   weigthed: FALSE

data {
  int<lower=1> V_net; // number of nodes
  int<lower=1> T_net; // number of time-steps in the observed network
  int<lower=1> K_net; // number of layers
  int<lower=1> H_dim; // number of dimensions in shared latent space
  int<lower=1> R_dim; // number of dimensions in layer-specific latent space
  
  int<lower=0,upper=1> y_tkij[T_net,K_net,V_net,V_net];
  
  real<lower=0> a_1;
  real<lower=0> a_2;
  
  real<lower=0> k_mu;
  real<lower=0> k_x;
  
  int<lower=1> T_all; // Total number of time-steps (incluiding forecast)
  vector[T_all] time_all; // Time stamps of all time-steps
  int<lower=1> time_all_idx_net[T_net]; // Indices of time_all that correspond to observed networks
}

transformed data {
  cov_matrix[T_all] mu_t_cov_prior;
  cov_matrix[T_all] x_t_cov_prior;
  cholesky_factor_cov[T_all] mu_t_cov_prior_chol;
  
  for(t1 in 1:T_all) {
    for(t2 in 1:t1) {
      mu_t_cov_prior[t1,t2] = exp(-k_mu*(time_all[t1]-time_all[t2])^2);
      mu_t_cov_prior[t2,t1] = exp(-k_mu*(time_all[t1]-time_all[t2])^2);
      x_t_cov_prior[t1,t2] = exp(-k_x*(time_all[t1]-time_all[t2])^2);
      x_t_cov_prior[t2,t1] = exp(-k_x*(time_all[t1]-time_all[t2])^2);
    }
  }
  
  mu_t_cov_prior_chol = cholesky_decompose(mu_t_cov_prior);
}

parameters {
  matrix[K_net,T_all] z_mu_kt;
  
  matrix[V_net,T_all] z_x_it_h_shared[2,H_dim];
  matrix[V_net,T_all] z_x_it_hk[2,R_dim,K_net];
  
  vector<lower=0>[H_dim] nu_h_shared[2];
  matrix<lower=0>[H_dim,K_net] nu_hk[2];
}

transformed parameters {
  matrix[V_net,V_net] s_ij_tk[T_all,K_net];
  
  matrix[T_all,K_net] mu_tk;
  
  matrix[T_all,V_net] x_ti_h_shared[2,H_dim];
  matrix[T_all,V_net] x_ti_hk[2,R_dim,K_net];
  
  matrix[V_net,H_dim] x_ih_t_shared[2,T_all];
  matrix[V_net,R_dim] x_ih_tk[2,T_all,K_net];
  
  vector[H_dim] tau_h_shared[2];
  matrix[H_dim,K_net] tau_hk[2];
  cholesky_factor_cov[T_all] x_t_shared_cov_chol[2,H_dim];
  cholesky_factor_cov[T_all] x_t_k_cov_chol[2,H_dim,K_net];
  
  // baseline process //
  mu_tk = mu_t_cov_prior_chol * z_mu_kt';
  
  // Shrinkage factors and covariance matrices //
  for(dir in 1:2){
    tau_h_shared[dir][1]=nu_h_shared[dir][1];
    x_t_shared_cov_chol[dir,1] = cholesky_decompose( (1/tau_h_shared[dir][1])*x_t_cov_prior );
    tau_hk[dir][1]=nu_hk[dir][1];
    for (k in 1:K_net) {
      x_t_k_cov_chol[dir,1,k] = cholesky_decompose( (1/tau_hk[dir][1,k])*x_t_cov_prior );
    }
    for(h in 2:H_dim) {
      tau_h_shared[dir][h] = tau_h_shared[dir][h-1]*nu_h_shared[dir][h];
      x_t_shared_cov_chol[dir,h] = cholesky_decompose( (1/tau_h_shared[dir][h])*x_t_cov_prior );
      for (k in 1:K_net) {
        tau_hk[dir][h,k] = tau_hk[dir][h-1,k] * nu_hk[dir][h,k];
        x_t_k_cov_chol[dir,h,k] = cholesky_decompose( (1/tau_hk[dir][h,k])*x_t_cov_prior );
      }
    }
  
    // Shared latent coordinates //
    for(h in 1:H_dim) {
      x_ti_h_shared[dir,h] = x_t_shared_cov_chol[dir,h] * z_x_it_h_shared[dir,h]';
    }
  
    // rearrange 
    for(t in 1:T_all) {
      for(h in 1:R_dim) {
        for(i in 1:V_net) {
          x_ih_t_shared[dir,t][i,h] = x_ti_h_shared[dir,h][t,i];
        }
      }
    }
  
    // Layer-specific latent coordinates //
    for(h in 1:R_dim) {
      for(k in 1:K_net) {
        x_ti_hk[dir,h,k] = x_t_k_cov_chol[dir,h,k] * z_x_it_hk[dir,h,k]';
      }
    }
  
    // rearrange
    for(t in 1:T_all) {
      for(k in 1:K_net) {
        for(h in 1:R_dim) {
          for(i in 1:V_net) {
            x_ih_tk[dir,t,k][i,h] = x_ti_hk[dir,h,k][t,i];
          }
        }
      }
    }
  }
  
  // Linear predictor //
  for(t in 1:T_all) {
    for(k in 1:K_net) {
      s_ij_tk[t,k] = mu_tk[t,k] + x_ih_t_shared[1,t]*x_ih_t_shared[2,t]' + x_ih_tk[1,t,k]*x_ih_tk[2,t,k]';
    }
  }
  
}

model {
  // The likelihood
  for (t in 1:T_net) {
    for (k in 1:K_net) {
      for (i in 1:V_net) {
        if(i>1){
          y_tkij[t,k,i,1:(i-1)] ~ bernoulli_logit( s_ij_tk[time_all_idx_net[t],k][i,1:(i-1)] );
        }
        if(i<V_net){
          y_tkij[t,k,i,(i+1):V_net] ~ bernoulli_logit( s_ij_tk[time_all_idx_net[t],k][i,(i+1):V_net] );
        }
      }
    }
  }
  
  // baseline process //
  for (k in 1:K_net) {
    z_mu_kt[k] ~ normal( 0 , 1 );
  }
  
  for(dir in 1:2){
    // Shared latent coordinates //
    for (h in 1:H_dim) {
      for (i in 1:V_net) {
        z_x_it_h_shared[dir,h][i] ~ normal( 0 , 1 );
      }
    }
    
    // Layer-specific latent coordinates //
    for (h in 1:H_dim) {
      for (k in 1:K_net) {
        for (i in 1:V_net) {
          z_x_it_hk[dir,h,k][i] ~ normal( 0 , 1 );
        }
      }
    }
    
    // Shrinkage factors //
    nu_h_shared[dir][1] ~ gamma(a_1,1);
    nu_hk[dir][1]~gamma(a_1,1);
    for(h in 2:H_dim) {
      nu_h_shared[dir][h] ~ gamma(a_2,1);
      nu_hk[dir][h] ~ gamma(a_2,1);
    }
  }
}

generated quantities {
  matrix[V_net,V_net] pi_ij_tk[T_all,K_net];
  for (t in 1:T_all) {
    for (k in 1:K_net) {
      pi_ij_tk[t,k] = inv_logit( s_ij_tk[t,k] );
    }
  }
}
