
// saved as net_d_0_w_0.stan

data {
  int<lower=1> V_net; // number of nodes
  int<lower=1> T_net; // number of time-steps
  int<lower=1> K_net; // number of layers
  int<lower=1> H_dim; // number of dimensions in shared latent space
  int<lower=1> R_dim; // number of dimensions in layer-specific latent space
  
  int<lower=0,upper=1> y_ijtk[V_net,V_net,T_net,K_net];
  
  cov_matrix[T_net] mu_t_cov_prior;
  cov_matrix[T_net] x_t_cov_prior;
  
  real<lower=0> a_1;
  real<lower=0> a_2;
}

transformed data {
  cholesky_factor_cov[T_net] mu_t_cov_prior_chol = cholesky_decompose(mu_t_cov_prior);
}

parameters {
  matrix[V_net,T_net] z_x_it_h_shared[H_dim];
  matrix[V_net,T_net] z_x_it_hk[R_dim,K_net];
  matrix[K_net,T_net] z_mu_kt;
  vector<lower=0>[H_dim] nu_h_shared;
  matrix<lower=0>[H_dim,K_net] nu_hk;
}

transformed parameters {
  // matrix[V_net,V_net] pi_ij_tk[T_net,K_net];
  matrix[V_net,V_net] s_ij_tk[T_net,K_net];
  
  matrix[T_net,K_net] mu_tk;
  
  matrix[T_net,V_net] x_ti_h_shared[H_dim];
  matrix[T_net,V_net] x_ti_hk[R_dim,K_net];
  
  matrix[V_net,H_dim] x_ih_t_shared[T_net];
  matrix[V_net,R_dim] x_ih_tk[T_net,K_net];
  
  vector[H_dim] tau_h_shared;
  matrix[H_dim,K_net] tau_hk;
  cholesky_factor_cov[T_net] x_t_shared_cov_chol[H_dim];
  cholesky_factor_cov[T_net] x_t_k_cov_chol[H_dim,K_net];
  
  // baseline process //
  mu_tk = mu_t_cov_prior_chol * z_mu_kt';
  
  // Shrinkage factors and covariance matrices //
  tau_h_shared[1]=nu_h_shared[1];
  x_t_shared_cov_chol[1] = cholesky_decompose( (1/tau_h_shared[1])*x_t_cov_prior );
  tau_hk[1]=nu_hk[1];
  for (k in 1:K_net) {
    x_t_k_cov_chol[1,k] = cholesky_decompose( (1/tau_hk[1,k])*x_t_cov_prior );
  }
  for(h in 2:H_dim) {
    tau_h_shared[h] = tau_h_shared[h-1]*nu_h_shared[h];
    x_t_shared_cov_chol[h] = cholesky_decompose( (1/tau_h_shared[h])*x_t_cov_prior );
    for (k in 1:K_net) {
      tau_hk[h,k] = tau_hk[h-1,k] .* nu_hk[h,k];
      x_t_k_cov_chol[h,k] = cholesky_decompose( (1/tau_hk[h,k])*x_t_cov_prior );
    }
  }
  
  // Shared latent coordinates //
  for(h in 1:H_dim) {
    x_ti_h_shared[h] = x_t_shared_cov_chol[h] * z_x_it_h_shared[h]';
  }
  // rearrange 
  for(t in 1:T_net) {
    for(h in 1:R_dim) {
      for(i in 1:V_net) {
        x_ih_t_shared[t][i,h] = x_ti_h_shared[h][t,i];
      }
    }
  }
  
  // Layer-specific latent coordinates //
  for(h in 1:R_dim) {
    for(k in 1:K_net) {
      x_ti_hk[h,k] = x_t_k_cov_chol[h,k] * z_x_it_hk[h,k]';
    }
  }
  // rearrange
  for(t in 1:T_net) {
    for(k in 1:K_net) {
      for(h in 1:R_dim) {
        for(i in 1:V_net) {
          x_ih_tk[t,k][i,h] = x_ti_hk[h,k][t,i];
        }
      }
    }
  }
  
  // Linear predictor //
  for(t in 1:T_net) {
    for(k in 1:K_net) {
      s_ij_tk[t,k] = mu_tk[t,k] + tcrossprod(x_ih_t_shared[t]) + tcrossprod(x_ih_tk[t,k]);
    }
  }
}

model {
  // The likelihood
  for (t in 1:T_net) {
    for (k in 1:K_net) {
      for (i in 2:V_net) {
        for (j in 1:(i-1)) {
          y_ijtk[i,j,t,k] ~ bernoulli_logit( s_ij_tk[t,k][i,j] );
        }
      }
    }
  }
  
  // baseline process //
  for (k in 1:K_net) {
    z_mu_kt[k] ~ normal( 0 , 1 );
  }
  
  // Shared latent coordinates //
  for (h in 1:H_dim) {
    for (i in 1:V_net) {
      z_x_it_h_shared[h][i] ~ normal( 0 , 1 );
    }
  }
  
  // Layer-specific latent coordinates //
  for (h in 1:H_dim) {
    for (k in 1:K_net) {
      for (i in 1:V_net) {
        z_x_it_hk[h,k][i] ~ normal( 0 , 1 );
      }
    }
  }
  
  // Shrinkage factors //
  nu_h_shared[1] ~ gamma(a_1,1);
  nu_hk[1]~gamma(a_1,1);
  for(h in 2:H_dim) {
    nu_h_shared[h] ~ gamma(a_2,1);
    nu_hk[h] ~ gamma(a_2,1);
  }
}

generated quantities {
  matrix[V_net,V_net] pi_ij_tk[T_net,K_net];
  for (t in 1:T_net) {
    for (k in 1:K_net) {
      pi_ij_tk[t,k] = inv_logit( s_ij_tk[t,k] );
    }
  }
}
