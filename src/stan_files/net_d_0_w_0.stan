
// saved as net_d_0_w_0.stan

data {
  int<lower=0> V_net; // number of nodes
  int<lower=0> T_net; // number of time-steps
  int<lower=0> K_net; // number of layers
  int<lower=0> H_dim; // number of dimensions in shared latent space
  int<lower=0> R_dim; // number of dimensions in layer-specific latent space
  
  int<lower=0> y_ijtk[V_net,V_net,T_net,K_net];
  
  vector[T_net] mu_t_mean_prior;
  corr_matrix[T_net] mu_t_cov_prior;
  
  vector[T_net] x_t_mean_prior;
  corr_matrix[T_net] x_t_cov_prior;
  
  real<lower=0> a_1;
  real<lower=0> a_2;
}

parameters {
  vector[T_net] x_t_ih_shared[V_net,H_dim];
  vector[T_net] x_t_ihk[V_net,R_dim,K_net];
  vector[T_net] mu_kt[K_net];
  vector<lower=0>[H_dim] nu_h_shared;
  matrix<lower=0>[H_dim,K_net] nu_hk;
}

transformed parameters {
  // real pi_ijtk[V_net,V_net,T_net,K_net];
  vector[H_dim] x_h_it_shared[V_net,T_net];
  vector[H_dim] x_h_itk[V_net,T_net,K_net];
  real s_ijtk[V_net,V_net,T_net,K_net];
  vector[H_dim] tau_h_shared;
  matrix[H_dim,K_net] tau_hk;
  
  for (t in 1:T_net) {
    for (i in 1:V_net) {
      for (h in 1:H_dim) {
        x_h_it_shared[i,t][h] = x_t_ih_shared[i,h][t];
        for (k in 1:K_net) {
          x_h_itk[i,t,k][h] = x_t_ihk[i,h,k][t];
        }
      }
    }
  }
  
  for (k in 1:K_net) {
    for (t in 1:T_net) {
      for (i in 2:V_net) {
        for (j in 1:(i-1)) {
          s_ijtk[i,j,t,k] = mu_kt[k][t] + x_h_it_shared[i,t]' * x_h_it_shared[j,t] + x_h_itk[i,t,k]' * x_h_itk[j,t,k];
        }
      }
    }
  }
  
  tau_h_shared[1]=nu_h_shared[1];
  tau_hk[1]=nu_hk[1];
  for(h in 2:H_dim) {
    tau_h_shared[h] = tau_h_shared[h-1]*nu_h_shared[h];
    tau_hk[h] = tau_hk[h-1] .* nu_hk[h];
  }
  
}

model {
  // The likelihood
  for (k in 1:K_net) {
    for (t in 1:T_net) {
      for (i in 2:V_net) {
        for (j in 1:(i-1)) {
          y_ijtk[i,j,t,k] ~ bernoulli_logit( s_ijtk[i,j,t,k] );
        }
      }
    }
  }
  // y_ijtk ~ bernoulli_logit( s_ijtk );
  
  for (k in 1:K_net) {
    mu_kt[k] ~ multi_normal( mu_t_mean_prior , mu_t_cov_prior );
  }
  
  for (i in 1:V_net) {
    for (h in 1:H_dim) {
      x_t_ih_shared[i,h] ~ multi_normal( x_t_mean_prior , (1/tau_h_shared[h])*x_t_cov_prior );
    }
    for (k in 1:K_net) {
      for (h in 1:R_dim) {
        x_t_ihk[i,h,k] ~ multi_normal( x_t_mean_prior , (1/tau_hk[h,k])*x_t_cov_prior );
      }
    }
  }
  
  nu_h_shared[1] ~ gamma(a_1,1);
  nu_hk[1]~gamma(a_1,1);
  for(h in 2:H_dim) {
    nu_h_shared[h] ~ gamma(a_2,1);
    nu_hk[h] ~ gamma(a_2,1);
  }
}
