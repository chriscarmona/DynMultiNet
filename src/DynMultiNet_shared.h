
arma::mat sample_mu_t_DynMultiNet_bin_cpp( arma::colvec mu_t,
                                           const arma::mat mu_t_cov_prior_inv,
                                           const arma::cube y_ijt,
                                           const arma::cube w_ijt,
                                           const arma::cube s_ijt );

arma::mat sample_beta_z_layer_DynMultiNet_bin_cpp( arma::colvec beta_t,
                                                   arma::colvec z_t,
                                                   const arma::mat beta_t_cov_prior_inv,
                                                   const arma::cube y_ijt,
                                                   const arma::cube w_ijt,
                                                   const arma::cube s_ijt );

arma::mat sample_x_iht_mat_DynMultiNet_bin_cpp( arma::mat x_iht_mat,
                                                const arma::mat x_t_sigma_prior_inv,
                                                const arma::mat tau_h,
                                                const arma::cube y_ijt,
                                                const arma::cube w_ijt,
                                                const arma::cube s_ijt,
                                                const arma::mat mu_t );

