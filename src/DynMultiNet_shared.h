
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

arma::mat sample_x_ith_shared_DynMultiNet_bin_cpp( arma::cube x_ith_shared,
                                                   const arma::mat x_t_sigma_prior_inv,
                                                   const arma::mat tau_h,
                                                   const arma::field<arma::cube> y_ijtk,
                                                   const arma::field<arma::cube> w_ijtk,
                                                   const arma::field<arma::cube> s_ijtk );

Rcpp::List sample_x_ith_shared_DynMultiNet_bin_dir_cpp( arma::mat x_ith_shared_send,
                                                        arma::mat x_ith_shared_receive,
                                                        const arma::mat x_t_sigma_prior_inv,
                                                        const arma::mat tau_h_send,
                                                        const arma::mat tau_h_receive,
                                                        const arma::field<arma::cube> y_ijtk,
                                                        const arma::field<arma::cube> w_ijtk,
                                                        const arma::field<arma::cube> s_ijtk );

Rcpp::List sample_x_ith_shared_DynMultiNet_bin_dir_cpp( arma::cube x_ith_shared_send,
                                                        arma::cube x_ith_shared_receive,
                                                        const arma::mat x_t_sigma_prior_inv,
                                                        const arma::colvec tau_h_shared_send,
                                                        const arma::colvec tau_h_shared_receive,
                                                        const arma::field<arma::cube> y_ijtk,
                                                        const arma::field<arma::cube> w_ijtk,
                                                        const arma::field<arma::cube> s_ijtk );
