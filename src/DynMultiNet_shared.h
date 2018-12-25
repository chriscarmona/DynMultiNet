
Rcpp::List sample_baseline_t_link_cpp( arma::colvec mu_t,
                                       const arma::mat mu_t_cov_prior_inv,
                                       const arma::cube y_ijt,
                                       const arma::cube w_ijt,
                                       arma::cube gamma_ijt,
                                       const bool directed );
