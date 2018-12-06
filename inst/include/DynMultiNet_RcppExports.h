// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_DynMultiNet_RCPPEXPORTS_H_GEN_
#define RCPP_DynMultiNet_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace DynMultiNet {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("DynMultiNet", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("DynMultiNet", "_DynMultiNet_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in DynMultiNet");
            }
        }
    }

    inline arma::mat sample_mu_t_DynMultiNet_bin_cpp(arma::colvec mu_t, const arma::mat mu_t_cov_prior_inv, const arma::cube y_ijt, const arma::cube w_ijt, const arma::cube s_ijt, const bool directed = false) {
        typedef SEXP(*Ptr_sample_mu_t_DynMultiNet_bin_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_mu_t_DynMultiNet_bin_cpp p_sample_mu_t_DynMultiNet_bin_cpp = NULL;
        if (p_sample_mu_t_DynMultiNet_bin_cpp == NULL) {
            validateSignature("arma::mat(*sample_mu_t_DynMultiNet_bin_cpp)(arma::colvec,const arma::mat,const arma::cube,const arma::cube,const arma::cube,const bool)");
            p_sample_mu_t_DynMultiNet_bin_cpp = (Ptr_sample_mu_t_DynMultiNet_bin_cpp)R_GetCCallable("DynMultiNet", "_DynMultiNet_sample_mu_t_DynMultiNet_bin_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_mu_t_DynMultiNet_bin_cpp(Shield<SEXP>(Rcpp::wrap(mu_t)), Shield<SEXP>(Rcpp::wrap(mu_t_cov_prior_inv)), Shield<SEXP>(Rcpp::wrap(y_ijt)), Shield<SEXP>(Rcpp::wrap(w_ijt)), Shield<SEXP>(Rcpp::wrap(s_ijt)), Shield<SEXP>(Rcpp::wrap(directed)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat sample_mu_t_DynMultiNet_bin_v2_cpp(arma::colvec mu_t, const arma::mat mu_t_cov_prior_inv, const arma::cube y_ijt, const arma::cube w_ijt, const arma::cube s_ijt, const bool directed = false) {
        typedef SEXP(*Ptr_sample_mu_t_DynMultiNet_bin_v2_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_mu_t_DynMultiNet_bin_v2_cpp p_sample_mu_t_DynMultiNet_bin_v2_cpp = NULL;
        if (p_sample_mu_t_DynMultiNet_bin_v2_cpp == NULL) {
            validateSignature("arma::mat(*sample_mu_t_DynMultiNet_bin_v2_cpp)(arma::colvec,const arma::mat,const arma::cube,const arma::cube,const arma::cube,const bool)");
            p_sample_mu_t_DynMultiNet_bin_v2_cpp = (Ptr_sample_mu_t_DynMultiNet_bin_v2_cpp)R_GetCCallable("DynMultiNet", "_DynMultiNet_sample_mu_t_DynMultiNet_bin_v2_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_mu_t_DynMultiNet_bin_v2_cpp(Shield<SEXP>(Rcpp::wrap(mu_t)), Shield<SEXP>(Rcpp::wrap(mu_t_cov_prior_inv)), Shield<SEXP>(Rcpp::wrap(y_ijt)), Shield<SEXP>(Rcpp::wrap(w_ijt)), Shield<SEXP>(Rcpp::wrap(s_ijt)), Shield<SEXP>(Rcpp::wrap(directed)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::colvec sample_beta_z_layer_DynMultiNet_bin_cpp(arma::colvec beta_t, arma::colvec z_t, const arma::mat beta_t_cov_prior_inv, const arma::cube y_ijt, const arma::cube w_ijt, const arma::cube s_ijt, const bool directed = false) {
        typedef SEXP(*Ptr_sample_beta_z_layer_DynMultiNet_bin_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_beta_z_layer_DynMultiNet_bin_cpp p_sample_beta_z_layer_DynMultiNet_bin_cpp = NULL;
        if (p_sample_beta_z_layer_DynMultiNet_bin_cpp == NULL) {
            validateSignature("arma::colvec(*sample_beta_z_layer_DynMultiNet_bin_cpp)(arma::colvec,arma::colvec,const arma::mat,const arma::cube,const arma::cube,const arma::cube,const bool)");
            p_sample_beta_z_layer_DynMultiNet_bin_cpp = (Ptr_sample_beta_z_layer_DynMultiNet_bin_cpp)R_GetCCallable("DynMultiNet", "_DynMultiNet_sample_beta_z_layer_DynMultiNet_bin_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_beta_z_layer_DynMultiNet_bin_cpp(Shield<SEXP>(Rcpp::wrap(beta_t)), Shield<SEXP>(Rcpp::wrap(z_t)), Shield<SEXP>(Rcpp::wrap(beta_t_cov_prior_inv)), Shield<SEXP>(Rcpp::wrap(y_ijt)), Shield<SEXP>(Rcpp::wrap(w_ijt)), Shield<SEXP>(Rcpp::wrap(s_ijt)), Shield<SEXP>(Rcpp::wrap(directed)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::colvec >(rcpp_result_gen);
    }

    inline arma::cube sample_x_ith_DynMultiNet_bin_cpp(arma::cube x_ith, const arma::mat x_t_sigma_prior_inv, const arma::colvec tau_h, const arma::cube y_ijt, const arma::cube w_ijt, const arma::cube s_ijt) {
        typedef SEXP(*Ptr_sample_x_ith_DynMultiNet_bin_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_x_ith_DynMultiNet_bin_cpp p_sample_x_ith_DynMultiNet_bin_cpp = NULL;
        if (p_sample_x_ith_DynMultiNet_bin_cpp == NULL) {
            validateSignature("arma::cube(*sample_x_ith_DynMultiNet_bin_cpp)(arma::cube,const arma::mat,const arma::colvec,const arma::cube,const arma::cube,const arma::cube)");
            p_sample_x_ith_DynMultiNet_bin_cpp = (Ptr_sample_x_ith_DynMultiNet_bin_cpp)R_GetCCallable("DynMultiNet", "_DynMultiNet_sample_x_ith_DynMultiNet_bin_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_x_ith_DynMultiNet_bin_cpp(Shield<SEXP>(Rcpp::wrap(x_ith)), Shield<SEXP>(Rcpp::wrap(x_t_sigma_prior_inv)), Shield<SEXP>(Rcpp::wrap(tau_h)), Shield<SEXP>(Rcpp::wrap(y_ijt)), Shield<SEXP>(Rcpp::wrap(w_ijt)), Shield<SEXP>(Rcpp::wrap(s_ijt)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline arma::cube sample_x_ith_shared_DynMultiNet_bin_cpp(arma::cube x_ith_shared, const arma::mat x_t_sigma_prior_inv, const arma::colvec tau_h, const arma::field<arma::cube> y_ijtk, const arma::field<arma::cube> w_ijtk, const arma::field<arma::cube> s_ijtk) {
        typedef SEXP(*Ptr_sample_x_ith_shared_DynMultiNet_bin_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_x_ith_shared_DynMultiNet_bin_cpp p_sample_x_ith_shared_DynMultiNet_bin_cpp = NULL;
        if (p_sample_x_ith_shared_DynMultiNet_bin_cpp == NULL) {
            validateSignature("arma::cube(*sample_x_ith_shared_DynMultiNet_bin_cpp)(arma::cube,const arma::mat,const arma::colvec,const arma::field<arma::cube>,const arma::field<arma::cube>,const arma::field<arma::cube>)");
            p_sample_x_ith_shared_DynMultiNet_bin_cpp = (Ptr_sample_x_ith_shared_DynMultiNet_bin_cpp)R_GetCCallable("DynMultiNet", "_DynMultiNet_sample_x_ith_shared_DynMultiNet_bin_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_x_ith_shared_DynMultiNet_bin_cpp(Shield<SEXP>(Rcpp::wrap(x_ith_shared)), Shield<SEXP>(Rcpp::wrap(x_t_sigma_prior_inv)), Shield<SEXP>(Rcpp::wrap(tau_h)), Shield<SEXP>(Rcpp::wrap(y_ijtk)), Shield<SEXP>(Rcpp::wrap(w_ijtk)), Shield<SEXP>(Rcpp::wrap(s_ijtk)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline Rcpp::List sample_x_ith_DynMultiNet_bin_dir_cpp(arma::cube x_ith_send, arma::cube x_ith_receive, const arma::mat x_t_sigma_prior_inv, const arma::colvec tau_h_send, const arma::colvec tau_h_receive, const arma::cube y_ijt, const arma::cube w_ijt, const arma::cube s_ijt) {
        typedef SEXP(*Ptr_sample_x_ith_DynMultiNet_bin_dir_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_x_ith_DynMultiNet_bin_dir_cpp p_sample_x_ith_DynMultiNet_bin_dir_cpp = NULL;
        if (p_sample_x_ith_DynMultiNet_bin_dir_cpp == NULL) {
            validateSignature("Rcpp::List(*sample_x_ith_DynMultiNet_bin_dir_cpp)(arma::cube,arma::cube,const arma::mat,const arma::colvec,const arma::colvec,const arma::cube,const arma::cube,const arma::cube)");
            p_sample_x_ith_DynMultiNet_bin_dir_cpp = (Ptr_sample_x_ith_DynMultiNet_bin_dir_cpp)R_GetCCallable("DynMultiNet", "_DynMultiNet_sample_x_ith_DynMultiNet_bin_dir_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_x_ith_DynMultiNet_bin_dir_cpp(Shield<SEXP>(Rcpp::wrap(x_ith_send)), Shield<SEXP>(Rcpp::wrap(x_ith_receive)), Shield<SEXP>(Rcpp::wrap(x_t_sigma_prior_inv)), Shield<SEXP>(Rcpp::wrap(tau_h_send)), Shield<SEXP>(Rcpp::wrap(tau_h_receive)), Shield<SEXP>(Rcpp::wrap(y_ijt)), Shield<SEXP>(Rcpp::wrap(w_ijt)), Shield<SEXP>(Rcpp::wrap(s_ijt)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::List sample_x_ith_shared_DynMultiNet_bin_dir_cpp(arma::cube x_ith_shared_send, arma::cube x_ith_shared_receive, const arma::mat x_t_sigma_prior_inv, const arma::colvec tau_h_shared_send, const arma::colvec tau_h_shared_receive, const arma::field<arma::cube> y_ijtk, const arma::field<arma::cube> w_ijtk, const arma::field<arma::cube> s_ijtk) {
        typedef SEXP(*Ptr_sample_x_ith_shared_DynMultiNet_bin_dir_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_x_ith_shared_DynMultiNet_bin_dir_cpp p_sample_x_ith_shared_DynMultiNet_bin_dir_cpp = NULL;
        if (p_sample_x_ith_shared_DynMultiNet_bin_dir_cpp == NULL) {
            validateSignature("Rcpp::List(*sample_x_ith_shared_DynMultiNet_bin_dir_cpp)(arma::cube,arma::cube,const arma::mat,const arma::colvec,const arma::colvec,const arma::field<arma::cube>,const arma::field<arma::cube>,const arma::field<arma::cube>)");
            p_sample_x_ith_shared_DynMultiNet_bin_dir_cpp = (Ptr_sample_x_ith_shared_DynMultiNet_bin_dir_cpp)R_GetCCallable("DynMultiNet", "_DynMultiNet_sample_x_ith_shared_DynMultiNet_bin_dir_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_x_ith_shared_DynMultiNet_bin_dir_cpp(Shield<SEXP>(Rcpp::wrap(x_ith_shared_send)), Shield<SEXP>(Rcpp::wrap(x_ith_shared_receive)), Shield<SEXP>(Rcpp::wrap(x_t_sigma_prior_inv)), Shield<SEXP>(Rcpp::wrap(tau_h_shared_send)), Shield<SEXP>(Rcpp::wrap(tau_h_shared_receive)), Shield<SEXP>(Rcpp::wrap(y_ijtk)), Shield<SEXP>(Rcpp::wrap(w_ijtk)), Shield<SEXP>(Rcpp::wrap(s_ijtk)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline SEXP c_initialize1(SEXP Data, SEXP DIMS, SEXP Yy, SEXP XSCALE, SEXP BETAIN, SEXP BETAOUT, SEXP WW) {
        typedef SEXP(*Ptr_c_initialize1)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_c_initialize1 p_c_initialize1 = NULL;
        if (p_c_initialize1 == NULL) {
            validateSignature("SEXP(*c_initialize1)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP)");
            p_c_initialize1 = (Ptr_c_initialize1)R_GetCCallable("DynMultiNet", "_DynMultiNet_c_initialize1");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_c_initialize1(Shield<SEXP>(Rcpp::wrap(Data)), Shield<SEXP>(Rcpp::wrap(DIMS)), Shield<SEXP>(Rcpp::wrap(Yy)), Shield<SEXP>(Rcpp::wrap(XSCALE)), Shield<SEXP>(Rcpp::wrap(BETAIN)), Shield<SEXP>(Rcpp::wrap(BETAOUT)), Shield<SEXP>(Rcpp::wrap(WW)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline SEXP c_initialize1_grad(SEXP Data, SEXP DIMS, SEXP Yy, SEXP XSCALE, SEXP BETAIN, SEXP BETAOUT, SEXP WW) {
        typedef SEXP(*Ptr_c_initialize1_grad)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_c_initialize1_grad p_c_initialize1_grad = NULL;
        if (p_c_initialize1_grad == NULL) {
            validateSignature("SEXP(*c_initialize1_grad)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP)");
            p_c_initialize1_grad = (Ptr_c_initialize1_grad)R_GetCCallable("DynMultiNet", "_DynMultiNet_c_initialize1_grad");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_c_initialize1_grad(Shield<SEXP>(Rcpp::wrap(Data)), Shield<SEXP>(Rcpp::wrap(DIMS)), Shield<SEXP>(Rcpp::wrap(Yy)), Shield<SEXP>(Rcpp::wrap(XSCALE)), Shield<SEXP>(Rcpp::wrap(BETAIN)), Shield<SEXP>(Rcpp::wrap(BETAOUT)), Shield<SEXP>(Rcpp::wrap(WW)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline SEXP c_update2(SEXP Xitm1, SEXP DIMS, SEXP TUNEX, SEXP Yy, SEXP BETAIN, SEXP BETAOUT, SEXP TUNEBIO, SEXP WW, SEXP t2X, SEXP s2X, SEXP xiBIN, SEXP xiBOUT, SEXP nuBIN, SEXP nuBOUT, SEXP CAUCHY, SEXP RNORMS, SEXP RNORMSBIO, SEXP ELOUT, SEXP ELIN, SEXP SUBSEQ, SEXP DEG) {
        typedef SEXP(*Ptr_c_update2)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_c_update2 p_c_update2 = NULL;
        if (p_c_update2 == NULL) {
            validateSignature("SEXP(*c_update2)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP)");
            p_c_update2 = (Ptr_c_update2)R_GetCCallable("DynMultiNet", "_DynMultiNet_c_update2");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_c_update2(Shield<SEXP>(Rcpp::wrap(Xitm1)), Shield<SEXP>(Rcpp::wrap(DIMS)), Shield<SEXP>(Rcpp::wrap(TUNEX)), Shield<SEXP>(Rcpp::wrap(Yy)), Shield<SEXP>(Rcpp::wrap(BETAIN)), Shield<SEXP>(Rcpp::wrap(BETAOUT)), Shield<SEXP>(Rcpp::wrap(TUNEBIO)), Shield<SEXP>(Rcpp::wrap(WW)), Shield<SEXP>(Rcpp::wrap(t2X)), Shield<SEXP>(Rcpp::wrap(s2X)), Shield<SEXP>(Rcpp::wrap(xiBIN)), Shield<SEXP>(Rcpp::wrap(xiBOUT)), Shield<SEXP>(Rcpp::wrap(nuBIN)), Shield<SEXP>(Rcpp::wrap(nuBOUT)), Shield<SEXP>(Rcpp::wrap(CAUCHY)), Shield<SEXP>(Rcpp::wrap(RNORMS)), Shield<SEXP>(Rcpp::wrap(RNORMSBIO)), Shield<SEXP>(Rcpp::wrap(ELOUT)), Shield<SEXP>(Rcpp::wrap(ELIN)), Shield<SEXP>(Rcpp::wrap(SUBSEQ)), Shield<SEXP>(Rcpp::wrap(DEG)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline SEXP c_update1(SEXP Xitm1, SEXP DIMS, SEXP TUNEX, SEXP Yy, SEXP BETAIN, SEXP BETAOUT, SEXP TUNEBIO, SEXP WW, SEXP t2X, SEXP s2X, SEXP xiBIN, SEXP xiBOUT, SEXP nuBIN, SEXP nuBOUT, SEXP CAUCHY, SEXP RNORMS, SEXP RNORMSBIO) {
        typedef SEXP(*Ptr_c_update1)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_c_update1 p_c_update1 = NULL;
        if (p_c_update1 == NULL) {
            validateSignature("SEXP(*c_update1)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP)");
            p_c_update1 = (Ptr_c_update1)R_GetCCallable("DynMultiNet", "_DynMultiNet_c_update1");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_c_update1(Shield<SEXP>(Rcpp::wrap(Xitm1)), Shield<SEXP>(Rcpp::wrap(DIMS)), Shield<SEXP>(Rcpp::wrap(TUNEX)), Shield<SEXP>(Rcpp::wrap(Yy)), Shield<SEXP>(Rcpp::wrap(BETAIN)), Shield<SEXP>(Rcpp::wrap(BETAOUT)), Shield<SEXP>(Rcpp::wrap(TUNEBIO)), Shield<SEXP>(Rcpp::wrap(WW)), Shield<SEXP>(Rcpp::wrap(t2X)), Shield<SEXP>(Rcpp::wrap(s2X)), Shield<SEXP>(Rcpp::wrap(xiBIN)), Shield<SEXP>(Rcpp::wrap(xiBOUT)), Shield<SEXP>(Rcpp::wrap(nuBIN)), Shield<SEXP>(Rcpp::wrap(nuBOUT)), Shield<SEXP>(Rcpp::wrap(CAUCHY)), Shield<SEXP>(Rcpp::wrap(RNORMS)), Shield<SEXP>(Rcpp::wrap(RNORMSBIO)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline SEXP c_t2s2Parms(SEXP DATA, SEXP DIMS, SEXP THETAT, SEXP THETAS, SEXP PHIT, SEXP PHIS) {
        typedef SEXP(*Ptr_c_t2s2Parms)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_c_t2s2Parms p_c_t2s2Parms = NULL;
        if (p_c_t2s2Parms == NULL) {
            validateSignature("SEXP(*c_t2s2Parms)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP)");
            p_c_t2s2Parms = (Ptr_c_t2s2Parms)R_GetCCallable("DynMultiNet", "_DynMultiNet_c_t2s2Parms");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_c_t2s2Parms(Shield<SEXP>(Rcpp::wrap(DATA)), Shield<SEXP>(Rcpp::wrap(DIMS)), Shield<SEXP>(Rcpp::wrap(THETAT)), Shield<SEXP>(Rcpp::wrap(THETAS)), Shield<SEXP>(Rcpp::wrap(PHIT)), Shield<SEXP>(Rcpp::wrap(PHIS)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline SEXP c_WAccProb2(SEXP Data, SEXP DIMS, SEXP Yy, SEXP BETAIN, SEXP BETAOUT, SEXP TUNEW, SEXP WWOld, SEXP WWNew, SEXP ELOUT, SEXP ELIN, SEXP SUBSEQ, SEXP DEG) {
        typedef SEXP(*Ptr_c_WAccProb2)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_c_WAccProb2 p_c_WAccProb2 = NULL;
        if (p_c_WAccProb2 == NULL) {
            validateSignature("SEXP(*c_WAccProb2)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP)");
            p_c_WAccProb2 = (Ptr_c_WAccProb2)R_GetCCallable("DynMultiNet", "_DynMultiNet_c_WAccProb2");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_c_WAccProb2(Shield<SEXP>(Rcpp::wrap(Data)), Shield<SEXP>(Rcpp::wrap(DIMS)), Shield<SEXP>(Rcpp::wrap(Yy)), Shield<SEXP>(Rcpp::wrap(BETAIN)), Shield<SEXP>(Rcpp::wrap(BETAOUT)), Shield<SEXP>(Rcpp::wrap(TUNEW)), Shield<SEXP>(Rcpp::wrap(WWOld)), Shield<SEXP>(Rcpp::wrap(WWNew)), Shield<SEXP>(Rcpp::wrap(ELOUT)), Shield<SEXP>(Rcpp::wrap(ELIN)), Shield<SEXP>(Rcpp::wrap(SUBSEQ)), Shield<SEXP>(Rcpp::wrap(DEG)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline SEXP c_WAccProb1(SEXP Data, SEXP DIMS, SEXP Yy, SEXP BETAIN, SEXP BETAOUT, SEXP TUNEW, SEXP WWOld, SEXP WWNew) {
        typedef SEXP(*Ptr_c_WAccProb1)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_c_WAccProb1 p_c_WAccProb1 = NULL;
        if (p_c_WAccProb1 == NULL) {
            validateSignature("SEXP(*c_WAccProb1)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP)");
            p_c_WAccProb1 = (Ptr_c_WAccProb1)R_GetCCallable("DynMultiNet", "_DynMultiNet_c_WAccProb1");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_c_WAccProb1(Shield<SEXP>(Rcpp::wrap(Data)), Shield<SEXP>(Rcpp::wrap(DIMS)), Shield<SEXP>(Rcpp::wrap(Yy)), Shield<SEXP>(Rcpp::wrap(BETAIN)), Shield<SEXP>(Rcpp::wrap(BETAOUT)), Shield<SEXP>(Rcpp::wrap(TUNEW)), Shield<SEXP>(Rcpp::wrap(WWOld)), Shield<SEXP>(Rcpp::wrap(WWNew)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline SEXP c_missing(SEXP Data, SEXP DIMS, SEXP MMM, SEXP Yy, SEXP Ttt, SEXP BETAIN, SEXP BETAOUT, SEXP WW) {
        typedef SEXP(*Ptr_c_missing)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_c_missing p_c_missing = NULL;
        if (p_c_missing == NULL) {
            validateSignature("SEXP(*c_missing)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP)");
            p_c_missing = (Ptr_c_missing)R_GetCCallable("DynMultiNet", "_DynMultiNet_c_missing");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_c_missing(Shield<SEXP>(Rcpp::wrap(Data)), Shield<SEXP>(Rcpp::wrap(DIMS)), Shield<SEXP>(Rcpp::wrap(MMM)), Shield<SEXP>(Rcpp::wrap(Yy)), Shield<SEXP>(Rcpp::wrap(Ttt)), Shield<SEXP>(Rcpp::wrap(BETAIN)), Shield<SEXP>(Rcpp::wrap(BETAOUT)), Shield<SEXP>(Rcpp::wrap(WW)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline SEXP c_postzeroprob(SEXP Xi1, SEXP Xi2, SEXP Xj1, SEXP Xj2, SEXP SS2, SEXP LAM, SEXP PP0) {
        typedef SEXP(*Ptr_c_postzeroprob)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_c_postzeroprob p_c_postzeroprob = NULL;
        if (p_c_postzeroprob == NULL) {
            validateSignature("SEXP(*c_postzeroprob)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP)");
            p_c_postzeroprob = (Ptr_c_postzeroprob)R_GetCCallable("DynMultiNet", "_DynMultiNet_c_postzeroprob");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_c_postzeroprob(Shield<SEXP>(Rcpp::wrap(Xi1)), Shield<SEXP>(Rcpp::wrap(Xi2)), Shield<SEXP>(Rcpp::wrap(Xj1)), Shield<SEXP>(Rcpp::wrap(Xj2)), Shield<SEXP>(Rcpp::wrap(SS2)), Shield<SEXP>(Rcpp::wrap(LAM)), Shield<SEXP>(Rcpp::wrap(PP0)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

    inline SEXP c_prediction(SEXP EX, SEXP SIG2, SEXP X1T, SEXP X2T, SEXP BIN, SEXP BOUT, SEXP WW) {
        typedef SEXP(*Ptr_c_prediction)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_c_prediction p_c_prediction = NULL;
        if (p_c_prediction == NULL) {
            validateSignature("SEXP(*c_prediction)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP)");
            p_c_prediction = (Ptr_c_prediction)R_GetCCallable("DynMultiNet", "_DynMultiNet_c_prediction");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_c_prediction(Shield<SEXP>(Rcpp::wrap(EX)), Shield<SEXP>(Rcpp::wrap(SIG2)), Shield<SEXP>(Rcpp::wrap(X1T)), Shield<SEXP>(Rcpp::wrap(X2T)), Shield<SEXP>(Rcpp::wrap(BIN)), Shield<SEXP>(Rcpp::wrap(BOUT)), Shield<SEXP>(Rcpp::wrap(WW)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<SEXP >(rcpp_result_gen);
    }

}

#endif // RCPP_DynMultiNet_RCPPEXPORTS_H_GEN_
