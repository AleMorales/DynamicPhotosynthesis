// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_RcppSundials_RCPPEXPORTS_H_GEN_
#define RCPP_RcppSundials_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace RcppSundials {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("RcppSundials", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("RcppSundials", "_RcppSundials_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in RcppSundials");
            }
        }
    }

    inline NumericVector wrap_cvodes(NumericVector times, NumericVector states_, NumericVector parameters_, List forcings_data_, List settings, SEXP model_, SEXP jacobian_) {
        typedef SEXP(*Ptr_wrap_cvodes)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_wrap_cvodes p_wrap_cvodes = NULL;
        if (p_wrap_cvodes == NULL) {
            validateSignature("NumericVector(*wrap_cvodes)(NumericVector,NumericVector,NumericVector,List,List,SEXP,SEXP)");
            p_wrap_cvodes = (Ptr_wrap_cvodes)R_GetCCallable("RcppSundials", "_RcppSundials_wrap_cvodes");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_wrap_cvodes(Shield<SEXP>(Rcpp::wrap(times)), Shield<SEXP>(Rcpp::wrap(states_)), Shield<SEXP>(Rcpp::wrap(parameters_)), Shield<SEXP>(Rcpp::wrap(forcings_data_)), Shield<SEXP>(Rcpp::wrap(settings)), Shield<SEXP>(Rcpp::wrap(model_)), Shield<SEXP>(Rcpp::wrap(jacobian_)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

    inline List cvode_calc_derivs(SEXP model_, NumericVector t, NumericVector states, NumericVector parameters, List forcings_data_) {
        typedef SEXP(*Ptr_cvode_calc_derivs)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvode_calc_derivs p_cvode_calc_derivs = NULL;
        if (p_cvode_calc_derivs == NULL) {
            validateSignature("List(*cvode_calc_derivs)(SEXP,NumericVector,NumericVector,NumericVector,List)");
            p_cvode_calc_derivs = (Ptr_cvode_calc_derivs)R_GetCCallable("RcppSundials", "_RcppSundials_cvode_calc_derivs");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cvode_calc_derivs(Shield<SEXP>(Rcpp::wrap(model_)), Shield<SEXP>(Rcpp::wrap(t)), Shield<SEXP>(Rcpp::wrap(states)), Shield<SEXP>(Rcpp::wrap(parameters)), Shield<SEXP>(Rcpp::wrap(forcings_data_)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline NumericMatrix cvode_calc_jac(SEXP jacobian_, NumericVector t, NumericVector states, NumericVector parameters, List forcings_data_) {
        typedef SEXP(*Ptr_cvode_calc_jac)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cvode_calc_jac p_cvode_calc_jac = NULL;
        if (p_cvode_calc_jac == NULL) {
            validateSignature("NumericMatrix(*cvode_calc_jac)(SEXP,NumericVector,NumericVector,NumericVector,List)");
            p_cvode_calc_jac = (Ptr_cvode_calc_jac)R_GetCCallable("RcppSundials", "_RcppSundials_cvode_calc_jac");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cvode_calc_jac(Shield<SEXP>(Rcpp::wrap(jacobian_)), Shield<SEXP>(Rcpp::wrap(t)), Shield<SEXP>(Rcpp::wrap(states)), Shield<SEXP>(Rcpp::wrap(parameters)), Shield<SEXP>(Rcpp::wrap(forcings_data_)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericMatrix >(rcpp_result_gen);
    }

}

#endif // RCPP_RcppSundials_RCPPEXPORTS_H_GEN_
