// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Sim_SP
SEXP Sim_SP(SEXP x_, SEXP cstar_, SEXP model_, SEXP corrFn_, SEXP keepPsi_, SEXP covPar_, SEXP alpha_);
RcppExport SEXP _MSP_Sim_SP(SEXP x_SEXP, SEXP cstar_SEXP, SEXP model_SEXP, SEXP corrFn_SEXP, SEXP keepPsi_SEXP, SEXP covPar_SEXP, SEXP alpha_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x_(x_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type cstar_(cstar_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type corrFn_(corrFn_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type keepPsi_(keepPsi_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type covPar_(covPar_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type alpha_(alpha_SEXP);
    rcpp_result_gen = Rcpp::wrap(Sim_SP(x_, cstar_, model_, corrFn_, keepPsi_, covPar_, alpha_));
    return rcpp_result_gen;
END_RCPP
}
// Sim_Reich
SEXP Sim_Reich(SEXP x_, SEXP knots_, SEXP Sigma_, SEXP keepPsi_, SEXP bw_, SEXP alpha_);
RcppExport SEXP _MSP_Sim_Reich(SEXP x_SEXP, SEXP knots_SEXP, SEXP Sigma_SEXP, SEXP keepPsi_SEXP, SEXP bw_SEXP, SEXP alpha_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x_(x_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type knots_(knots_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Sigma_(Sigma_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type keepPsi_(keepPsi_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type bw_(bw_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type alpha_(alpha_SEXP);
    rcpp_result_gen = Rcpp::wrap(Sim_Reich(x_, knots_, Sigma_, keepPsi_, bw_, alpha_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MSP_Sim_SP", (DL_FUNC) &_MSP_Sim_SP, 7},
    {"_MSP_Sim_Reich", (DL_FUNC) &_MSP_Sim_Reich, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_MSP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}