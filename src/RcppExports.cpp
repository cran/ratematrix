// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logLikPrunningMCMC_C
double logLikPrunningMCMC_C(arma::mat X, int k, int p, arma::uvec nodes, arma::uvec des, arma::uvec anc, arma::uvec names_anc, arma::mat mapped_edge, arma::cube R, arma::vec mu);
RcppExport SEXP _ratematrix_logLikPrunningMCMC_C(SEXP XSEXP, SEXP kSEXP, SEXP pSEXP, SEXP nodesSEXP, SEXP desSEXP, SEXP ancSEXP, SEXP names_ancSEXP, SEXP mapped_edgeSEXP, SEXP RSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type des(desSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type anc(ancSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type names_anc(names_ancSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mapped_edge(mapped_edgeSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikPrunningMCMC_C(X, k, p, nodes, des, anc, names_anc, mapped_edge, R, mu));
    return rcpp_result_gen;
END_RCPP
}
// cov2cor_C
arma::mat cov2cor_C(arma::mat V);
RcppExport SEXP _ratematrix_cov2cor_C(SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(cov2cor_C(V));
    return rcpp_result_gen;
END_RCPP
}
// priorCorr_C
double priorCorr_C(arma::cube corr, arma::vec nu, arma::cube sigma);
RcppExport SEXP _ratematrix_priorCorr_C(SEXP corrSEXP, SEXP nuSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type corr(corrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(priorCorr_C(corr, nu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// makePropIWish_C
arma::mat makePropIWish_C(arma::mat vcv, double k, double v);
RcppExport SEXP _ratematrix_makePropIWish_C(SEXP vcvSEXP, SEXP kSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type vcv(vcvSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(makePropIWish_C(vcv, k, v));
    return rcpp_result_gen;
END_RCPP
}
// runRatematrixMCMC_C
std::string runRatematrixMCMC_C(arma::mat X, int k, int p, arma::uvec nodes, arma::uvec des, arma::uvec anc, arma::uvec names_anc, arma::mat mapped_edge, arma::cube R, arma::vec mu, arma::mat sd, arma::cube Rcorr, arma::vec w_mu, arma::mat par_prior_mu, std::string den_mu, arma::mat w_sd, arma::mat par_prior_sd, std::string den_sd, arma::vec nu, arma::cube sigma, arma::vec v, std::string log_file, std::string mcmc_file, double prob_sample_root, double prob_sample_sd, int gen, int write_header);
RcppExport SEXP _ratematrix_runRatematrixMCMC_C(SEXP XSEXP, SEXP kSEXP, SEXP pSEXP, SEXP nodesSEXP, SEXP desSEXP, SEXP ancSEXP, SEXP names_ancSEXP, SEXP mapped_edgeSEXP, SEXP RSEXP, SEXP muSEXP, SEXP sdSEXP, SEXP RcorrSEXP, SEXP w_muSEXP, SEXP par_prior_muSEXP, SEXP den_muSEXP, SEXP w_sdSEXP, SEXP par_prior_sdSEXP, SEXP den_sdSEXP, SEXP nuSEXP, SEXP sigmaSEXP, SEXP vSEXP, SEXP log_fileSEXP, SEXP mcmc_fileSEXP, SEXP prob_sample_rootSEXP, SEXP prob_sample_sdSEXP, SEXP genSEXP, SEXP write_headerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type des(desSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type anc(ancSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type names_anc(names_ancSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mapped_edge(mapped_edgeSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Rcorr(RcorrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w_mu(w_muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type par_prior_mu(par_prior_muSEXP);
    Rcpp::traits::input_parameter< std::string >::type den_mu(den_muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w_sd(w_sdSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type par_prior_sd(par_prior_sdSEXP);
    Rcpp::traits::input_parameter< std::string >::type den_sd(den_sdSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< std::string >::type log_file(log_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type mcmc_file(mcmc_fileSEXP);
    Rcpp::traits::input_parameter< double >::type prob_sample_root(prob_sample_rootSEXP);
    Rcpp::traits::input_parameter< double >::type prob_sample_sd(prob_sample_sdSEXP);
    Rcpp::traits::input_parameter< int >::type gen(genSEXP);
    Rcpp::traits::input_parameter< int >::type write_header(write_headerSEXP);
    rcpp_result_gen = Rcpp::wrap(runRatematrixMCMC_C(X, k, p, nodes, des, anc, names_anc, mapped_edge, R, mu, sd, Rcorr, w_mu, par_prior_mu, den_mu, w_sd, par_prior_sd, den_sd, nu, sigma, v, log_file, mcmc_file, prob_sample_root, prob_sample_sd, gen, write_header));
    return rcpp_result_gen;
END_RCPP
}
// runRatematrixMultiMCMC_C
std::string runRatematrixMultiMCMC_C(arma::mat X, int k, int p, arma::umat nodes, arma::umat des, arma::umat anc, arma::umat names_anc, arma::cube mapped_edge, arma::cube R, arma::vec mu, arma::mat sd, arma::cube Rcorr, arma::vec w_mu, arma::mat par_prior_mu, std::string den_mu, arma::mat w_sd, arma::mat par_prior_sd, std::string den_sd, arma::vec nu, arma::cube sigma, arma::vec v, std::string log_file, std::string mcmc_file, double prob_sample_root, double prob_sample_sd, int gen, int write_header);
RcppExport SEXP _ratematrix_runRatematrixMultiMCMC_C(SEXP XSEXP, SEXP kSEXP, SEXP pSEXP, SEXP nodesSEXP, SEXP desSEXP, SEXP ancSEXP, SEXP names_ancSEXP, SEXP mapped_edgeSEXP, SEXP RSEXP, SEXP muSEXP, SEXP sdSEXP, SEXP RcorrSEXP, SEXP w_muSEXP, SEXP par_prior_muSEXP, SEXP den_muSEXP, SEXP w_sdSEXP, SEXP par_prior_sdSEXP, SEXP den_sdSEXP, SEXP nuSEXP, SEXP sigmaSEXP, SEXP vSEXP, SEXP log_fileSEXP, SEXP mcmc_fileSEXP, SEXP prob_sample_rootSEXP, SEXP prob_sample_sdSEXP, SEXP genSEXP, SEXP write_headerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type des(desSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type anc(ancSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type names_anc(names_ancSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type mapped_edge(mapped_edgeSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Rcorr(RcorrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w_mu(w_muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type par_prior_mu(par_prior_muSEXP);
    Rcpp::traits::input_parameter< std::string >::type den_mu(den_muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type w_sd(w_sdSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type par_prior_sd(par_prior_sdSEXP);
    Rcpp::traits::input_parameter< std::string >::type den_sd(den_sdSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< std::string >::type log_file(log_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type mcmc_file(mcmc_fileSEXP);
    Rcpp::traits::input_parameter< double >::type prob_sample_root(prob_sample_rootSEXP);
    Rcpp::traits::input_parameter< double >::type prob_sample_sd(prob_sample_sdSEXP);
    Rcpp::traits::input_parameter< int >::type gen(genSEXP);
    Rcpp::traits::input_parameter< int >::type write_header(write_headerSEXP);
    rcpp_result_gen = Rcpp::wrap(runRatematrixMultiMCMC_C(X, k, p, nodes, des, anc, names_anc, mapped_edge, R, mu, sd, Rcorr, w_mu, par_prior_mu, den_mu, w_sd, par_prior_sd, den_sd, nu, sigma, v, log_file, mcmc_file, prob_sample_root, prob_sample_sd, gen, write_header));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ratematrix_logLikPrunningMCMC_C", (DL_FUNC) &_ratematrix_logLikPrunningMCMC_C, 10},
    {"_ratematrix_cov2cor_C", (DL_FUNC) &_ratematrix_cov2cor_C, 1},
    {"_ratematrix_priorCorr_C", (DL_FUNC) &_ratematrix_priorCorr_C, 3},
    {"_ratematrix_makePropIWish_C", (DL_FUNC) &_ratematrix_makePropIWish_C, 3},
    {"_ratematrix_runRatematrixMCMC_C", (DL_FUNC) &_ratematrix_runRatematrixMCMC_C, 27},
    {"_ratematrix_runRatematrixMultiMCMC_C", (DL_FUNC) &_ratematrix_runRatematrixMultiMCMC_C, 27},
    {NULL, NULL, 0}
};

RcppExport void R_init_ratematrix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}