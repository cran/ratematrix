# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

makeSimmapMappedEdge <- function(n_nodes, n_tips, n_states, edge_len, edge_mat, parents, X, Q, root_node, root_type, sims_limit) {
    .Call('_ratematrix_makeSimmapMappedEdge', PACKAGE = 'ratematrix', n_nodes, n_tips, n_states, edge_len, edge_mat, parents, X, Q, root_node, root_type, sims_limit)
}

makeSimmapMaps <- function(n_nodes, n_tips, n_states, edge_len, edge_mat, parents, X, Q, root_node, root_type, max_nshifts) {
    .Call('_ratematrix_makeSimmapMaps', PACKAGE = 'ratematrix', n_nodes, n_tips, n_states, edge_len, edge_mat, parents, X, Q, root_node, root_type, max_nshifts)
}

logLikMk_C <- function(n_nodes, n_tips, n_states, edge_len, edge_mat, parents, X, Q, root_node, root_type) {
    .Call('_ratematrix_logLikMk_C', PACKAGE = 'ratematrix', n_nodes, n_tips, n_states, edge_len, edge_mat, parents, X, Q, root_node, root_type)
}

logLikPrunningMCMC_C <- function(X, k, p, nodes, des, anc, names_anc, mapped_edge, R, mu) {
    .Call('_ratematrix_logLikPrunningMCMC_C', PACKAGE = 'ratematrix', X, k, p, nodes, des, anc, names_anc, mapped_edge, R, mu)
}

runRatematrixMCMC_C <- function(X, k, p, nodes, des, anc, names_anc, mapped_edge, R, mu, sd, Rcorr, w_mu, par_prior_mu, den_mu, w_sd, par_prior_sd, den_sd, nu, sigma, v, log_file, mcmc_file, prob_sample_root, prob_sample_sd, gen, post_seq, write_header) {
    .Call('_ratematrix_runRatematrixMCMC_C', PACKAGE = 'ratematrix', X, k, p, nodes, des, anc, names_anc, mapped_edge, R, mu, sd, Rcorr, w_mu, par_prior_mu, den_mu, w_sd, par_prior_sd, den_sd, nu, sigma, v, log_file, mcmc_file, prob_sample_root, prob_sample_sd, gen, post_seq, write_header)
}

runRatematrixMultiMCMC_C <- function(X, k, p, nodes, des, anc, names_anc, mapped_edge, R, mu, sd, Rcorr, w_mu, par_prior_mu, den_mu, w_sd, par_prior_sd, den_sd, nu, sigma, v, log_file, mcmc_file, prob_sample_root, prob_sample_sd, gen, post_seq, write_header) {
    .Call('_ratematrix_runRatematrixMultiMCMC_C', PACKAGE = 'ratematrix', X, k, p, nodes, des, anc, names_anc, mapped_edge, R, mu, sd, Rcorr, w_mu, par_prior_mu, den_mu, w_sd, par_prior_sd, den_sd, nu, sigma, v, log_file, mcmc_file, prob_sample_root, prob_sample_sd, gen, post_seq, write_header)
}

buildQ <- function(vec_Q, size, model_Q) {
    .Call('_ratematrix_buildQ', PACKAGE = 'ratematrix', vec_Q, size, model_Q)
}

runRatematrixMCMC_jointMk_C <- function(X, datMk, k, p, nodes, n_tips, des, anc, names_anc, mapped_edge, edge_mat, n_nodes, Q, w_Q, model_Q, root_type, den_Q, par_prior_Q, R, mu, sd, Rcorr, w_mu, par_prior_mu, den_mu, w_sd, par_prior_sd, den_sd, nu, sigma, v, log_file, mcmc_file, Q_mcmc_file, par_prob, gen, post_seq, write_header, sims_limit) {
    .Call('_ratematrix_runRatematrixMCMC_jointMk_C', PACKAGE = 'ratematrix', X, datMk, k, p, nodes, n_tips, des, anc, names_anc, mapped_edge, edge_mat, n_nodes, Q, w_Q, model_Q, root_type, den_Q, par_prior_Q, R, mu, sd, Rcorr, w_mu, par_prior_mu, den_mu, w_sd, par_prior_sd, den_sd, nu, sigma, v, log_file, mcmc_file, Q_mcmc_file, par_prob, gen, post_seq, write_header, sims_limit)
}

