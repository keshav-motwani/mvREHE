#' compute_w_diag
#'
#' LDSC-style weights for the same-SNP (diagonal) moment equations.
#' The weight for SNP i is 1 / (l2_i * Var(z_a z_b)_i), where l2_i (to the
#' first power) corrects for the ~l2_i near-duplicate equations of SNP i's LD
#' partners, and the variance term uses plug-in values shared across trait
#' pairs. Note the cross-SNP weights in compute_W_sparse use the product
#' l2_r * l2_c instead: pair equation (r, c) has ~l2_r * l2_c near-duplicates.
#'
#' @param R_diag  diagonal of R (typically all 1)
#' @param l2      LD scores, diag(R^2)
#' @param h2      plug-in heritability
#' @param rho_g   plug-in genetic correlation
#' @param rho_e   plug-in residual correlation
#' @param N       GWAS sample size
#' @param M       number of reference SNPs
#'
#' @return numeric vector of per-SNP weights
compute_w_diag = function(R_diag, l2, h2, rho_g, rho_e, N, M) {
  term_s  = (N / M) * l2 * h2 + R_diag * (1 - h2)
  term_st = (N / M) * l2 * rho_g * h2 + R_diag * rho_e * (1 - h2)
  # Floor the LD weight at R_diag as in ldsc's max(l2, 1): keeps SNPs with
  # l2 = 0 (possible when l2 = diag(RsRsT)) at finite weight.
  1 / (pmax(l2, R_diag) * (term_s^2 + term_st^2))
}

#' mvHE_SNP
#'
#' @param Y
#' @param D_list
#' @param GWAS_N
#' @param M
#' @param truncate
#'
#' @return
#' @export
mvHE_SNP = function(Y, D_list, GWAS_N, M = NULL, init_h2 = 0.1, init_rho_g = 0.3, init_rho_e = 0.05, truncate = TRUE, cross_snp = FALSE) {

  stopifnot(Matrix::nnzero(D_list[[2]]) >= Matrix::nnzero(D_list[[1]]))
  stopifnot(all(sapply(D_list, function(x) is(x, "dsCMatrix"))))

  if (is.null(M)) M = nrow(D_list[[1]])

  R  = D_list[[1]]
  R2 = D_list[[2]] * M / GWAS_N

  col_vars  = matrixStats::colVars(Y)
  col_vars[col_vars < 1e-10] = 1
  w_columns = 1 / col_vars

  if (cross_snp) {
    W_row_pairs = compute_W_sparse(R, R2, Matrix::diag(R), Matrix::diag(R2),
                                    h2 = init_h2, rho_g = init_rho_g, rho_e = init_rho_e,
                                    N = GWAS_N, M = M)
    fit = mvHE(Y, D_list, W_row_pairs = W_row_pairs, w_columns = w_columns,
               truncate = truncate)
  } else {
    w = compute_w_diag(Matrix::diag(R), Matrix::diag(R2),
                       h2 = init_h2, rho_g = init_rho_g, rho_e = init_rho_e,
                       N = GWAS_N, M = M)
    fit = mvHE_diag(Y, D_list, W_row_pairs = w, w_columns = w_columns,
                    truncate = truncate)
  }
  return(fit)
}

#' mvREHE_SNP
#'
#' @param Y
#' @param D_list
#' @param GWAS_N
#' @param M
#' @param init
#' @param truncate
#'
#' @return
#' @export
mvREHE_SNP = function(Y, D_list, GWAS_N, M = NULL, init_h2 = 0.1, init_rho_g = 0.3, init_rho_e = 0.05, truncate = TRUE, cross_snp = FALSE) {

  stopifnot(Matrix::nnzero(D_list[[2]]) >= Matrix::nnzero(D_list[[1]]))
  stopifnot(all(sapply(D_list, function(x) is(x, "dsCMatrix"))))

  if (is.null(M)) M = nrow(D_list[[1]])

  R  = D_list[[1]]
  R2 = D_list[[2]] * M / GWAS_N

  col_vars  = matrixStats::colVars(Y)
  col_vars[col_vars < 1e-10] = 1
  w_columns = 1 / col_vars

  if (cross_snp) {
    W_row_pairs = compute_W_sparse(R, R2, Matrix::diag(R), Matrix::diag(R2),
                                    h2 = init_h2, rho_g = init_rho_g, rho_e = init_rho_e,
                                    N = GWAS_N, M = M)
    fit = mvREHE(Y, D_list, W_row_pairs = W_row_pairs, w_columns = w_columns,
                 truncate = truncate)
  } else {
    w = compute_w_diag(Matrix::diag(R), Matrix::diag(R2),
                       h2 = init_h2, rho_g = init_rho_g, rho_e = init_rho_e,
                       N = GWAS_N, M = M)
    fit = mvREHE_diag(Y, D_list, W_row_pairs = w, w_columns = w_columns,
                      truncate = truncate)
  }
  return(fit)

}

