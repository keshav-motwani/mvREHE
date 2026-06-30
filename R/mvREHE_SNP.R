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
#'
#' @examples
mvHE_SNP = function(Y, D_list, GWAS_N, M = NULL, truncate = TRUE) {

  stopifnot(Matrix::nnzero(D_list[[2]]) > Matrix::nnzero(D_list[[1]]))
  stopifnot(all(sapply(D_list, function(x) is(x, "dsCMatrix"))))

  if (is.null(M)) M = nrow(D_list[[1]])

  R  = D_list[[1]]
  R2 = D_list[[2]] * M / GWAS_N  # recover unscaled R2_overlap for compute_W_sparse

  # Step 1: single-shot unweighted OLS for hyperparameters
  fit_uw = mvHE(Y, D_list, truncate = TRUE)
  Sg_psd = fit_uw$Sigma_hat[[2]]
  Se_psd = fit_uw$Sigma_hat[[1]]
  Sg_pd  = Sg_psd + diag(1e-10, nrow(Sg_psd))
  Se_pd  = Se_psd + diag(1e-10, nrow(Se_psd))

  # h2: D_list[[2]] is pre-scaled by N/M, so Sg is already on the individual phenotype scale
  Sg_ind  = diag(Sg_psd)
  avg_h2  = mean(Sg_ind / pmax(diag(Se_psd) + Sg_ind, 1e-10), na.rm = TRUE)

  # rho_g, rho_e
  avg_rho_g = mean(cov2cor(Sg_pd)[lower.tri(Sg_pd, diag = TRUE)], na.rm = TRUE)
  avg_rho_e = mean(cov2cor(Se_pd)[lower.tri(Se_pd, diag = TRUE)], na.rm = TRUE)

  # Step 2: compute weights using data-driven hyperparameters
  W_row_pairs = compute_W_sparse(R, R2, Matrix::diag(R), Matrix::diag(R2),
                                  h2 = avg_h2, rho_g = avg_rho_g, rho_e = avg_rho_e,
                                  N = GWAS_N, M = M)
  col_vars  = matrixStats::colVars(Y)
  col_vars[col_vars < 1e-10] = 1
  w_columns = 1 / col_vars

  # Step 3: weighted single-shot OLS
  fit = mvHE(Y, D_list, W_row_pairs = W_row_pairs, w_columns = w_columns,
             truncate = truncate)

  return(fit)

}

#' mvREHE_SNP
#'
#' @param Y
#' @param D_list
#' @param GWAS_N
#' @param M
#' @param init "mvHE" (single-shot OLS) or "mvREHE" (iterative) for the unweighted
#'   first pass that derives weighting hyperparameters.
#'   D_list[[2]] must be R2_overlap * GWAS_N / M (pre-scaled by N/M).
#' @param truncate
#'
#' @return
#' @export
#'
#' @examples
mvREHE_SNP = function(Y, D_list, GWAS_N, M = NULL, init = c("mvHE", "mvREHE"), truncate = TRUE) {

  init = match.arg(init)

  stopifnot(Matrix::nnzero(D_list[[2]]) > Matrix::nnzero(D_list[[1]]))
  stopifnot(all(sapply(D_list, function(x) is(x, "dsCMatrix"))))

  if (is.null(M)) M = nrow(D_list[[1]])

  R  = D_list[[1]]
  R2 = D_list[[2]] * M / GWAS_N  # recover unscaled R2_overlap for compute_W_sparse

  # Step 1: unweighted fit to derive data-driven hyperparameters for W
  if (init == "mvHE") {
    fit_uw = mvHE(Y, D_list, truncate = TRUE)
  } else {
    fit_uw = mvREHE(Y, D_list, truncate = truncate)
  }
  Sg_psd = fit_uw$Sigma_hat[[2]]
  Se_psd = fit_uw$Sigma_hat[[1]]
  Sg_pd  = Sg_psd + diag(1e-10, nrow(Sg_psd))
  Se_pd  = Se_psd + diag(1e-10, nrow(Se_psd))

  # h2: D_list[[2]] is pre-scaled by N/M, so Sg is already on the individual phenotype scale
  Sg_ind  = diag(Sg_psd)
  avg_h2  = mean(Sg_ind / pmax(diag(Se_psd) + Sg_ind, 1e-10), na.rm = TRUE)

  # rho_g, rho_e
  avg_rho_g = mean(cov2cor(Sg_pd)[lower.tri(Sg_pd, diag = TRUE)], na.rm = TRUE)
  avg_rho_e = mean(cov2cor(Se_pd)[lower.tri(Se_pd, diag = TRUE)], na.rm = TRUE)

  # Step 2: compute weights using data-driven hyperparameters
  W_row_pairs = compute_W_sparse(R, R2, Matrix::diag(R), Matrix::diag(R2),
                                  h2 = avg_h2, rho_g = avg_rho_g, rho_e = avg_rho_e,
                                  N = GWAS_N, M = M)
  col_vars  = matrixStats::colVars(Y)
  col_vars[col_vars < 1e-10] = 1
  w_columns = 1 / col_vars

  # Step 3: weighted iterative fit
  fit = mvREHE(Y, D_list, W_row_pairs = W_row_pairs, w_columns = w_columns,
               truncate = truncate)

  return(fit)

}
