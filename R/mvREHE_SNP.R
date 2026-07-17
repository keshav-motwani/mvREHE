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
mvHE_SNP = function(Y, D_list, GWAS_N, M = NULL, truncate = TRUE, cross_snp = FALSE) {

  stopifnot(Matrix::nnzero(D_list[[2]]) >= Matrix::nnzero(D_list[[1]]))
  stopifnot(all(sapply(D_list, function(x) is(x, "dsCMatrix"))))

  if (is.null(M)) M = nrow(D_list[[1]])

  R  = D_list[[1]]
  R2 = D_list[[2]] * M / GWAS_N

  fit_uw = mvHE_diag(Y, D_list, truncate = TRUE)
  Sg_psd = fit_uw$Sigma_hat[[2]]
  Se_psd = fit_uw$Sigma_hat[[1]]
  Sg_pd  = Sg_psd + diag(1e-10, nrow(Sg_psd))
  Se_pd  = Se_psd + diag(1e-10, nrow(Se_psd))
  Sg_ind    = diag(Sg_psd)
  avg_h2    = mean(Sg_ind / pmax(diag(Se_psd) + Sg_ind, 1e-10), na.rm = TRUE)
  avg_rho_g = mean(cov2cor(Sg_pd)[lower.tri(Sg_pd, diag = TRUE)], na.rm = TRUE)
  avg_rho_e = mean(cov2cor(Se_pd)[lower.tri(Se_pd, diag = TRUE)], na.rm = TRUE)

  W_row_pairs = compute_W_sparse(R, R2, Matrix::diag(R), Matrix::diag(R2),
                                  h2 = avg_h2, rho_g = avg_rho_g, rho_e = avg_rho_e,
                                  N = GWAS_N, M = M)
  col_vars  = matrixStats::colVars(Y)
  col_vars[col_vars < 1e-10] = 1
  w_columns = 1 / col_vars

  if (cross_snp) {
    fit = mvHE(Y, D_list, W_row_pairs = W_row_pairs, w_columns = w_columns,
               truncate = truncate)
  } else {
    fit = mvHE_diag(Y, D_list, W_row_pairs = W_row_pairs, w_columns = w_columns,
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
mvREHE_SNP = function(Y, D_list, GWAS_N, M = NULL, init = c("mvHE", "mvREHE"), truncate = TRUE, cross_snp = FALSE) {

  init = match.arg(init)

  stopifnot(Matrix::nnzero(D_list[[2]]) >= Matrix::nnzero(D_list[[1]]))
  stopifnot(all(sapply(D_list, function(x) is(x, "dsCMatrix"))))

  if (is.null(M)) M = nrow(D_list[[1]])

  R  = D_list[[1]]
  R2 = D_list[[2]] * M / GWAS_N

  if (init == "mvHE") {
    fit_uw = mvHE_diag(Y, D_list, truncate = TRUE)
  } else {
    fit_uw = mvREHE_diag(Y, D_list, truncate = truncate)
  }
  Sg_psd = fit_uw$Sigma_hat[[2]]
  Se_psd = fit_uw$Sigma_hat[[1]]
  Sg_pd  = Sg_psd + diag(1e-10, nrow(Sg_psd))
  Se_pd  = Se_psd + diag(1e-10, nrow(Se_psd))
  Sg_ind    = diag(Sg_psd)
  avg_h2    = mean(Sg_ind / pmax(diag(Se_psd) + Sg_ind, 1e-10), na.rm = TRUE)
  avg_rho_g = mean(cov2cor(Sg_pd)[lower.tri(Sg_pd, diag = TRUE)], na.rm = TRUE)
  avg_rho_e = mean(cov2cor(Se_pd)[lower.tri(Se_pd, diag = TRUE)], na.rm = TRUE)

  W_row_pairs = compute_W_sparse(R, R2, Matrix::diag(R), Matrix::diag(R2),
                                  h2 = avg_h2, rho_g = avg_rho_g, rho_e = avg_rho_e,
                                  N = GWAS_N, M = M)
  col_vars  = matrixStats::colVars(Y)
  col_vars[col_vars < 1e-10] = 1
  w_columns = 1 / col_vars

  if (cross_snp) {
    fit = mvREHE(Y, D_list, W_row_pairs = W_row_pairs, w_columns = w_columns,
                 truncate = truncate)
  } else {
    fit = mvREHE_diag(Y, D_list, W_row_pairs = W_row_pairs, w_columns = w_columns,
                      truncate = truncate)
  }
  return(fit)
}

