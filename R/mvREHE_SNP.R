#' mvREHE_SNP
#'
#' @param Y
#' @param D_list
#' @param tolerance
#' @param max_iter
#' @param Sigma_init_list
#' @param W_list
#' @param Q
#' @param row_indices
#' @param col_indices
#'
#' @return
#' @export
#'
#' @examples
mvREHE_SNP = function(Y, D_list, GWAS_N, M = NULL, init_h2 = 0.1, init_rho_g = 0.3, init_rho_e = 0.05, refit = TRUE, tolerance = 1e-6, max_iter = 1000, return_full = TRUE, Sigma_init_list = NULL) {

  stopifnot(Matrix::nnzero(D_list[[2]]) > Matrix::nnzero(D_list[[1]]))
  stopifnot(all(sapply(D_list, function(x) is(x, "dsCMatrix"))))

  if (is.null(M)) M = nrow(D_list[[1]])

  R = D_list[[1]]
  R2 = D_list[[2]] * M / GWAS_N

  W_row_pairs = compute_W_sparse(R, R2, Matrix::diag(R), Matrix::diag(R2), h2 = init_h2, rho_g = init_rho_g, rho_e = init_rho_e, N = GWAS_N, M = M)
  col_vars = matrixStats::colVars(Y)
  w_columns = 1 / col_vars
  w_columns[col_vars < 1e-10] = 1
  fit = mvREHE(Y, D_list, W_row_pairs = W_row_pairs, w_columns = w_columns, tolerance, max_iter, return_full, Sigma_init_list)

  if (refit) {

    avg_h2 = mean(diag(fit$Sigma_hat[[2]]) / (diag(fit$Sigma_hat[[1]]) + diag(fit$Sigma_hat[[2]])), na.rm = TRUE)

    avg_rho_g = mean(cov2cor(fit$Sigma_hat[[2]])[lower.tri(fit$Sigma_hat[[1]])], na.rm = TRUE)
    avg_rho_e = mean(cov2cor(fit$Sigma_hat[[1]])[lower.tri(fit$Sigma_hat[[2]])], na.rm = TRUE)

    print(c(avg_h2, avg_rho_g, avg_rho_e))

    W_row_pairs = compute_W_sparse(R, R2, Matrix::diag(R), Matrix::diag(R2), h2 = avg_h2, rho_g = avg_rho_g, rho_e = avg_rho_e, N = GWAS_N, M = M)
    total_vars = diag(fit$Sigma_hat[[1]] + fit$Sigma_hat[[2]])
    w_columns = 1 / total_vars
    w_columns[total_vars < 1e-10] = 1
    fit = mvREHE(Y, D_list, W_row_pairs = W_row_pairs, w_columns = w_columns, tolerance, max_iter, return_full, Sigma_init_list)

  }

  return(fit)

}
