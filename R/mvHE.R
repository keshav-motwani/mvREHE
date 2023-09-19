#' mvHE
#'
#' @param Y
#' @param D_list
#'
#' @return
#' @export
#'
#' @examples
mvHE = function(Y, D_list) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  q = ncol(Y)

  Sigma_hat = replicate(length(D_list), matrix(NA, q, q), simplify = FALSE)

  indices = do.call(`+`, D_list) > 0
  row_indices = which(indices, arr.ind = TRUE)[, 1]
  col_indices = which(indices, arr.ind = TRUE)[, 2]
  indices = which(indices)

  X_tilde = do.call(cbind, lapply(D_list, c))
  X_tilde = X_tilde[indices, ]

  XtXinv = solve(crossprod(X_tilde))

  for (j in 1:q) {

    for (m in 1:j) {

      Y_tilde = compute_Y_tilde(Y, row_indices - 1, col_indices - 1, j - 1, m - 1)

      sigma_hat = XtXinv %*% crossprod(X_tilde, Y_tilde)
      for (k in 1:length(D_list)) {
        Sigma_hat[[k]][j, m] = Sigma_hat[[k]][m, j] = sigma_hat[k]
      }

    }

  }

  for (k in 1:length(D_list)) {

    Sigma_k_hat = Sigma_hat[[k]]
    eigen_Sigma_k_hat = eigen(Sigma_k_hat)
    if (any(eigen_Sigma_k_hat$values < 0)) {
      Sigma_hat[[k]] = eigen_Sigma_k_hat$vectors %*% diag(pmax(eigen_Sigma_k_hat$values, 0)) %*% t(eigen_Sigma_k_hat$vectors)
      attr(Sigma_hat[[k]], "truncated") = TRUE
    } else {
      attr(Sigma_hat[[k]], "truncated") = FALSE
    }
    attr(Sigma_hat[[k]], "min_eigenvalue") = min(eigen_Sigma_k_hat$values)

  }

  return(list(Sigma_hat = Sigma_hat))

}
