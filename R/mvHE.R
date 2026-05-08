#' mvHE
#'
#' @param Y
#' @param D_list
#'
#' @return
#' @export
#'
#' @examples
mvHE = function(Y, D_list, truncate = TRUE) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  q = ncol(Y)
  n = nrow(Y)

  highdim = q > n

  if (highdim) {

    s = svd(Y)
    Y = Y %*% s$v
    q = ncol(Y)

  }

  Sigma_hat = replicate(length(D_list), matrix(NA, q, q), simplify = FALSE)

  if (all(sapply(D_list, function(x) is(x, "sparseMatrix")))) {

    sum_D = Reduce(`+`, D_list)
    sum_T = as(sum_D, "dgTMatrix")

    row_indices = sum_T@i + 1
    col_indices = sum_T@j + 1

    X_tilde_list = lapply(D_list, function(m) {
      m[cbind(row_indices, col_indices)]
    })

    X_tilde = do.call(cbind, X_tilde_list)

  } else {
    mask = Reduce(`+`, D_list) > 0
    coord_matrix = which(mask, arr.ind = TRUE)
    row_indices = coord_matrix[, 1]
    col_indices = coord_matrix[, 2]
    indices = which(mask)
    X_tilde = sapply(D_list, function(m) m[indices])
  }

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

  if (truncate) {
    for (k in 1:length(D_list)) {
      Sigma_k_hat = Sigma_hat[[k]]
      eigen_Sigma_k_hat = eigen(Sigma_k_hat)
      Sigma_hat[[k]] = eigen_Sigma_k_hat$vectors %*% diag(c(pmax(eigen_Sigma_k_hat$values, 0)), ncol(eigen_Sigma_k_hat$vectors), ncol(eigen_Sigma_k_hat$vectors)) %*% t(eigen_Sigma_k_hat$vectors)
      if (any(eigen_Sigma_k_hat$values < 0)) {
        attr(Sigma_hat[[k]], "truncated") = TRUE
      } else {
        attr(Sigma_hat[[k]], "truncated") = FALSE
      }
    }
  }

  if (highdim) {

    Sigma_hat = lapply(Sigma_hat, function(Sigma_r) s$v %*% Sigma_r %*% t(s$v))

  }

  return(list(Sigma_hat = Sigma_hat))

}
