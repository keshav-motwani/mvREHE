#' mvHE
#'
#' @param Y
#' @param D_list
#' @param W_row_pairs
#' @param w_columns
#' @param truncate
#'
#' @return
#' @export
#'
#' @examples
mvHE = function(Y, D_list, W_row_pairs = NULL, w_columns = NULL, truncate = TRUE) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  q = ncol(Y)
  n = nrow(Y)
  K = length(D_list)

  sparse = all(sapply(D_list, function(x) is(x, "dsCMatrix")))

  if (!is.null(w_columns)) {
    Y = t(t(Y) * sqrt(w_columns))
  }

  highdim = q > n
  if (highdim) {
    s = svd(Y)
    Y = Y %*% s$v
    q = ncol(Y)
  }

  # Compute G matrices: G[[k]][j,m] = Y^T (W * D_k) Y  (W = identity if unweighted)
  if (is.null(W_row_pairs)) {
    G_list = lapply(D_list, function(D) crossprod(Y, as.matrix(D %*% Y)))
  } else if (sparse) {
    G_list = lapply(D_list, function(D) {
      as.matrix(crossprod(Y, compute_WDY(D, W_row_pairs, Y)))
    })
  } else {
    G_list = lapply(D_list, function(D) crossprod(Y, (W_row_pairs * D) %*% Y))
  }

  Q = compute_Q(D_list, W_row_pairs, sparse)

  # Solve Q * Sigma_mat = G_mat for all q² elements simultaneously (direct OLS)
  G_mat = do.call(rbind, lapply(G_list, as.vector))  # K × q²
  Sigma_mat = solve(Q, G_mat)                        # K × q²

  Sigma_hat = lapply(seq_len(K), function(k) {
    S = matrix(Sigma_mat[k, ], q, q)
    (S + t(S)) / 2
  })

  if (truncate) {
    for (k in 1:K) {
      eig = eigen(Sigma_hat[[k]], symmetric = TRUE)
      Sigma_hat[[k]] = eig$vectors %*% (t(eig$vectors) * pmax(eig$values, 0))
      attr(Sigma_hat[[k]], "truncated") = any(eig$values < 0)
    }
  }

  if (highdim) {
    if (!is.null(w_columns)) {
      V = t(t(s$v) / sqrt(w_columns))
    } else {
      V = s$v
    }
    Sigma_hat = lapply(Sigma_hat, function(S) V %*% S %*% t(V))
  } else {
    if (!is.null(w_columns)) {
      W_inv_sqrt = 1 / sqrt(w_columns)
      Sigma_hat = lapply(Sigma_hat, function(S) t(t(S * W_inv_sqrt) * W_inv_sqrt))
    }
  }

  return(list(Sigma_hat = Sigma_hat))

}
