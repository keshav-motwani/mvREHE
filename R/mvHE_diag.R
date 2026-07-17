#' mvHE_diag
#'
#' Like mvHE but uses only same-observation pairs (diagonal of D_k) in the
#' moment equations: G_k = Y' diag(D_k) Y, Q[k,z] = sum(diag(D_k)*diag(D_z)).
#' Equivalent to equation 13 of the mvREHE paper with i=l only.
#'
#' @param Y
#' @param D_list
#' @param W_row_pairs
#' @param w_columns
#' @param truncate
#' @return
#' @export
mvHE_diag = function(Y, D_list, W_row_pairs = NULL, w_columns = NULL, truncate = TRUE) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  q = ncol(Y); n = nrow(Y); K = length(D_list)
  sparse = all(sapply(D_list, function(x) is(x, "dsCMatrix")))

  if (!is.null(w_columns)) Y = t(t(Y) * sqrt(w_columns))

  highdim = q > n
  if (highdim) { s = svd(Y); Y = Y %*% s$v; q = ncol(Y) }

  # Diagonal of each D_k (replaces full D_k %*% Y)
  d_list = lapply(D_list, function(D) if (sparse) Matrix::diag(D) else diag(D))

  # Diagonal of W_row_pairs (replaces full W * D product)
  w_diag = if (!is.null(W_row_pairs)) {
    if (sparse) Matrix::diag(W_row_pairs) else diag(W_row_pairs)
  } else rep(1.0, n)

  # G_k = Y' diag(w * d_k) Y
  G_list = lapply(d_list, function(d) crossprod(Y, (w_diag * d) * Y))

  # Q[k,z] = sum(w * d_k * d_z)
  Q = matrix(NA, K, K)
  for (i in 1:K) for (j in 1:i)
    Q[i, j] = Q[j, i] = sum(w_diag * d_list[[i]] * d_list[[j]])

  G_mat    = do.call(rbind, lapply(G_list, as.vector))
  Sigma_mat = solve(Q, G_mat)

  Sigma_hat = lapply(seq_len(K), function(k) {
    S = matrix(Sigma_mat[k, ], q, q); (S + t(S)) / 2
  })

  if (truncate) {
    for (k in 1:K) {
      eig = eigen(Sigma_hat[[k]], symmetric = TRUE)
      Sigma_hat[[k]] = eig$vectors %*% (t(eig$vectors) * pmax(eig$values, 0))
      attr(Sigma_hat[[k]], "truncated") = any(eig$values < 0)
    }
  }

  if (highdim) {
    V = if (!is.null(w_columns)) t(t(s$v) / sqrt(w_columns)) else s$v
    Sigma_hat = lapply(Sigma_hat, function(S) V %*% S %*% t(V))
  } else if (!is.null(w_columns)) {
    W_inv = 1 / sqrt(w_columns)
    Sigma_hat = lapply(Sigma_hat, function(S) t(t(S * W_inv) * W_inv))
  }

  return(list(Sigma_hat = Sigma_hat))
}

#' mvREHE_diag
#'
#' Like mvREHE but uses only same-observation pairs (diagonal of D_k) in the
#' moment equations: G_k = Y' diag(D_k) Y, Q[k,z] = sum(diag(D_k)*diag(D_z)).
#'
#' @param Y
#' @param D_list
#' @param W_row_pairs
#' @param w_columns
#' @param truncate
#' @param tolerance
#' @param max_iter
#' @param Sigma_init_list
#' @param track_loss
#' @return
#' @export
mvREHE_diag = function(Y, D_list, W_row_pairs = NULL, w_columns = NULL,
                        tolerance = 1e-6, max_iter = 1000, return_full = TRUE,
                        Sigma_init_list = NULL, track_loss = FALSE, truncate = TRUE) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  sparse = all(sapply(D_list, function(x) is(x, "dsCMatrix")))

  n = nrow(Y); q = ncol(Y); K = length(D_list)
  difference   = numeric(max_iter)
  loss_history = if (track_loss) numeric(max_iter) else NULL

  if (!is.null(w_columns)) Y = t(t(Y) * sqrt(w_columns))

  highdim = q > n
  if (highdim) { s = svd(Y); Y = Y %*% s$v; q = ncol(Y) }

  if (is.null(Sigma_init_list)) {
    Sigma_list = lapply(1:K, function(i) matrix(0, q, q))
  } else if (is.character(Sigma_init_list) && Sigma_init_list == "mvHE") {
    Sigma_list = mvHE_diag(Y, D_list)$Sigma_hat
  } else {
    Sigma_list = Sigma_init_list
  }

  d_list = lapply(D_list, function(D) if (sparse) Matrix::diag(D) else diag(D))
  w_diag = if (!is.null(W_row_pairs)) {
    if (sparse) Matrix::diag(W_row_pairs) else diag(W_row_pairs)
  } else rep(1.0, n)

  G_list = lapply(d_list, function(d) crossprod(Y, (w_diag * d) * Y))

  Q = matrix(NA, K, K)
  for (i in 1:K) for (j in 1:i)
    Q[i, j] = Q[j, i] = sum(w_diag * d_list[[i]] * d_list[[j]])

  for (iter in 1:max_iter) {

    Sigma_list_old = Sigma_list

    for (z in seq(K, 1)) {
      mat = G_list[[z]]
      for (k in setdiff(1:K, z)) mat = mat - Sigma_list[[k]] * Q[k, z]
      if (truncate) {
        eig = eigen(mat, symmetric = TRUE)
        Sigma_list[[z]] = eig$vectors %*% (t(eig$vectors) * pmax(c(eig$values) / Q[z, z], 0))
      } else {
        Sigma_list[[z]] = mat / Q[z, z]
      }
    }

    if (track_loss) {
      lin  = sum(sapply(1:K, function(k) sum(Sigma_list[[k]] * G_list[[k]])))
      quad = sum(sapply(1:K, function(k) sum(sapply(1:K, function(z) Q[k,z] * sum(Sigma_list[[k]] * Sigma_list[[z]])))))
      loss_history[iter] = quad - 2 * lin
    }

    if (!is.null(tolerance)) {
      difference[iter] = mean(mapply(Sigma_list_old, Sigma_list,
        FUN = function(x, y) norm(x - y, "F") / (norm(x, "F") + 1e-10)), na.rm = TRUE)
      if (iter > 1 && (is.na(difference[iter]) || difference[iter] < tolerance)) break
    }
  }

  result = list(difference = difference[difference != 0])
  if (track_loss) result$loss = loss_history[loss_history != 0]

  if (highdim) {
    result$Sigma_hat = Sigma_list
    result$V = if (!is.null(w_columns)) t(t(s$v) / sqrt(w_columns)) else s$v
    if (return_full) {
      result$Sigma_hat = lapply(Sigma_list, function(S) result$V %*% S %*% t(result$V))
      result$V = NULL
    }
  } else {
    if (!is.null(w_columns)) {
      W_inv = 1 / sqrt(w_columns)
      Sigma_list = lapply(Sigma_list, function(S) t(t(S * W_inv) * W_inv))
    }
    result$Sigma_hat = Sigma_list
  }

  return(result)
}
