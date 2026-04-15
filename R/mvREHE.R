#' mvREHE
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
mvREHE = function(Y, D_list, W_row_pairs = NULL, w_columns = NULL, tolerance = 1e-6, max_iter = 1000, return_full = TRUE, Sigma_init_list = NULL) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  if (all(sapply(D_list, function(x) is(x, "dsCMatrix")))){
    sparse = TRUE
    if (!is.null(W_row_pairs)) {
      stopifnot(is(W_row_pairs, "dsCMatrix"))
    }
  } else {
    sparse = FALSE
    stopifnot(all(sapply(D_list, is.matrix)))
    if (!is.null(W_row_pairs)) {
      stopifnot(is.matrix(W_row_pairs))
    }
  }

  n = nrow(Y)
  q = ncol(Y)
  K = length(D_list)
  difference = numeric(max_iter)

  if (!is.null(w_columns)) {
    Y = t(t(Y) * sqrt(w_columns))
  }

  highdim = q > n
  if (highdim) {
    s = svd(Y)
    Y = Y %*% s$v
    q = ncol(Y)
  }

  if (is.null(Sigma_init_list)) {
    Sigma_list = lapply(1:length(D_list), function(i) matrix(0, q, q))
  } else if (is.character(Sigma_init_list) && Sigma_init_list == "mvHE") {
    Sigma_list = mvHE(Y, D_list)$Sigma_hat
  } else {
    Sigma_list = Sigma_init_list
  }

  if (is.null(W_row_pairs)) {
    G_list = lapply(D_list, function(D) crossprod(Y, as.matrix(D %*% Y)))
  } else {
    if (!sparse) {
      G_list = lapply(D_list, function(D) crossprod(Y, (W_row_pairs * D) %*% Y))
    } else {
      G_list = lapply(D_list, function(D) {
        WDY = compute_WDY(D, W_row_pairs, Y)
        as.matrix(crossprod(Y, WDY))
      })
    }
  }

  Q = compute_Q(D_list, W_row_pairs, sparse)

  for (iter in 1:max_iter) {

    Sigma_list_old = Sigma_list

    for (z in 1:K) {
      mat = G_list[[z]]
      for (k in setdiff(1:K, z)) {
        mat = mat - Sigma_list[[k]] * Q[k, z]
      }
      eig = eigen(mat, symmetric = TRUE)
      Sigma_list[[z]] = eig$vectors %*% (t(eig$vectors) * pmax(c(eig$values) / Q[z, z], 0))
    }

    if (!is.null(tolerance)) {
      difference[iter] = mean(mapply(Sigma_list_old, Sigma_list, FUN = function(x, y) norm(x - y, "F") / norm(x, "F")), na.rm = TRUE)
      if (iter > 1 && is.na(difference[iter])) {
        break
      }
      if (iter > 1 && difference[iter] < tolerance) {
        break
      }
    }

  }

  result = list(difference = difference[difference != 0])

  if (highdim) {

    result$Sigma_hat = Sigma_list
    result$V = s$v

    if (!is.null(w_columns)) {
      result$V = t(t(result$V) / sqrt(w_columns))
    }

    if (return_full) {
      result$Sigma_hat = lapply(Sigma_list, function(Sigma_r) result$V %*% Sigma_r %*% t(result$V))
      result$V = NULL
    }

  } else {

    if (!is.null(w_columns)) {
      Sigma_list = lapply(Sigma_list, function(Sigma) diag(1 / sqrt(w_columns)) %*% Sigma %*% diag(1 / sqrt(w_columns)))
    }

    result$Sigma_hat = Sigma_list

  }

  return(result)

}

compute_Q = function(D_list, W_row_pairs, sparse) {

  K = length(D_list)

  Q = matrix(NA, K, K)

  if (sparse) {

    for (i in 1:K) {
      for (j in 1:i) {
        Q[i, j] = Q[j, i] = frobenius_inner_product(D_list[[i]], D_list[[j]], W_row_pairs)
      }
    }

  } else {

    for (i in 1:K) {
      for (j in 1:i) {
        if (!is.null(W_row_pairs)) {
          Q[i, j] = Q[j, i] = sum(D_list[[i]] * D_list[[j]] * W_row_pairs)
        } else {
          Q[i, j] = Q[j, i] = sum(D_list[[i]] * D_list[[j]])
        }
      }
    }

  }

  return(Q)

}
