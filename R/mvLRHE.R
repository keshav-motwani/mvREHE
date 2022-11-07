#' mvLRHE
#'
#' @param Y
#' @param D_list
#' @param tolerance
#' @param max_iter
#' @param Sigma_init_list
#'
#' @return
#' @export
#'
#' @examples
mvLRHE = function(Y, D_list, r, tolerance = 1e-9, max_iter = 10000, Sigma_init_list = NULL) {

  q = ncol(Y)
  K = length(D_list)
  objective = numeric(max_iter)

  if (is.null(Sigma_init_list)) {
    Sigma_list = lapply(1:length(D_list), function(i) clusterGeneration::rcorrmatrix(q))
  } else if (is.character(Sigma_init_list) && Sigma_init_list == "mvHE") {
    Sigma_list = mvHE(Y, D_list)$Sigma_hat
  } else {
    Sigma_list = Sigma_init_list
  }

  Y_tilde_list = lapply(1:q, function(j) lapply(1:j, function(m) c(tcrossprod(Y[, j], Y[, m]))))
  X_tilde = do.call(cbind, lapply(D_list, function(D) c(D)))

  W_list = compute_W_list(Y, D_list)
  Q = compute_Q(D_list)

  for (iter in 1:max_iter) {

    for (z in 1:K) {

      mat = W_list[[z]]
      for (k in setdiff(1:K, z)) {
        mat = mat - Sigma_list[[k]] * Q[k, z]
      }
      mat = mat / Q[z, z]

      eig = RSpectra::eigs_sym(mat, r)
      Sigma_list[[z]] = eig$vectors %*% (t(eig$vectors) * pmax(eig$values, 0))

    }

    objective[iter] = loss2(Y_tilde_list, X_tilde, Sigma_list)

    if (iter > 1 && objective[iter] > objective[iter - 1]) {
      stop("Not decreasing")
    }

    if (iter > 1 && abs(objective[iter - 1] - objective[iter]) / objective[iter - 1] < tolerance) {
      break
    }

  }

  return(list(Sigma_hat = Sigma_list,
              objective = objective[objective != 0]))

}

compute_W_list = function(Y, D_list) {

  q = ncol(Y)
  n = nrow(Y)
  K = length(D_list)

  W_list = vector(mode = "list", length = K)

  for (k in 1:K) {

    W = matrix(0, nrow = q, ncol = q)

    for (i in 1:n) {
      for (l in 1:n) {
        W = W + D_list[[k]][i, l] * crossprod(Y[i, , drop = FALSE], Y[l, , drop = FALSE])
      }
    }

    W_list[[k]] = W

  }

  return(W_list)

}

compute_Q = function(D_list) {

  K = length(D_list)

  Q = matrix(NA, K, K)

  for (i in 1:K) {
    for (j in 1:i) {
      Q[i, j] = Q[j, i] = sum(D_list[[i]] * D_list[[j]])
    }
  }

  return(Q)

}
