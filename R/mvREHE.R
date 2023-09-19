#' cv_mvREHE_L2
#'
#' @param Y
#' @param D_list
#' @param K
#' @param n_rank
#' @param rank_max
#' @param rank_min
#' @param tolerance
#' @param max_iter
#' @param Sigma_init_list
#'
#' @return
#' @export
#'
#' @examples
cv_mvREHE_L2 = function(Y, D_list, K = 5, folds = NULL, grid = FALSE, n_lambda = 10, lambda_max = 1e-3, lambda_min = 1e-6, tolerance = NULL, max_iter = 100, Sigma_init_list = NULL) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y)/K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }

  q = ncol(Y)

  lambda_seq = seq(log10(lambda_max), log10(lambda_min), length.out = n_lambda)
  lambda_grid = do.call(expand.grid, replicate(length(D_list), lambda_seq, simplify = F))

  loss = numeric(nrow(lambda_grid))

  D_list_mk_list = vector(mode = "list", length = K)
  W_list_mk_list = vector(mode = "list", length = K)
  Q_mk_list = vector(mode = "list", length = K)
  D_list_k_list = vector(mode = "list", length = K)
  row_indices_mk_list = vector(mode = "list", length = K)
  col_indices_mk_list = vector(mode = "list", length = K)
  row_indices_k_list = vector(mode = "list", length = K)
  col_indices_k_list = vector(mode = "list", length = K)
  for (k in 1:K) {
    D_list_k_list[[k]] = lapply(D_list, function(D) D[folds[[k]], folds[[k]]])
    D_list_mk_list[[k]] = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]])
    Q_mk_list[[k]] = compute_Q(D_list_mk_list[[k]])
    indices = do.call(`+`, D_list_k_list[[k]]) > 0
    indices = indices & lower.tri(indices, diag = TRUE)
    row_indices_k_list[[k]] = which(indices, arr.ind = TRUE)[, 1] - 1
    col_indices_k_list[[k]] = which(indices, arr.ind = TRUE)[, 2] - 1
    indices = do.call(`+`, D_list_mk_list[[k]]) > 0
    indices = indices & lower.tri(indices, diag = TRUE)
    row_indices_mk_list[[k]] = which(indices, arr.ind = TRUE)[, 1] - 1
    col_indices_mk_list[[k]] = which(indices, arr.ind = TRUE)[, 2] - 1
    W_list_mk_list[[k]] = lapply(1:length(D_list), function(x) matrix(0, q, q))
    compute_W_list(Y[-folds[[k]], , drop = FALSE], D_list_mk_list[[k]], W_list_mk_list[[k]], row_indices_mk_list[[k]], col_indices_mk_list[[k]])
  }

  cv_loss = function(lambda) {
    loss_l = 0
    for (k in 1:K) {
      cat(k)
      Sigma_hat = mvREHE(Y[-folds[[k]], , drop = FALSE], D_list_mk_list[[k]], lambda = 10^lambda, tolerance, max_iter, Sigma_init_list, W_list = W_list_mk_list[[k]], Q = Q_mk_list[[k]], row_indices = row_indices_mk_list[[k]], col_indices = col_indices_mk_list[[k]])$Sigma_hat
      loss_l = loss_l + length(folds[[k]])^2 * loss(Y[folds[[k]], , drop = FALSE], D_list_k_list[[k]], Sigma_hat, row_indices_k_list[[k]], col_indices_k_list[[k]], rep(0, length(D_list)))
    }
    return(loss_l)
  }

  if (grid) {

    for (l in 1:nrow(lambda_grid)) {
      cat(l)
      loss[l] = cv_loss(as.numeric(lambda_grid[l, ]))
    }

    l = which.min(loss)
    lambda = 10^(lambda_grid[l, ])

  } else {

    lambda = 10^coef(optimx::optimx(rep(-5, length(D_list)), cv_loss, method = "L-BFGS-B", control = list(trace = 0, kkt =FALSE, starttests = FALSE)))

  }
  cat("\n")
  print(lambda)

  result = mvREHE(Y, D_list, lambda = as.numeric(lambda), tolerance, max_iter, Sigma_init_list)
  result$cv_loss = loss
  result$lambda = lambda
  result$lambda_grid = lambda_grid

  result

}

#' mvREHE
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
mvREHE = function(Y, D_list, lambda = NULL, tolerance = NULL, max_iter = 100, Sigma_init_list = NULL, W_list = NULL, Q = NULL, row_indices = NULL, col_indices = NULL) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  n = nrow(Y)
  q = ncol(Y)
  K = length(D_list)
  objective = numeric(max_iter)

  if (is.null(lambda)) lambda = rep(0, K)
  if (is.null(Sigma_init_list)) {
    Sigma_list = lapply(1:length(D_list), function(i) diag(1, q, q))
  } else if (is.character(Sigma_init_list) && Sigma_init_list == "mvHE") {
    Sigma_list = mvHE(Y, D_list)$Sigma_hat
  } else {
    Sigma_list = Sigma_init_list
  }
  if (is.null(row_indices) && (!is.null(tolerance) | is.null(W_list))) {
    indices = do.call(`+`, D_list) > 0
    indices = indices & lower.tri(indices, diag = TRUE)
    row_indices = which(indices, arr.ind = TRUE)[, 1] - 1
    col_indices = which(indices, arr.ind = TRUE)[, 2] - 1
  }
  if (is.null(W_list)) {
    W_list = lapply(1:K, function(x) matrix(0, q, q))
    compute_W_list(Y, D_list, W_list, row_indices, col_indices)
  }
  if (is.null(Q)) {
    Q = compute_Q(D_list)
  }

  for (iter in 1:max_iter) {

    for (z in 1:K) {
      mat = W_list[[z]]
      for (k in setdiff(1:K, z)) {
        mat = mat - Sigma_list[[k]] * Q[k, z]
      }
      eig = eigen(mat)
      Sigma_list[[z]] = eig$vectors %*% (t(eig$vectors) * pmax(eig$values / (Q[z, z] + (lambda[z] * n^2)), 0))
    }

    if (!is.null(tolerance)) {
      objective[iter] = loss(Y, D_list, Sigma_list, lambda, row_indices, col_indices)
      if (iter > 1 && abs(objective[iter - 1] - objective[iter]) / objective[iter - 1] < tolerance) {
        break
      }
    }

  }

  return(list(Sigma_hat = Sigma_list,
              objective = objective[objective != 0]))

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
