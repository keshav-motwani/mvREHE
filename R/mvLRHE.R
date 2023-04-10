#' oracle_mvLRHE
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
oracle_mvLRHE = function(Y, D_list, Sigma_list, n_rank = 5, rank_max = ncol(Y), rank_min = 1, tolerance = 1e-7, max_iter = 10000, Sigma_init_list = NULL) {

  q = ncol(Y)

  rank_seq = floor(seq(rank_min, rank_max, length.out = n_rank))
  rank_grid = do.call(expand.grid, replicate(length(D_list), rank_seq, simplify = F))

  loss = numeric(nrow(rank_grid))
  loss[] = NA

  for (l in 1:nrow(rank_grid)) {

    print(l)

    fit = mvLRHE(Y, D_list, as.numeric(rank_grid[l, ]), tolerance, max_iter, Sigma_init_list)
    Sigma_hat = fit$Sigma_hat

    loss[l] = sum(sapply(1:length(Sigma_hat), function(k) norm(Sigma_hat[[k]] - Sigma_list[[k]], "2")))

    if (loss[l] == min(loss, na.rm = T)) best_fit = fit

  }

  l = which.min(loss)
  rank = rank_grid[l, ]

  result = best_fit
  result$oracle_loss = loss
  result$rank = rank
  result$rank_grid = rank_grid

  result

}

#' cv_mvLRHE_rank
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
cv_mvLRHE_rank = function(Y, D_list, K, n_rank = 5, rank_max = ncol(Y), rank_min = 1, tolerance = 1e-7, max_iter = 10000, Sigma_init_list = NULL) {

  folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y)/K)))

  q = ncol(Y)

  rank_seq = floor(seq(rank_min, rank_max, length.out = n_rank))
  rank_grid = do.call(expand.grid, replicate(length(D_list), rank_seq, simplify = F))

  Sigma_init_list = Sigma_init_list

  loss = numeric(nrow(rank_grid))

  for (l in 1:nrow(rank_grid)) {

    print(l)

    loss_l = 0

    for (k in 1:K) {

      print(k)

      D_list_mk = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]])
      Sigma_hat = mvLRHE(Y[-folds[[k]], ], D_list_mk, r = as.numeric(rank_grid[l, ]), lambda = NULL, tolerance, max_iter, Sigma_init_list)$Sigma_hat

      X_tilde_k = do.call(cbind, lapply(D_list, function(D) c(D[folds[[k]], folds[[k]]])))

      loss_l = loss_l + length(folds[[k]])^2 * loss3(Y[folds[[k]], ], X_tilde_k, Sigma_hat, rep(0, length(D_list)))

    }

    loss[l] = loss_l

  }

  l = which.min(loss)
  rank = rank_grid[l, ]
  print(rank)

  result = mvLRHE(Y, D_list, r = rank, lambda = NULL, tolerance, max_iter, Sigma_init_list)
  result$cv_loss = loss
  result$rank = rank
  result$rank_grid = rank_grid

  result

}

#' cv_mvLRHE_rank
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
cv_mvLRHE_L2 = function(Y, D_list, K, n_lambda = 5, lambda_max = 1e-4, lambda_min = 1e-12, tolerance = 1e-7, max_iter = 10000, Sigma_init_list = NULL) {

  folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y)/K)))

  q = ncol(Y)

  lambda_seq = c(log_seq(lambda_max, lambda_min, n_lambda), 0)
  lambda_grid = do.call(expand.grid, replicate(length(D_list), lambda_seq, simplify = F))

  Sigma_init_list = Sigma_init_list

  loss = numeric(nrow(lambda_grid))

  for (l in 1:nrow(lambda_grid)) {

    print(l)

    loss_l = 0

    for (k in 1:K) {

      print(k)

      D_list_mk = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]])
      Sigma_hat = mvLRHE(Y[-folds[[k]], ], D_list_mk, r = NULL, lambda = as.numeric(lambda_grid[l, ]), tolerance, max_iter, Sigma_init_list)$Sigma_hat

      X_tilde_k = do.call(cbind, lapply(D_list, function(D) c(D[folds[[k]], folds[[k]]])))

      loss_l = loss_l + length(folds[[k]])^2 * loss3(Y[folds[[k]], ], X_tilde_k, Sigma_hat, rep(0, length(D_list)))

    }

    loss[l] = loss_l

  }

  l = which.min(loss)
  lambda = lambda_grid[l, ]
  print(lambda)

  result = mvLRHE(Y, D_list, r = NULL, lambda = as.numeric(lambda), tolerance, max_iter, Sigma_init_list)
  result$cv_loss = loss
  result$lambda = lambda
  result$lambda_grid = lambda_grid

  result

}

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
mvLRHE = function(Y, D_list, r = NULL, lambda = NULL, tolerance = 1e-7, max_iter = 10000, Sigma_init_list = NULL) {

  n = nrow(Y)
  q = ncol(Y)
  K = length(D_list)
  objective = numeric(max_iter)

  if (is.null(r)) r = rep(q, K)
  if (is.null(lambda)) lambda = rep(0, K)

  if (is.null(Sigma_init_list)) {
    Sigma_list = lapply(1:length(D_list), function(i) clusterGeneration::rcorrmatrix(q))
  } else if (is.character(Sigma_init_list) && Sigma_init_list == "mvHE") {
    Sigma_list = mvHE(Y, D_list)$Sigma_hat
  } else {
    Sigma_list = Sigma_init_list
  }

  X_tilde = do.call(cbind, lapply(D_list, function(D) c(D)))

  W_list = lapply(1:K, function(x) matrix(0, q, q))
  compute_W_list(Y, D_list, W_list)
  Q = compute_Q(D_list)

  for (iter in 1:max_iter) {

    for (z in 1:K) {

      mat = W_list[[z]]
      for (k in setdiff(1:K, z)) {
        mat = mat - Sigma_list[[k]] * Q[k, z]
      }

      eig = RSpectra::eigs_sym(mat, r[z])
      Sigma_list[[z]] = eig$vectors %*% (t(eig$vectors) * pmax(eig$values / (Q[z, z] + (lambda[z] * n^2 * q^2)), 0))

    }

    if (iter %% 10 == 0) {

      it = (iter %/% 10)

      objective[it] = loss3(Y, X_tilde, Sigma_list, lambda)

      if (it > 1 && objective[it] > objective[it - 1]) {
        print(c(objective[1:it]))
      }

      if (it > 1 && abs(objective[it - 1] - objective[it]) / objective[it - 1] < tolerance) {
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
