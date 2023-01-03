#' mvREHE
#'
#' @param Y
#' @param D_list
#' @param Sigma_list
#' @param n_lambda
#' @param lambda_max
#' @param lambda_min
#' @param tolerance
#' @param max_iter
#' @param L_init_list
#' @param algorithm
#'
#' @return
#' @export
#'
#' @examples
oracle_mvREHE = function(Y, D_list, Sigma_list, n_lambda = 10, lambda_max = 1, lambda_min = 1e-8, tolerance = 1e-9, max_iter = 10000, L_init_list = NULL, algorithm = "L-BFGS-B") {

  q = ncol(Y)

  lambda_seq = c(log_seq(lambda_max, lambda_min, n_lambda), 0)
  lambda_grid = do.call(expand.grid, replicate(length(D_list), lambda_seq, simplify = F))

  L_init_list_orig = L_init_list

  loss = numeric(nrow(lambda_grid))
  loss[] = NA

  for (l in 1:nrow(lambda_grid)) {

    print(l)

    fit = mvREHE(Y, D_list, lambda_grid[l, ], tolerance, max_iter, L_init_list, algorithm)
    Sigma_hat = fit$Sigma_hat
    L_hat = fit$L_hat

    loss[l] = sum(sapply(1:length(Sigma_hat), function(k) norm(Sigma_hat[[k]] - Sigma_list[[k]], "2")))

    if (loss[l] = min(loss, na.rm = T)) best_fit = fit

    L_init_list = L_hat

  }

  lambda = lambda_grid[which.min(loss), ]

  result = best_fit
  result$oracle_loss = loss
  result$lambda = lambda

  result

}

#' mvREHE
#'
#' @param Y
#' @param D_list
#' @param K
#' @param n_lambda
#' @param lambda_max
#' @param lambda_min
#' @param tolerance
#' @param max_iter
#' @param L_init_list
#' @param algorithm
#'
#' @return
#' @export
#'
#' @examples
cv_mvREHE = function(Y, D_list, K, n_lambda = 10, lambda_max = 1, lambda_min = 1e-8, tolerance = 1e-9, max_iter = 10000, L_init_list = NULL, algorithm = "L-BFGS-B") {

  folds = split(sample(1:nrow(Y), nrow(Y)), 1:K)

  q = ncol(Y)

  lambda_seq = c(log_seq(lambda_max, lambda_min, n_lambda), 0)
  lambda_grid = do.call(expand.grid, replicate(length(D_list), lambda_seq, simplify = F))

  L_init_list_orig = L_init_list

  loss = numeric(nrow(lambda_grid))

  for (l in 1:nrow(lambda_grid)) {

    print(l)

    loss_l = 0

    for (k in 1:K) {

      print(k)

      D_list_mk = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]])
      fit = mvREHE(Y[-folds[[k]], ], D_list_mk, lambda_grid[l, ], tolerance, max_iter, L_init_list, algorithm)
      Sigma_hat = fit$Sigma_hat
      L_hat = fit$L_hat

      Y_tilde_list_k = lapply(1:q, function(j) lapply(1:j, function(m) c(tcrossprod(Y[folds[[k]], j], Y[folds[[k]], m]))))
      X_tilde_k = do.call(cbind, lapply(D_list, function(D) c(D[folds[[k]], folds[[k]]])))

      loss_l = loss_l + length(folds[[k]])^2 * loss2(Y_tilde_list_k, X_tilde_k, Sigma_hat)

      L_init_list = L_hat

    }

    loss[l] = loss_l

  }

  lambda = lambda_grid[which.min(loss), ]

  result = mvREHE(Y, D_list, lambda, tolerance, max_iter, L_init_list_orig, algorithm)
  result$cv_loss = loss
  result$lambda = lambda
  result$lambda_grid = lambda_grid

  result

}

log_seq = function(from, to, length) {

  sequence = exp(seq(log(from), log(to), length.out = length))

  sequence[1] = from
  if (length > 1) sequence[length] = to

  sequence

}

#' mvREHE
#'
#' @param Y
#' @param D_list
#' @param tolerance
#' @param max_iter
#' @param L_init_list
#' @param algorithm
#'
#' @return
#' @export
#'
#' @examples
mvREHE = function(Y, D_list, lambda = NULL, tolerance = 1e-9, max_iter = 10000, L_init_list = NULL, algorithm = "L-BFGS-B") {

  q = ncol(Y)
  objective = numeric(max_iter)

  if (is.null(L_init_list)) {
    L_list = lapply(1:length(D_list), function(i) t(chol(clusterGeneration::rcorrmatrix(q))))
  } else {
    L_list = L_init_list
  }

  if (is.null(lambda)) {
    lambda = rep(0, length(D_list))
  }

  Y_tilde_list = lapply(1:q, function(j) lapply(1:j, function(m) c(tcrossprod(Y[, j], Y[, m]))))
  X_tilde = do.call(cbind, lapply(D_list, function(D) c(D)))
  gradient_list = replicate(length(D_list), matrix(0, q, q), simplify = FALSE)

  fit = fit_optim(Y_tilde_list, X_tilde, lambda, max_iter, tolerance, L_list, algorithm)
  L_list = fit$L_list
  objective = fit$objective

  return(list(Sigma_hat = lapply(L_list, function(L) tcrossprod(L)),
              L_hat = L_list,
              objective = objective[objective != 0]))

}
