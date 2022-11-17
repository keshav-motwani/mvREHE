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
  lambda_seq = c(log_seq(lambda_max, lambda_min, n_lambda), 0) # / (q * (q - 1) / 2)

  L_init_list_orig = L_init_list

  loss = numeric(length(lambda_seq))

  for (l in 1:length(lambda_seq)) {

    print(l)

    loss_l = 0

    for (k in 1:K) {

      print(k)

      D_list_mk = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]])
      Sigma_hat = mvREHE(Y[-folds[[k]], ], D_list_mk, lambda_seq[l], tolerance, max_iter, L_init_list, algorithm)$Sigma_hat

      Y_tilde_list_k = lapply(1:q, function(j) lapply(1:j, function(m) c(tcrossprod(Y[folds[[k]], j], Y[folds[[k]], m]))))
      X_tilde_k = do.call(cbind, lapply(D_list, function(D) c(D[folds[[k]], folds[[k]]])))

      loss_l = loss_l + length(folds[[k]])^2 * mvREHE:::loss2(Y_tilde_list_k, X_tilde_k, Sigma_hat)

      # L_init_list = lapply(Sigma_hat, function(x) t(chol(x, pivot = TRUE)))

    }

    loss[l] = loss_l

  }

  lambda = lambda_seq[which.min(loss)]

  result = mvREHE(Y, D_list, lambda, tolerance, max_iter, L_init_list_orig, algorithm)
  result$cv_loss = loss
  result$lambda = lambda

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
mvREHE = function(Y, D_list, lambda = 0, tolerance = 1e-9, max_iter = 10000, L_init_list = NULL, algorithm = "L-BFGS-B") {

  q = ncol(Y)
  objective = numeric(max_iter)

  if (is.null(L_init_list)) {
    L_list = lapply(1:length(D_list), function(i) t(chol(clusterGeneration::rcorrmatrix(q))))
  } else if (is.character(L_init_list) && L_init_list == "mvHE") {
    mvHE_estimate = mvHE(Y, D_list)$Sigma_hat
    if (all(sapply(mvHE_estimate, function(x) attr(x, "truncated")) == 0)) return(list(Sigma_hat = mvHE_estimate))
    L_list = lapply(1:length(mvHE_estimate), function(i) t(chol(mvHE_estimate[[i]], pivot = TRUE)))
  } else {
    L_list = L_init_list
  }

  Y_tilde_list = lapply(1:q, function(j) lapply(1:j, function(m) c(tcrossprod(Y[, j], Y[, m]))))
  X_tilde = do.call(cbind, lapply(D_list, function(D) c(D)))
  gradient_list = replicate(length(D_list), matrix(0, q, q), simplify = FALSE)

  if (algorithm == "GD") {
    objective = fit_GD(Y_tilde_list, X_tilde, lambda, max_iter, tolerance, L_list, gradient_list)
  } else if (algorithm %in% c("L-BFGS-B")) {
    fit = fit_optim(Y_tilde_list, X_tilde, lambda, max_iter, L_list, algorithm)
    L_list = fit$L_list
    objective = fit$objective
  }

  return(list(Sigma_hat = lapply(L_list, function(L) tcrossprod(L)),
              objective = objective[objective != 0]))

}

# loss = function(Y_tilde_list, X_tilde, L_list) {
#
#   value = 0
#   for (j in 1:length(Y_tilde_list)) {
#     for (m in 1:length(Y_tilde_list[[j]])) {
#       sigma = sapply(L_list, function(L) tcrossprod(L[j, , drop = FALSE], L[m, , drop = FALSE]))
#       value = value + sum((Y_tilde_list[[j]][[m]] - (X_tilde %*% sigma)) ^ 2)
#     }
#   }
#
#   return(value)
#
# }

# gradient = function(Y_tilde_list, X_tilde, L_old_list, a, b, z) {
#
#   q = length(Y_tilde_list)
#
#   sigma = sapply(L_old_list, function(L) tcrossprod(L[a, , drop = FALSE]))
#   value = crossprod(-4 * L_old_list[[z]][a, b] * X_tilde[, z, drop = FALSE], Y_tilde_list[[a]][[a]] - (X_tilde %*% sigma))
#
#   if ((a - 1) >= b) {
#     for (m in b:(a - 1)) {
#       sigma = sapply(L_old_list, function(L) tcrossprod(L[a, , drop = FALSE], L[m, , drop = FALSE]))
#       value = value + crossprod(-2 * L_old_list[[z]][m, b] * X_tilde[, z, drop = FALSE], Y_tilde_list[[a]][[m]] - (X_tilde %*% sigma))
#     }
#   }
#
#   if (q >= (a + 1)) {
#     for (j in (a + 1):q) {
#       sigma = sapply(L_old_list, function(L) tcrossprod(L[j, , drop = FALSE], L[a, , drop = FALSE]))
#       value = value + crossprod(-2 * L_old_list[[z]][j, b] * X_tilde[, z, drop = FALSE], Y_tilde_list[[j]][[a]] - (X_tilde %*% sigma))
#     }
#   }
#
#   return(value[1, 1])
#
# }
#
# gradient_full = function(Y_tilde_list, X_tilde, L_old_list, gradient_list) {
#
#   q = length(Y_tilde_list)
#   K = ncol(X_tilde)
#
#   for (a in 1:q) {
#     sigma = sapply(L_old_list, function(L) tcrossprod(L[a, , drop = FALSE]))
#     U = Y_tilde_list[[a]][[a]] - X_tilde %*% sigma
#     for (z in 1:K) {
#       ZtU = crossprod(X_tilde[, z], U)
#       for (b in 1:a) {
#         gradient_list[[z]][a, b] = -4 * L_old_list[[z]][a, b] * ZtU
#       }
#     }
#   }
#
#   for (a in 1:q) {
#     if ((a - 1) >= 1) {
#       for (m in 1:(a - 1)) {
#         sigma = sapply(L_old_list, function(L) tcrossprod(L[a, , drop = FALSE], L[m, , drop = FALSE]))
#         V = Y_tilde_list[[a]][[m]] - X_tilde %*% sigma
#         for (z in 1:K) {
#           ZtV = crossprod(X_tilde[, z], V)
#           for (b in 1:a) {
#             if ((m >= b) & (m <= (a - 1))) {
#               gradient_list[[z]][a, b] = gradient_list[[z]][a, b] - 2 * L_old_list[[z]][m, b] * ZtV
#             }
#           }
#         }
#       }
#     }
#   }
#
#   for (a in 1:q) {
#     if (q >= a + 1) {
#       for (j in (a + 1):q) {
#         sigma = sapply(L_old_list, function(L) tcrossprod(L[j, , drop = FALSE], L[a, , drop = FALSE]))
#         W = Y_tilde_list[[j]][[a]] - X_tilde %*% sigma
#         for (z in 1:K) {
#           ZtW = crossprod(X_tilde[, z], W)
#           for (b in 1:a) {
#             gradient_list[[z]][a, b] = gradient_list[[z]][a, b] - 2 * L_old_list[[z]][j, b] * ZtW
#           }
#         }
#       }
#     }
#   }
#
#   return(gradient_list)
#
# }

# mvREHE = function(Y, D_list, tolerance = 1e-9, max_iter = 10000, init = FALSE) {
#
#   q = ncol(Y)
#   objective = numeric(max_iter)
#
#   if (init) {
#     mvHE_estimate = mvHE(Y, D_list[[1]], D_list[[2]])$Sigma_hat
#     if (all(sapply(mvHE_estimate, function(x) attr(x, "truncated")) == 0)) return(list(Sigma_hat = mvHE_estimate))
#     L_list = lapply(1:length(mvHE_estimate), function(i) t(chol(mvHE_estimate[[i]], pivot = TRUE)))
#   } else {
#     L_list = lapply(1:length(D_list), function(i) t(chol(clusterGeneration::rcorrmatrix(q))))
#   }
#
#   Y_tilde_list = lapply(1:q, function(j) lapply(1:j, function(m) c(tcrossprod(Y[, j], Y[, m]))))
#   X_tilde = do.call(cbind, lapply(D_list, function(D) c(D)))
#
#   grad = replicate(length(L_list), matrix(0, q, q), simplify = FALSE)
#
#   step_size = 1e-4
#
#   for (iter in 1:max_iter) {
#
#     L_old_list = L_list
#     loss_old = loss(Y_tilde_list, X_tilde, L_old_list)
#     gradient_full(Y_tilde_list, X_tilde, L_old_list, grad)
#
#     step_size = step_size * 2
#     line_search = TRUE
#
#     while(line_search) {
#
#       rhs = loss_old
#       for (z in 1:length(L_list)) {
#         L_list[[z]] = L_old_list[[z]] - step_size * grad[[z]]
#         rhs = rhs - 0.5 * step_size * sum(grad[[z]] ^ 2)
#       }
#       lhs = loss(Y_tilde_list, X_tilde, L_list)
#
#       if (is.na(lhs) | is.na(rhs) | lhs > rhs) {
#         step_size = 0.5 * step_size
#       } else {
#         line_search = FALSE
#       }
#
#     }
#
#     objective[iter] = loss(Y_tilde_list, X_tilde, L_list)
#
#     if (iter > 1 && objective[iter] > objective[iter - 1]) {
#       stop("Not decreasing")
#     }
#
#     if (iter > 1 && (objective[iter - 1] - objective[iter]) / objective[iter - 1] < tolerance) {
#       break
#     }
#
#   }
#
#   return(list(Sigma_hat = lapply(L_list, function(L) tcrossprod(L)),
#               objective = objective[objective != 0]))
#
# }

