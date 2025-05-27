matrix_regression_from_cov = function(Sigma_hat, V, rank, lambda, max_iter = 1000, tolerance = 1e-6, n_init = 5) {

  if (!is.null(V)) {
    p = sqrt(nrow(V) - 1)
  } else {
    p = sqrt(nrow(Sigma_hat) - 1)
  }

  best_objective = Inf
  best_solution = NULL

  for (init in 1:n_init) {

    beta1 = matrix(rnorm(p * rank, sd = 1), p, rank)
    beta2 = beta1

    objective = numeric(2 * max_iter)
    difference = numeric(max_iter)

    for (iter in 1:max_iter) {

      beta1_old = beta1
      beta2_old = beta2

      if (!is.null(V)) {
        # A = compute_A_for_Sigma_y_Xbeta2(beta2, V)
        A = rbind(c(1, rep(0, p^2)), cbind(rep(0, p * rank), kronecker(t(beta2), diag(1, p, p)))) %*% V
      } else {
        A = rbind(c(1, rep(0, p^2)), cbind(rep(0, p * rank), kronecker(t(beta2), diag(1, p, p))))
      }
      Sigma = A %*% Sigma_hat %*% t(A)
      beta1[] = solve(Sigma[-1, -1] + diag(lambda, ncol(Sigma) - 1, ncol(Sigma) - 1), Sigma[1, -1])

      if (!is.null(V)) {
        # A = compute_A_for_objective(beta1, beta2, V)
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2)))) %*% V
      } else {
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2))))
      }
      Sigma = A %*% Sigma_hat %*% t(A)
      objective[2 * iter - 1] = Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2] + lambda * (sum(beta1^2) + sum(beta2^2))

      if (!is.null(V)) {
        # A = compute_A_for_Sigma_y_Xtbeta1(beta1, V)
        A = rbind(c(1, rep(0, p^2)), cbind(rep(0, p * rank), kronecker(diag(1, p, p), t(beta1)))) %*% V
      } else {
        A = rbind(c(1, rep(0, p^2)), cbind(rep(0, p * rank), kronecker(diag(1, p, p), t(beta1))))
      }
      Sigma = A %*% Sigma_hat %*% t(A)
      beta2[] = matrix(solve(Sigma[-1, -1] + diag(lambda, ncol(Sigma) - 1, ncol(Sigma) - 1), Sigma[1, -1]), ncol = rank, byrow = TRUE)

      if (!is.null(V)) {
        # A = compute_A_for_objective(beta1, beta2, V)
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2)))) %*% V
      } else {
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2))))
      }
      Sigma = A %*% Sigma_hat %*% t(A)
      objective[2 * iter] = Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2] + lambda * (sum(beta1^2) + sum(beta2^2))
      difference[iter] = norm(beta1 %*% t(beta2) - beta1_old %*% t(beta2_old), "F") / (norm(beta1 %*% t(beta2), "F") + 1e-12)

      if (iter > 5 && difference[iter] < tolerance) {
        break
      }
    }

    final_objective = objective[max(which(objective != 0))]
    if (final_objective < best_objective) {
      best_objective = final_objective
      best_solution = list(
        beta1 = beta1,
        beta2 = beta2,
        difference = difference[difference != 0],
        objective = objective[objective != 0],
        init_number = init,
        final_objective = final_objective
      )
    }
  }

  return(best_solution)
}

latent_matrix_regression = function(fit, component, outcome, covariates, rank, lambda, ...) {

  Sigma_hat = fit$Sigma_hat[[component]]

  covariates = vech_to_vec_indices(covariates)

  if (!is.null(fit$V)) {
    V = fit$V[c(outcome, covariates), ]
    Sigma_hat = V %*% Sigma_hat %*% t(V)
    V = NULL
  } else {
    V = NULL
    Sigma_hat = Sigma_hat[c(outcome, covariates), c(outcome, covariates)]
  }

  Sigma_hat = cov2cor(Sigma_hat)
  Sigma_hat[is.na(Sigma_hat)] = 0

  matrix_regression_from_cov(Sigma_hat, V, rank, lambda, ...)

}

latent_ridge_regression = function(fit, component, outcome, covariates, lambda) {

  if (is.null(fit$V)) {
    Sigma_hat = cov2cor(fit$Sigma_hat[[component]])
  } else {
    Sigma_hat = cov2cor(fit$V %*% fit$Sigma_hat[[component]] %*% t(fit$V))
  }
  Sigma_hat[is.na(Sigma_hat)] = 0
  Sigma_hat = Sigma_hat[c(outcome, covariates), c(outcome, covariates)]

  return(solve(Sigma_hat[-1, -1] + diag(lambda, length(covariates)), Sigma_hat[-1, 1]))

}

raw_matrix_regression = function(Y, outcome, covariates, rank, lambda, ...) {

  Sigma_hat = cor(Y)
  Sigma_hat[is.na(Sigma_hat)] = 0

  covariates = vech_to_vec_indices(covariates)

  V = NULL
  Sigma_hat = Sigma_hat[c(outcome, covariates), c(outcome, covariates)]

  matrix_regression_from_cov(Sigma_hat, V, rank, lambda, ...)

}

raw_ridge_regression = function(Y, outcome, covariates, lambda) {

  Sigma_hat = cor(Y)
  Sigma_hat[is.na(Sigma_hat)] = 0

  Sigma_hat = Sigma_hat[c(outcome, covariates), c(outcome, covariates)]

  return(solve(Sigma_hat[-1, -1] + diag(lambda, length(covariates)), Sigma_hat[-1, 1]))

}

vech_to_vec_indices <- function(vech_indices) {

  q <- length(vech_indices)
  p <- (sqrt(8*q + 1) - 1)/2

  if(p != round(p)) {
    stop("Input vector length is not valid for a symmetric matrix vech representation")
  }

  mat = matrix(0, p, p)

  mat[lower.tri(mat, diag = TRUE)] = vech_indices

  c(mat + t(mat) - diag(diag(mat), p, p))

}

cv_latent_matrix_regression <- function(Y, D_list, outcomes, covariates, rank_seq, lambda_seq, estimator,
                                        K = 2, folds = NULL, cores = 8, ...) {

  require(parallel)

  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y) / K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }

  covariates = vech_to_vec_indices(covariates)

  cv_r2 = array(0, dim = c(length(D_list), length(outcomes), length(rank_seq), length(lambda_seq)))

  for (k in 1:K) {

    fit_train = estimator(Y[-folds[[k]], ], D_list = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]]))
    fit_test = estimator(Y[folds[[k]], ], D_list = lapply(D_list, function(D) D[folds[[k]], folds[[k]]]))

    for (c in 1:length(D_list)) {

      # Create full grid of outcome × rank × lambda
      grid <- expand.grid(
        o_idx = seq_along(outcomes),
        r_idx = seq_along(rank_seq),
        l_idx = seq_along(lambda_seq)
      )

      results_grid <- mclapply(seq_len(nrow(grid)), function(idx) {
        row <- grid[idx, ]
        o_idx <- row$o_idx
        r_idx <- row$r_idx
        l_idx <- row$l_idx

        o <- outcomes[o_idx]
        r <- rank_seq[r_idx]
        l <- lambda_seq[l_idx]

        # Training
        Sigma_hat_train <- fit_train$Sigma_hat[[c]]
        if (!is.null(fit_train$V)) {
          V_train <- fit_train$V[c(o, covariates), ]
          Sigma_hat_train = V_train %*% Sigma_hat_train %*% t(V_train)
          V_train = NULL
        } else {
          V_train = NULL
          Sigma_hat_train = Sigma_hat_train[c(o, covariates), c(o, covariates)]
        }
        Sigma_hat_train = cov2cor(Sigma_hat_train)
        Sigma_hat_train[is.na(Sigma_hat_train)] = 0
        fit <- matrix_regression_from_cov(Sigma_hat_train, V_train, r, l, ...)

        # Testing
        Sigma_hat_test <- fit_test$Sigma_hat[[c]]
        if (!is.null(fit_test$V)) {
          V_test <- fit_test$V[c(o, covariates), ]
          Sigma_hat_test = V_test %*% Sigma_hat_test %*% t(V_test)
          V_test = NULL
        } else {
          V_test = NULL
          Sigma_hat_test = Sigma_hat_test[c(o, covariates), c(o, covariates)]
        }
        Sigma_hat_test = cov2cor(Sigma_hat_test)
        Sigma_hat_test[is.na(Sigma_hat_test)] = 0

        A <- rbind(
          c(1, rep(0, length(covariates))),
          c(0, c(tcrossprod(fit$beta1, fit$beta2)))
        )

        if (!is.null(V_test)) {
          A = A %*% V_test
        }

        Sigma <- A %*% Sigma_hat_test %*% t(A)
        r2 <- 1 - (Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2]) / Sigma[1, 1]

        list(o_idx = o_idx, r_idx = r_idx, l_idx = l_idx, r2 = r2)
      }, mc.cores = cores)

      # Accumulate results
      for (res in results_grid) {
        cv_r2[c, res$o_idx, res$r_idx, res$l_idx] <-
          cv_r2[c, res$o_idx, res$r_idx, res$l_idx] + res$r2 / K
      }
    }
  }

  # Select best rank/lambda based on max CV R²
  lambda = matrix(NA, length(D_list), length(outcomes))
  rank = matrix(NA, length(D_list), length(outcomes))

  for (c in 1:length(D_list)) {
    for (o in 1:length(outcomes)) {
      idx = which(cv_r2[c, o, , , drop = FALSE] == max(cv_r2[c, o, , ]), arr.ind = TRUE)[1, ]
      rank[c, o] = rank_seq[idx[3]]
      lambda[c, o] = lambda_seq[idx[4]]
    }
  }

  return(list(
    lambda = lambda,
    rank = rank,
    cv_r2 = pmax(apply(cv_r2, c(1, 2), max), 0),
    cv_r2_full = cv_r2
  ))

}

cv_latent_ridge_regression = function(Y, D_list, fit, outcomes, covariates, n_lambda, estimator, K = 2, folds = NULL) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y)/K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }

  components = 1:length(D_list)

  lambda_grid = array(NA, dim = c(length(components), n_lambda))
  cv_r2 = array(0, dim = c(length(components), length(outcomes), n_lambda))

  for (k in 1:K) {

    fit_train = estimator(Y[-folds[[k]], ], D_list = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]]))
    fit_test = estimator(Y[folds[[k]], ], D_list = lapply(D_list, function(D) D[folds[[k]], folds[[k]]]))

    for (c in 1:length(components)) {

      component = components[c]

      cov_hat = cov2cor(fit$Sigma_hat[[component]])
      cov_hat[is.na(cov_hat)] = 0
      cov_hat_train = cov2cor(fit_train$Sigma_hat[[component]])
      cov_hat_train[is.na(cov_hat_train)] = 0
      cov_hat_test = cov2cor(fit_test$Sigma_hat[[component]])
      cov_hat_test[is.na(cov_hat_test)] = 0

      max_eigenvalue = max(eigen(cov_hat)$val)
      lambda_grid[component, ] = 10^seq(log10(max_eigenvalue / 1000), log10(max_eigenvalue), length.out = 100)

      for (o in 1:length(outcomes)) {

        outcome = outcomes[o]

        for (l in 1:n_lambda) {

          print(l)

          beta_hat = solve(cov_hat_train[covariates, covariates] + diag(lambda_grid[c, l], length(covariates), length(covariates)), cov_hat_train[covariates, outcome])

          r2 = 1 - (cov_hat_test[outcome, outcome] - 2 * cov_hat_test[outcome, covariates] %*% beta_hat + t(beta_hat) %*% cov_hat_test[covariates, covariates] %*% beta_hat) / cov_hat_test[outcome, outcome]
          cv_r2[component, o, l] = cv_r2[component, o, l] + r2 / K

        }

      }

    }

  }

  lambda = matrix(NA, length(components), length(outcomes))
  for (c in 1:length(components)) {
    for (o in 1:length(outcomes)) {
      lambda[c, o] = lambda_grid[c, ][which.max(cv_r2[c, o, ])]
    }
  }

  return(list(lambda = lambda, cv_r2 = pmax(apply(cv_r2, c(1, 2), max), 0), cv_r2_full = cv_r2))

}


