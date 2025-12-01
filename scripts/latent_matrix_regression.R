Rcpp::sourceCpp("scripts/matrix_regression_from_cov.cpp")

# library(Rcpp)
# 
# cppFunction('
# // Compute A1 * B for any matrix B where A1 is matrix needed to update beta1
# NumericMatrix compute_A1B(NumericMatrix beta2, NumericMatrix B) {
#   int p = beta2.nrow();
#   int rank = beta2.ncol();
#   int numCols = B.ncol();
#   NumericMatrix result(1 + p * rank, numCols);
# 
#   // First row is B[0,]
#   for(int j = 0; j < numCols; j++) {
#     result(0, j) = B(0, j);
#   }
# 
#   // Compute A = rbind(c(1, rep(0, p^2)),
#   //                   cbind(rep(0, p * rank),
#   //                         kronecker(t(beta2), diag(1, p, p)))) %*% B
#   for(int r = 0; r < rank; r++) {
#     for(int i = 0; i < p; i++) {
#       for(int col = 0; col < numCols; col++) {
#         double sum = 0.0;
#         for(int k = 0; k < p; k++) {
#           sum += beta2(k, r) * B(1 + i*p + k, col);
#         }
#         result(1 + r*p + i, col) = sum;
#       }
#     }
#   }
# 
#   return result;
# }')
# 
# cppFunction('
# // Compute A2 * B for any matrix B where A2 is matrix needed to update beta2
# NumericMatrix compute_A2B(NumericMatrix beta1, NumericMatrix B) {
#   int p = beta1.nrow();
#   int rank = beta1.ncol();
#   int numCols = B.ncol();
#   NumericMatrix result(1 + p * rank, numCols);
#   // First row is B[0,]
#   for(int j = 0; j < numCols; j++) {
#     result(0, j) = B(0, j);
#   }
#   // Compute A = rbind(c(1, rep(0, p^2)),
#   //                   cbind(rep(0, p * rank),
#   //                         kronecker(diag(1, p, p), t(beta1)))) %*% B
#   for(int r = 0; r < rank; r++) {
#     for(int i = 0; i < p; i++) {
#       for(int col = 0; col < numCols; col++) {
#         double sum = 0.0;
#         for(int k = 0; k < p; k++) {
#           sum += beta1(k, r) * B(1 + i*p + k, col);
#         }
#         result(1 + i*rank + r, col) = sum;
#       }
#     }
#   }
#   return result;
# }')
# 
# matrix_regression_from_cov = function(Sigma_hat,
#                                       V,
#                                       rank,
#                                       lambda,
#                                       max_iter = 1000,
#                                       tolerance = 1e-6,
#                                       n_init = 5) {
#   if (!is.null(V)) {
#     p = sqrt(nrow(V) - 1)
#   } else {
#     p = sqrt(nrow(Sigma_hat) - 1)
#   }
# 
#   best_objective = Inf
#   best_solution = NULL
# 
#   for (init in 1:n_init) {
#     beta1 = matrix(rnorm(p * rank, sd = 1), p, rank)
#     beta2 = beta1
# 
#     objective = numeric(max_iter)
#     difference = numeric(max_iter)
# 
#     for (iter in 1:max_iter) {
#       beta1_old = beta1
#       beta2_old = beta2
# 
#       # A1 = rbind(c(1, rep(0, p ^ 2)), cbind(rep(0, p * rank), kronecker(t(beta2), diag(1, p, p))))
#       if (!is.null(V)) {
#         A = compute_A1B(beta2, V)
#         Sigma = A %*% Sigma_hat %*% t(A)
#       } else {
#         temp = compute_A1B(beta2, Sigma_hat)
#         Sigma = t(compute_A1B(beta2, t(temp)))
#       }
#       beta1[] = solve(Sigma[-1,-1] + diag(lambda, ncol(Sigma) - 1, ncol(Sigma) - 1), Sigma[1,-1])
# 
#       # A2 = rbind(c(1, rep(0, p ^ 2)), cbind(rep(0, p * rank), kronecker(diag(1, p, p), t(beta1))))
#       if (!is.null(V)) {
#         A = compute_A2B(beta1, V)
#         Sigma = A %*% Sigma_hat %*% t(A)
#       } else {
#         temp = compute_A2B(beta1, Sigma_hat)
#         Sigma = t(compute_A2B(beta1, t(temp)))
#       }
#       beta2[] = matrix(solve(Sigma[-1,-1] + diag(
#         lambda, ncol(Sigma) - 1, ncol(Sigma) - 1
#       ), Sigma[1,-1]),
#       ncol = rank,
#       byrow = TRUE)
# 
#       A = rbind(c(1, rep(0, p ^ 2)), c(0, c(tcrossprod(
#         beta1, beta2
#       ))))
# 
#       if (!is.null(V)) {
#         A = A %*% V
#       }
# 
#       Sigma = A %*% Sigma_hat %*% t(A)
#       objective[iter] = Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2] + lambda * (sum(beta1 ^
#                                                                                       2) + sum(beta2 ^ 2))
#       difference[iter] = norm(beta1 %*% t(beta2) - beta1_old %*% t(beta2_old), "F") / (norm(beta1 %*% t(beta2), "F") + 1e-12)
# 
#       if (iter > 5 && difference[iter] < tolerance) {
#         break
#       }
#     }
# 
#     if (iter == max_iter) {
#       warning(paste("Did not converge"))
#     }
# 
#     final_objective = objective[iter]
#     if (final_objective < best_objective) {
#       print(c(final_objective, best_objective))
#       best_objective = final_objective
#       best_solution = list(
#         beta1 = beta1,
#         beta2 = beta2,
#         difference = difference[difference != 0],
#         objective = objective[objective != 0],
#         init_number = init,
#         final_objective = final_objective
#       )
#     }
#   }
# 
#   return(best_solution)
# }

latent_matrix_regression = function(fit,
                                    component,
                                    outcome,
                                    covariates,
                                    rank,
                                    lambda,
                                    ...) {
  
  Sigma_hat = fit$Sigma_hat[[component]]
  
  # covariates = vech_to_vec_indices(covariates)
  
  if (!is.null(fit$V)) {
    V = fit$V[c(outcome, covariates),]
    Sigma_hat = V %*% Sigma_hat %*% t(V)
    V = NULL
  } else {
    V = NULL
    Sigma_hat = Sigma_hat[c(outcome, covariates), c(outcome, covariates)]
  }
  
  Sigma_hat = cov2cor_NA0(Sigma_hat)

  matrix_regression_from_cov(Sigma_hat, V, rank, lambda, ...)
  
}

raw_matrix_regression = function(Y, outcome, covariates, rank, lambda, ...) {
  Sigma_hat = cov2cor_NA0(cov(Y))

  # covariates = vech_to_vec_indices(covariates)
  
  V = NULL
  Sigma_hat = Sigma_hat[c(outcome, covariates), c(outcome, covariates)]
  
  matrix_regression_from_cov(Sigma_hat, V, rank, lambda, ...)
  
}

vech_to_vec_indices = function(vech_indices) {
  q = length(vech_indices)
  p = (sqrt(8 * q + 1) - 1) / 2
  
  if (p != round(p)) {
    stop("Input vector length is not valid for a symmetric matrix vech representation")
  }
  
  mat = matrix(0, p, p)
  
  mat[lower.tri(mat, diag = TRUE)] = vech_indices
  
  c(mat + t(mat) - diag(diag(mat), p, p))
  
}

cv_latent_matrix_regression = function(Y,
                                       D_list,
                                       outcomes,
                                       covariates,
                                       rank_seq,
                                       lambda_seq,
                                       estimator,
                                       K = 2,
                                       folds = NULL,
                                       cores = 8,
                                       ...) {
  require(parallel)

  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y) / K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }

  # covariates = vech_to_vec_indices(covariates)

  cv_r2 = array(0, dim = c(
    length(D_list),
    length(outcomes),
    length(rank_seq),
    length(lambda_seq)
  ))

  for (k in 1:K) {
    
    print(K)
    
    fit_train = estimator(Y[-folds[[k]],], D_list = lapply(D_list, function(D)
      D[-folds[[k]],-folds[[k]]]))
    fit_test = estimator(Y[folds[[k]],], D_list = lapply(D_list, function(D)
      D[folds[[k]], folds[[k]]]))

    for (c in 1:length(D_list)) {
      
      print(c)
      
      # Create full grid of outcome × rank × lambda
      grid = expand.grid(
        o_idx = seq_along(outcomes),
        r_idx = seq_along(rank_seq),
        l_idx = seq_along(lambda_seq)
      )

      Sigma_hat_train = fit_train$Sigma_hat[[c]]
      if (!is.null(fit_train$V)) {
        sds = rowSums((fit_train$V %*% Sigma_hat_train) * fit_train$V)
        fit_train$V[sds < .Machine$double.eps, ] = 0
        sds[sds < .Machine$double.eps] = 1
        V_train = fit_train$V / sqrt(sds)
      } else {
        Sigma_hat_train = cov2cor_NA0(Sigma_hat_train)
        V_train = NULL
      }
      
      Sigma_hat_test = fit_test$Sigma_hat[[c]]
      if (!is.null(fit_test$V)) {
        sds = rowSums((fit_test$V %*% Sigma_hat_test) * fit_test$V)
        fit_test$V[sds < .Machine$double.eps, ] = 0
        sds[sds < .Machine$double.eps] = 1
        V_test = fit_test$V / sqrt(sds)
      } else {
        Sigma_hat_test = cov2cor_NA0(Sigma_hat_test)
        V_test = NULL
      }

      results_grid = mclapply(seq_len(nrow(grid)), function(idx) {
        row = grid[idx,]
        o_idx = row$o_idx
        r_idx = row$r_idx
        l_idx = row$l_idx

        o = outcomes[o_idx]
        r = rank_seq[r_idx]
        l = lambda_seq[l_idx]

        # Training
        if (!is.null(V_train)) {
          V_train_sub = V_train[c(o, covariates),]
          Sigma_hat_train_sub = Sigma_hat_train
        } else {
          V_train_sub = NULL
          Sigma_hat_train_sub = Sigma_hat_train[c(o, covariates), c(o, covariates)]
        }
        fit = matrix_regression_from_cov(Sigma_hat_train_sub, V_train_sub, r, l, ...)

        # Testing
        if (!is.null(V_test)) {
          V_test_sub = V_test[c(o, covariates),]
          Sigma_hat_test_sub = Sigma_hat_test
        } else {
          V_test_sub = NULL
          Sigma_hat_test_sub = Sigma_hat_test[c(o, covariates), c(o, covariates)]
        }

        A = rbind(c(1, rep(0, length(covariates))),
                  c(0, c(
                    tcrossprod(fit$beta1, fit$beta2)
                  )))

        if (!is.null(V_test_sub)) {
          A = A %*% V_test_sub
        }

        Sigma = A %*% Sigma_hat_test_sub %*% t(A)
        r2 = 1 - (Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2]) / Sigma[1, 1]

        list(
          o_idx = o_idx,
          r_idx = r_idx,
          l_idx = l_idx,
          r2 = r2
        )
      }, mc.cores = cores)

      # Accumulate results
      for (res in results_grid) {
        cv_r2[c, res$o_idx, res$r_idx, res$l_idx] =
          cv_r2[c, res$o_idx, res$r_idx, res$l_idx] + res$r2 / K
      }
    }
  }

  # Select best rank/lambda based on max CV R²
  lambda = matrix(NA, length(D_list), length(outcomes))
  rank = matrix(NA, length(D_list), length(outcomes))
  
  for (c in 1:length(D_list)) {
    for (o in 1:length(outcomes)) {
      idxs = which(cv_r2[c, o, , , drop = FALSE] == max(cv_r2[c, o, , ]), arr.ind = TRUE)
      if (length(idxs) > 0) {
        idx = idxs[1, ]
        rank[c, o] = rank_seq[idx[3]]
        lambda[c, o] = lambda_seq[idx[4]]
      } else {
        rank[c, o] = NA
        lambda[c, o] = NA
      }
    }
  }

  return(list(
    lambda = lambda,
    rank = rank,
    cv_r2 = pmax(apply(cv_r2, c(1, 2), max), 0),
    cv_r2_full = cv_r2
  ))

}

cv_raw_matrix_regression = function(Y,
                                    outcomes,
                                    covariates,
                                    rank_seq,
                                    lambda_seq,
                                    K = 2,
                                    folds = NULL,
                                    cores = 8,
                                    ...) {
  require(parallel)
  
  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y) / K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }
  
  # covariates = vech_to_vec_indices(covariates)
  
  cv_r2 = array(0, dim = c(length(outcomes), length(rank_seq), length(lambda_seq)))
  
  for (k in 1:K) {
    Sigma_hat_train = cov2cor_NA0(cov(Y[-folds[[k]],]))
    Sigma_hat_test = cov2cor_NA0(cov(Y[folds[[k]],]))
    

    V_train = V_test = NULL
    
    grid = expand.grid(
      o_idx = seq_along(outcomes),
      r_idx = seq_along(rank_seq),
      l_idx = seq_along(lambda_seq)
    )
    
    results_grid = mclapply(seq_len(nrow(grid)), function(idx) {
      row = grid[idx,]

      o_idx = row$o_idx
      r_idx = row$r_idx
      l_idx = row$l_idx
      
      o = outcomes[o_idx]
      r = rank_seq[r_idx]
      l = lambda_seq[l_idx]
      
      fit = matrix_regression_from_cov(Sigma_hat_train[c(o, covariates), c(o, covariates)], V_train, r, l, ...)
      
      A = rbind(c(1, rep(0, length(covariates))),
                c(0, c(tcrossprod(
                  fit$beta1, fit$beta2
                ))))
      
      if (!is.null(V_test)) {
        A = A %*% V_test
      }
      
      Sigma = A %*% Sigma_hat_test[c(o, covariates), c(o, covariates)] %*% t(A)
      r2 = 1 - (Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2]) / Sigma[1, 1]
      
      list(
        o_idx = o_idx,
        r_idx = r_idx,
        l_idx = l_idx,
        r2 = r2
      )
    }, mc.cores = cores)
    
    for (res in results_grid) {
      cv_r2[res$o_idx, res$r_idx, res$l_idx] =
        cv_r2[res$o_idx, res$r_idx, res$l_idx] + res$r2 / K
    }
    
  }
  
  lambda = numeric(length(outcomes))
  rank = numeric(length(outcomes))

  for (o in 1:length(outcomes)) {
    idxs = which(cv_r2[o, , , drop = FALSE] == max(cv_r2[o, , ]), arr.ind = TRUE)
    if (length(idxs) > 0) {
      idx = idxs[1, ]
      rank[o] = rank_seq[idx[2]]
      lambda[o] = lambda_seq[idx[3]]
    } else {
      rank[o] = NA
      lambda[o] = NA
    }
  }
  
  return(list(
    lambda = lambda,
    rank = rank,
    cv_r2 = pmax(apply(cv_r2, 1, max), 0),
    cv_r2_full = cv_r2
  ))
  
}
