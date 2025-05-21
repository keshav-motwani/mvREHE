library(mvREHE)
library(Matrix)

expand_estimate = function(estimate, size) {
  kronecker(estimate, tcrossprod(rep(1, size / ncol(estimate))))
}

spectral_error = function(estimate, truth) {
  if (!is.null(estimate) & !is.null(truth)) {
    if (ncol(truth) > ncol(estimate)) {
      estimate = expand_estimate(estimate, ncol(truth))
    }
    norm(estimate - truth, "2")
  } else {
    NA
  }
}

squared_error = function(estimate, truth) {
  if (!is.null(estimate) & !is.null(truth)) {
    if (ncol(truth) > ncol(estimate)) {
      estimate = expand_estimate(estimate, ncol(truth))
    }
    norm(estimate - truth, "F")
  } else {
    NA
  }
}

diag_squared_error = function(estimate, truth) {
  if (!is.null(estimate) & !is.null(truth)) {
    if (ncol(truth) > ncol(estimate)) {
      estimate = expand_estimate(estimate, ncol(truth))
    }
    sqrt(sum(diag(estimate - truth)^2))
  } else {
    NA
  }
}

generate_uniform_Sigma = function(q) {
  matrix = clusterGeneration::rcorrmatrix(q)
  attr(matrix, "sqrt") = sqrt_matrix(matrix)
  matrix
}

generate_fast_Sigma = function(q) {
  V = pracma::randortho(q)
  val = 1/(1:q)^1.25
  matrix = V %*% diag(val) %*% t(V)
  attr(matrix, "sqrt") = V %*% diag(sqrt(val)) %*% t(V)
  matrix
}

generate_moderate_Sigma = function(q) {
  V = pracma::randortho(q)
  val = 1/(1:q)
  matrix = V %*% diag(val) %*% t(V)
  attr(matrix, "sqrt") = V %*% diag(sqrt(val)) %*% t(V)
  matrix
}

generate_slow_Sigma = function(q) {
  V = pracma::randortho(q)
  val = 1/(1:q)^1.25
  matrix = V %*% diag(val) %*% t(V)
  attr(matrix, "sqrt") = V %*% diag(sqrt(val)) %*% t(V)
  matrix
}

generate_smooth_Sigma = function(q, alpha, K = 50) {

  phi_k = function(t,k) cos(pi*k*t)
  lambda_sqrt_k = function(alpha,k) k^(-alpha)
  t_grid = seq(0, 1, length.out = q)

  eigs_sqrt_0 = sapply(1:K, function(kk) lambda_sqrt_k(alpha, kk))
  basis = sapply(1:K, function(kk) phi_k(t_grid,kk))

  matrix = basis %*% diag(eigs_sqrt_0^2) %*% t(basis)
  attr(matrix, "sqrt") = sqrt_matrix(matrix)

  matrix

}

generate_smooth_1_Sigma = function(q) generate_smooth_Sigma(q, 1)
generate_smooth_2_Sigma = function(q) generate_smooth_Sigma(q, 2)

sqrt_matrix = function(A) {
  eig = eigen(A)
  eig$vec %*% diag(sqrt(pmax(eig$val, 0))) %*% t(eig$vec)
}

modified_chol = function(K, n) {

  chol = chol(K[!duplicated(K), !duplicated(K)])
  expand = diag(1, nrow(K), nrow(K))
  expand = expand[, !duplicated(K)]
  expand[cbind(which(duplicated(K)), which(duplicated(K)) - 1:sum(duplicated(K)))] = 1
  chol = chol %*% t(expand)
  chol = (Matrix::bdiag(replicate(ceiling(n / 1000), chol, simplify = FALSE)))

  return(chol)

}

hcp_kinship = function(n) {

  kinship = R.matlab::readMat('data/kinship.mat'); # This is 2*K in the Solar-Eclipse notation
  K_G = as(kinship$K[[1]], "TsparseMatrix"); # Kinship matrix
  K_G = as.matrix(K_G)
  order = hclust(as.dist(-K_G))$order
  K_G = K_G[order, order]
  K_G = K_G[-327, -327] # otherwise 327 is related to two groups of unrelated individuals, which for some reason makes K_C have negative eigenvalues?
  K_G = K_G[1:min(n, 1000), 1:min(n, 1000)]

  K_C = (K_G > 0) * 1
  chol_C = modified_chol(K_C, n)
  K_C = as.matrix(Matrix::bdiag(replicate(ceiling(n / 1000), K_C, simplify = FALSE)))
  attr(K_C, "chol") = as.matrix(chol_C)

  chol_G = modified_chol(K_G, n)
  K_G = as.matrix(Matrix::bdiag(replicate(ceiling(n / 1000), K_G, simplify = FALSE)))
  attr(K_G, "chol") = as.matrix(chol_G)

  return(list(G = K_G, C = K_C))

}

smooth_cov = function(cov, diag = FALSE, output_size = 1000) {

  obsGrid = seq(0, 1, length.out = ncol(cov))
  rcov = list()
  rcov$dataType = "Dense"
  rcov$tPairs = as.matrix(expand.grid(obsGrid, obsGrid))
  if (!diag) {
    indices = rcov$tPairs[, 1] != rcov$tPairs[, 2]
  } else {
    indices = 1:nrow(rcov$tPairs)
  }
  rcov$tPairs = rcov$tPairs[indices, ]
  rcov$cxxn = c(cov)[indices]

  gcvObj = fdapace:::GCVLwls2DV2(obsGrid, obsGrid, kern = "epan", rcov = rcov, t = list(obsGrid))
  bwCov = gcvObj$h
  out = seq(0, 1, length.out = output_size)
  smoothCov = fdapace:::Lwls2D(bwCov, "epan", xin=rcov$tPairs, yin=rcov$cxxn,
                               xout1=out, xout2=out)

  smoothCov

}

matrix_regression_from_cov = function(Sigma_r_hat, V, rank, lambda, max_iter = 1000, tolerance = 1e-6, n_init = 5) {

  if (!is.null(V)) {
    p = sqrt(nrow(V) - 1)
  } else {
    p = sqrt(nrow(Sigma_r_hat) - 1)
  }
  best_objective = Inf
  best_solution = NULL

  for (init in 1:n_init) {

    print(init)

    if (length(p) == 0) {
      print("huh")
    }
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
      Sigma = A %*% Sigma_r_hat %*% t(A)
      beta1[] = solve(Sigma[-1, -1] + diag(lambda, ncol(Sigma) - 1, ncol(Sigma) - 1), Sigma[1, -1])

      if (!is.null(V)) {
        # A = compute_A_for_objective(beta1, beta2, V)
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2)))) %*% V
      } else {
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2))))
      }
      Sigma = A %*% Sigma_r_hat %*% t(A)
      objective[2 * iter - 1] = Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2] + lambda * (sum(beta1^2) + sum(beta2^2))

      if (!is.null(V)) {
        # A = compute_A_for_Sigma_y_Xtbeta1(beta1, V)
        A = rbind(c(1, rep(0, p^2)), cbind(rep(0, p * rank), kronecker(diag(1, p, p), t(beta1)))) %*% V
      } else {
        A = rbind(c(1, rep(0, p^2)), cbind(rep(0, p * rank), kronecker(diag(1, p, p), t(beta1))))
      }
      Sigma = A %*% Sigma_r_hat %*% t(A)
      beta2[] = matrix(solve(Sigma[-1, -1] + diag(lambda, ncol(Sigma) - 1, ncol(Sigma) - 1), Sigma[1, -1]), ncol = rank, byrow = TRUE)

      if (!is.null(V)) {
        # A = compute_A_for_objective(beta1, beta2, V)
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2)))) %*% V
      } else {
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2))))
      }
      Sigma = A %*% Sigma_r_hat %*% t(A)
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

latent_matrix_regression = function(Y, D_list, component, outcome, covariates, rank, lambda, estimator, ...) {

  fit = estimator(Y, D_list)

  Sigma_r_hat = fit$Sigma_hat[[component]]

  covariates = vech_to_vec_indices(covariates)

  if (!is.null(fit$V)) {
    V <- fit$V[c(outcome, covariates), ]
  } else {
    V = NULL
    Sigma_r_hat = Sigma_r_hat[c(outcome, covariates), c(outcome, covariates)]
  }

  matrix_regression_from_cov(Sigma_r_hat, V, rank, lambda, ...)

}


matrix_regression = function(Y, outcome, covariates, rank, lambda, ...) {

  Sigma_r_hat = cov(Y)

  covariates = vech_to_vec_indices(covariates)

  V = NULL
  Sigma_r_hat = Sigma_r_hat[c(outcome, covariates), c(outcome, covariates)]

  matrix_regression_from_cov(Sigma_r_hat, V, rank, lambda, ...)

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
                                        K = 2, folds = NULL, cores = parallel::detectCores() - 1, ...) {

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

    fit_train = estimator(Y[-folds[[k]], ], D_list = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]]), return_full = FALSE)
    fit_test = estimator(Y[folds[[k]], ], D_list = lapply(D_list, function(D) D[folds[[k]], folds[[k]]]), return_full = FALSE)

    for (c in 1:length(D_list)) {

      # Create full grid of outcome × rank × lambda
      grid <- expand.grid(
        o_idx = seq_along(outcomes),
        r_idx = seq_along(rank_seq),
        l_idx = seq_along(lambda_seq)
      )

      results_grid <- lapply(seq_len(nrow(grid)), function(idx) {
        row <- grid[idx, ]
        o_idx <- row$o_idx
        r_idx <- row$r_idx
        l_idx <- row$l_idx

        o <- outcomes[o_idx]
        r <- rank_seq[r_idx]
        l <- lambda_seq[l_idx]

        # Training
        Sigma_r_hat_train <- fit_train$Sigma_r_hat[[c]]
        if (!is.null(fit_train$V)) {
          V_train <- fit_train$V[c(o, covariates), ]
        } else {
          V_train = NULL
          Sigma_r_hat_train = Sigma_r_hat_train[c(o, covariates), c(o, covariates)]
        }
        fit <- matrix_regression_from_cov(Sigma_r_hat_train, V_train, r, l, ...)
        # Testing
        Sigma_r_hat_test <- fit_test$Sigma_r_hat[[c]]
        if (!is.null(fit_test$V)) {
          V_test <- fit_test$V[c(o, covariates), ]
        } else {
          V_test = NULL
          Sigma_r_hat_test = Sigma_r_hat_test[c(o, covariates), c(o, covariates)]
        }

        A <- rbind(
          c(1, rep(0, length(covariates))),
          c(0, c(tcrossprod(fit$beta1, fit$beta2)))
        )

        if (!is.null(V_test)) {
          A = A %*% V_test
        }

        Sigma <- A %*% Sigma_r_hat_test %*% t(A)
        r2 <- 1 - (Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2]) / Sigma[1, 1]

        list(o_idx = o_idx, r_idx = r_idx, l_idx = l_idx, r2 = r2 / K)
      })

      # Accumulate results
      for (res in results_grid) {
        cv_r2[c, res$o_idx, res$r_idx, res$l_idx] <-
          cv_r2[c, res$o_idx, res$r_idx, res$l_idx] + res$r2
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

get_vec_beta = function(beta1, beta2) {
  beta = tcrossprod(beta1, beta2) + tcrossprod(beta2, beta1) - diag(diag(tcrossprod(beta1, beta2)))
  beta[lower.tri(beta, diag = TRUE)]
}

beta_error = function(Y, D_list, Sigma_hat, Sigma_r_hat, V, estimator, Sigma_true, beta_true) {

  outcome = 1
  covariates = 2:56

  cv_fit = cv_latent_matrix_regression(Y, D_list, outcome, covariates, 1:3, 10^seq(0, -2, length.out = 10), estimator, 2)

  beta_hat = lapply(1:3, function(k) {
    covariates = vech_to_vec_indices(covariates)
    if (!is.null(V)) {
      V <- V[c(outcome, covariates), ]
      Sigma = Sigma_r_hat[[k]]
    } else {
      Sigma = Sigma_r_hat[[k]][c(outcome, covariates), c(outcome, covariates)]
    }
    fit = matrix_regression_from_cov(Sigma, V, cv_fit$rank[k, 1], cv_fit$lambda[k, 1])
    get_vec_beta(fit$beta1, fit$beta2)
  })

  beta_error = sapply(1:3, function(k) sqrt(sum((beta_hat[[k]] - beta_true[[k]])^2)))

  r2 = sapply(1:3, function(k) 1 - (Sigma_hat[[k]][outcome, outcome] - 2 * Sigma_hat[[k]][outcome, covariates] %*% beta_true[[k]] + t(beta_true[[k]]) %*% Sigma_hat[[k]][covariates, covariates] %*% beta_true[[k]]) / Sigma_hat[[k]][outcome, outcome])
  r2_true = sapply(1:3, function(k) 1 - (Sigma_true[[k]][outcome, outcome] - 2 * Sigma_true[[k]][outcome, covariates] %*% beta_true[[k]] + t(beta_true[[k]]) %*% Sigma_true[[k]][covariates, covariates] %*% beta_true[[k]]) / Sigma_true[[k]][outcome, outcome])

  r2_error = abs(r2 - r2_true)

  return(rbind(beta_error, r2_error))

}

simulation = function(components, n, q, Sigma, method, id, replicate, DATA_ANALYSIS_RESULT_PATH) {

  D_0 = diag(1, nrow = n, ncol = n)
  D_1_and_2 = hcp_kinship(n)
  D_1 = D_1_and_2[[1]]
  D_2 = D_1_and_2[[2]]

  colnames(D_0) = colnames(D_1) = colnames(D_2) = rownames(D_0) = rownames(D_1) = rownames(D_2) = as.character(1:n)

  fit = readRDS(file.path(DATA_ANALYSIS_RESULT_PATH, "fit.rds"))

  heritability_prop = sapply(fit$Sigma_hat, function(x) sum(diag(x)))/sum(sapply(fit$Sigma_hat, function(x) sum(diag(x))))

  if (!grepl("data", Sigma)) {

    outcome = 1
    covariates = setdiff(1:q, outcome)

    set.seed(123)
    Sigma_0 = heritability_prop[1] * get(paste0("generate_", Sigma, "_Sigma"))(q)
    Sigma_1 = heritability_prop[2] * get(paste0("generate_", Sigma, "_Sigma"))(q)
    Sigma_2 = heritability_prop[3] * get(paste0("generate_", Sigma, "_Sigma"))(q)
    sqrt_Sigma_0 = sqrt(heritability_prop[1]) * attr(Sigma_0, "sqrt")
    sqrt_Sigma_1 = sqrt(heritability_prop[2]) * attr(Sigma_1, "sqrt")
    sqrt_Sigma_2 = sqrt(heritability_prop[3]) * attr(Sigma_2, "sqrt")

  } else {

    Sigma_hat = fit$Sigma_hat_sim

    if (length(Sigma_hat) == 2) {
      Sigma_hat[[3]] = NA
    }

    Sigma_0 = Sigma_hat[[1]]
    Sigma_1 = Sigma_hat[[2]]
    Sigma_2 = Sigma_hat[[3]]
    sqrt_Sigma_0 = attr(Sigma_0, "sqrt")
    sqrt_Sigma_1 = attr(Sigma_1, "sqrt")
    sqrt_Sigma_2 = attr(Sigma_2, "sqrt")

    q = ncol(Sigma_hat[[1]])

  }

  chol_D_1 = attr(D_1, "chol")
  chol_D_2 = attr(D_2, "chol")

  if (id == "smooth") {
    sqrt_Sigma_0 = sqrt_matrix(Sigma_0 + diag(1, q, q))
    set.seed(123)
    Sigma_0 = heritability_prop[1] * get(paste0("generate_", Sigma, "_Sigma"))(1000)
    Sigma_1 = heritability_prop[2] * get(paste0("generate_", Sigma, "_Sigma"))(1000)
    Sigma_2 = heritability_prop[3] * get(paste0("generate_", Sigma, "_Sigma"))(1000)
  }

  if (components == 2) {
    Sigma_list_truth = list(Sigma_0, Sigma_1)
    D_list = list(D_0, D_1)
  } else if (components == 3) {
    Sigma_list_truth = list(Sigma_0, Sigma_1, Sigma_2)
    D_list = list(D_0, D_1, D_2)
  }

  set.seed(replicate)
  Epsilon = matrix(rnorm(n * q), nrow = n) %*% t(sqrt_Sigma_0)
  Gamma_1 = t(chol_D_1) %*% matrix(rnorm(nrow(chol_D_1) * q), nrow = nrow(chol_D_1)) %*% t(sqrt_Sigma_1)
  Gamma_2 = t(chol_D_2) %*% matrix(rnorm(nrow(chol_D_2) * q), nrow = nrow(chol_D_2)) %*% t(sqrt_Sigma_2)

  Y = Epsilon + Gamma_1
  if (components == 3) {
    Y = Y + Gamma_2
  }

  if (grepl("smoothed", method)) {
    smoothed = TRUE
    method = gsub("-smoothed", "", method)
  } else {
    smoothed = FALSE
  }

  if (method == "mvHE") {
    estimator = mvHE
  } else if (method == "mvREHE") {
    estimator = mvREHE
  } else if (method == "mvREHE_cvDR") {
    estimator = function(Y, D_list) {
      fit = mvREHE_cvDR(Y, D_list, K = 5, r_seq = floor(seq(5, q, length.out = 20)), compute_full_Sigma = TRUE)
      fit$Sigma_hat = lapply(fit$Sigma_r_hat, function(Sigma) fit$V %*% Sigma %*% t(fit$V))
      fit
    }
  } else if (method == "mvREML") {
    source("scripts/mvREML.R")
    estimator = mvREML
  }
  else if (method == "mvREML_DR5") {
    source("scripts/mvREML.R")
    estimator = mvREML_DR5
  } else if (method == "HE") {
    estimator = function(Y, D_list) {
      univariate(Y, D_list, mvHE)
    }
  } else if (method == "REHE") {
    estimator = function(Y, D_list) {
      univariate(Y, D_list, mvREHE)
    }
  } else if (method == "REML") {
    source("scripts/mvREML.R")
    estimator = function(Y, D_list) {
      univariate(Y, D_list, mvREML)
    }
  } else if (method == "GEMMA") {
    source("scripts/GEMMA.R")
    estimator = GEMMA
  }

  time = system.time({estimate = estimator(Y, D_list)})[3]

  # if (method == "mvHE") {
  #   truncated = 1 * sapply(1:length(estimate$Sigma_hat), function(i) attr(estimate$Sigma_hat[[i]], "truncated"))
  # } else {
  #   truncated = rep(0, length(estimate$Sigma_hat))
  # }

  min_eigenvalue = NA # sapply(estimate$Sigma_hat, function(Sigma) min(eigen(Sigma)$val))

  if (smoothed) {
    time = system.time({estimate$Sigma_hat = lapply(estimate$Sigma_hat, smooth_cov)})[3] + time
  }

  if (grepl("data", id) & !grepl("REML|cv", method) & grepl("mv", method)) {
    beta_error = beta_error(Y, D_list, estimate$Sigma_hat, estimate$Sigma_r_hat, estimate$V, estimator, Sigma_list_truth, fit$beta)
  } else {
    beta_error = rbind(NA, NA)
  }

  # rs = intersect(c(1, 3, 5), 1:(q-1))
  # max_principal_angle = sapply(rs, function(r) mapply(max_principal_angle, estimate$Sigma_hat, Sigma_list_truth, r = r))
  # rownames(max_principal_angle) = paste0("Sigma_", 1:length(D_list) - 1)
  # colnames(max_principal_angle) = rs
  # max_principal_angle = reshape2::melt(max_principal_angle, varnames = c("estimate", "r"))

  return(list(
    time = time,
    estimate = estimate,
    Sigma_list_truth = Sigma_list_truth,
    # truncated = truncated,
    min_eigenvalue = min_eigenvalue,
    spectral_error = mapply(spectral_error, estimate$Sigma_hat, Sigma_list_truth),
    squared_error = mapply(squared_error, estimate$Sigma_hat, Sigma_list_truth),
    diag_squared_error = mapply(diag_squared_error, estimate$Sigma_hat, Sigma_list_truth),
    # max_principal_angle = max_principal_angle,
    # h2_error = h2_error(estimate$Sigma_hat, Sigma_list_truth),
    beta_error = beta_error[1, ],
    r2_error = beta_error[2, ]
  ))

}
