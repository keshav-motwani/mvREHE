library(mvREHE)
library(Matrix)

source("scripts/latent_matrix_regression.R")
source("scripts/latent_ridge_regression.R")
source("scripts/latent_lasso_regression.R")

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
  chol = as.matrix(Matrix::bdiag(replicate(ceiling(n / 1000), chol, simplify = FALSE)))

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
  attr(K_C, "chol") = chol_C

  chol_G = modified_chol(K_G, n)
  K_G = as.matrix(Matrix::bdiag(replicate(ceiling(n / 1000), K_G, simplify = FALSE)))
  attr(K_G, "chol") = chol_G

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

h2_error = function(Sigma_list_estimate, Sigma_list_truth) {

  indices = which(diag(Sigma_list_truth[[1]]) > 1e-14)
  
  h2_estimate = sapply(indices, function(j) Sigma_list_estimate[[2]][j, j] / sum(sapply(Sigma_list_estimate, function(Sigma) Sigma[j, j])))
  h2_truth = sapply(indices, function(j) Sigma_list_truth[[2]][j, j] / sum(sapply(Sigma_list_truth, function(Sigma) Sigma[j, j])))

  sqrt(sum((h2_estimate - h2_truth)^2))

}

max_principal_angle = function(cov_estimate, cov_truth, r) {

  if (!is.null(cov_estimate) & !is.null(cov_truth) & !any(is.na(cov_estimate))) {

    if (ncol(cov_truth) > ncol(cov_estimate)) {
      cov_estimate = expand_estimate(cov_estimate, ncol(cov_truth))
    }

    cor_estimate = cov2cor_NA0(cov_estimate)
    cor_estimate[is.na(cor_estimate)] = 0

    X = eigen(cor_estimate)$vectors[, 1:r, drop = FALSE]
    Y = eigen(cov2cor_NA0(cov_truth))$vectors[, 1:r, drop = FALSE]

    pracma::subspace(X, Y) * 180 / pi

  } else {
    NA
  }

}


get_vec_beta = function(beta1, beta2) {
  beta = tcrossprod(beta1, beta2) + tcrossprod(beta2, beta1) - diag(diag(tcrossprod(beta1, beta2)))
  beta[lower.tri(beta, diag = TRUE)]
}

cov2cor_NA0 = function(cov) {
  cov_hat = cov2cor(cov)
  cov_hat[diag(cov) < .Machine$double.eps, ] = 0
  cov_hat[, diag(cov) < .Machine$double.eps] = 0
  cov_hat
}

r2_beta_hat = function(Y, D_list, fit, estimator, Sigma_true) {
  
  outcomes = setdiff(1:55, c(1, cumsum(10:1)[1:9] + 1))
  covariates = 56:110
  
  Sigma_true = lapply(Sigma_true, cov2cor_NA0)
  
  
  matrix_lambda_rank = cv_latent_matrix_regression(Y, D_list, outcomes, covariates, 1:3, 10^seq(2, -3, length.out = 25), estimator, tolerance = 1e-3)
  ridge_lambda = cv_latent_ridge_regression(Y, D_list, fit, outcomes, covariates, 100, estimator)
  lasso_lambda = cv_latent_lasso_regression(Y, D_list, fit, outcomes, covariates, 100, estimator)
  
  beta_hat_matrix = lapply(1:length(D_list), function(c) {
    sapply(1:length(outcomes), function(o_index) {
      if (matrix_lambda_rank$cv_r2[c, o_index] > 0) {
        reg_fit = latent_matrix_regression(fit, c, outcomes[o_index], covariates, matrix_lambda_rank$rank[c, o_index], matrix_lambda_rank$lambda[c, o_index])
        get_vec_beta(reg_fit$beta1, reg_fit$beta2)
      } else {
        rep(0, length(covariates))
      }
    })
  })
  
  beta_hat_ridge = lapply(1:length(D_list), function(c) {
    sapply(1:length(outcomes), function(o_index) {
      if (ridge_lambda$cv_r2[c, o_index] > 0) {
        latent_ridge_regression(fit, c, outcomes[o_index], covariates, ridge_lambda$lambda[c, o_index])
      } else {
        rep(0, length(covariates))
      }
    })
  })
  
  beta_hat_lasso = lapply(1:length(D_list), function(c) {
    sapply(1:length(outcomes), function(o_index) {
      if (lasso_lambda$cv_r2[c, o_index] > 0) {
        latent_lasso_regression(fit, c, outcomes[o_index], covariates, lasso_lambda$lambda[c, o_index])
      } else {
        rep(0, length(covariates))
      }
    })
  })
  
  r2_matrix = t(sapply(1:length(D_list), function(k) sapply(1:length(outcomes), function(o_index) 1 - (Sigma_true[[k]][outcomes[o_index], outcomes[o_index]] - 2 * t(beta_hat_matrix[[k]][, o_index]) %*% Sigma_true[[k]][covariates, outcomes[o_index]] + t(beta_hat_matrix[[k]][, o_index]) %*% Sigma_true[[k]][covariates, covariates] %*% beta_hat_matrix[[k]][, o_index]) / Sigma_true[[k]][outcomes[o_index], outcomes[o_index]])))
  r2_lasso = t(sapply(1:length(D_list), function(k) sapply(1:length(outcomes), function(o_index) 1 - (Sigma_true[[k]][outcomes[o_index], outcomes[o_index]] - 2 * t(beta_hat_lasso[[k]][, o_index]) %*% Sigma_true[[k]][covariates, outcomes[o_index]] + t(beta_hat_lasso[[k]][, o_index]) %*% Sigma_true[[k]][covariates, covariates] %*% beta_hat_lasso[[k]][, o_index]) / Sigma_true[[k]][outcomes[o_index], outcomes[o_index]])))
  r2_ridge = t(sapply(1:length(D_list), function(k) sapply(1:length(outcomes), function(o_index) 1 - (Sigma_true[[k]][outcomes[o_index], outcomes[o_index]] - 2 * t(beta_hat_ridge[[k]][, o_index]) %*% Sigma_true[[k]][covariates, outcomes[o_index]] + t(beta_hat_ridge[[k]][, o_index]) %*% Sigma_true[[k]][covariates, covariates] %*% beta_hat_ridge[[k]][, o_index]) / Sigma_true[[k]][outcomes[o_index], outcomes[o_index]])))
  
  result = list(matrix = r2_matrix, ridge = r2_ridge, lasso = r2_lasso)
  
  # Initialize empty lists to store values
  method_list <- character()
  component_list <- integer()
  outcome_list <- integer()
  value_list <- numeric()
  
  # Loop over each method
  for (method in names(result)) {
    mat <- result[[method]]
    n_row <- nrow(mat)
    n_col <- ncol(mat)
    
    method_list <- c(method_list, rep(method, n_row * n_col))
    component_list <- c(component_list, rep(1:n_row, times = n_col))
    outcome_list <- c(outcome_list, outcomes[rep(1:n_col, each = n_row)])
    value_list <- c(value_list, c(mat))
  }
  
  # Create data frame
  df <- data.frame(
    regression_method = method_list,
    estimate = paste0("Sigma_", component_list - 1),
    outcome = outcome_list,
    r2 = value_list
  )
  
  attr(df, "cv_r2") = list(matrix_lambda_rank$cv_r2_full, ridge_lambda$cv_r2_full, lasso_lambda$cv_r2_full)
  
  return(df)
  
}

simulation = function(components, n, q, Sigma, method, id, replicate, DATA_ANALYSIS_RESULT_PATH) {

  D_0 = diag(1, nrow = n, ncol = n)
  D_1_and_2 = hcp_kinship(n)
  D_1 = D_1_and_2[[1]]
  D_2 = D_1_and_2[[2]]

  colnames(D_0) = colnames(D_1) = colnames(D_2) = rownames(D_0) = rownames(D_1) = rownames(D_2) = as.character(1:n)

  fit = readRDS(file.path(DATA_ANALYSIS_RESULT_PATH, "fit.rds"))

  heritability_prop = sapply(fit$Sigma_hat, function(x) sum(diag(x)))/sum(sapply(fit$Sigma_hat, function(x) sum(diag(x))))

  if (Sigma != "data") {

    set.seed(123)
    Sigma_0 = heritability_prop[1] * get(paste0("generate_", Sigma, "_Sigma"))(q)
    Sigma_1 = heritability_prop[2] * get(paste0("generate_", Sigma, "_Sigma"))(q)
    Sigma_2 = heritability_prop[3] * get(paste0("generate_", Sigma, "_Sigma"))(q)
    sqrt_Sigma_0 = sqrt(heritability_prop[1]) * attr(Sigma_0, "sqrt")
    sqrt_Sigma_1 = sqrt(heritability_prop[2]) * attr(Sigma_1, "sqrt")
    sqrt_Sigma_2 = sqrt(heritability_prop[3]) * attr(Sigma_2, "sqrt")

  } else {

    Sigma_hat = fit$Sigma_hat
    Sigma_0 = Sigma_hat[[1]]
    Sigma_1 = Sigma_hat[[2]]
    Sigma_2 = Sigma_hat[[3]]
    sqrt_Sigma_0 = attr(Sigma_0, "sqrt")
    sqrt_Sigma_1 = attr(Sigma_1, "sqrt")
    sqrt_Sigma_2 = attr(Sigma_2, "sqrt")
    q = ncol(Sigma_0)

  }

  chol_D_0 = D_0
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
  Epsilon = t(chol_D_0) %*% matrix(rnorm(n * q), nrow = n) %*% t(sqrt_Sigma_0)
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
  } else if (method == "mvREML") {
    source("scripts/mvREML.R")
    estimator = mvREML
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

  if (method == "mvHE") {
    truncated = 1 * sapply(1:length(estimate$Sigma_hat), function(i) attr(estimate$Sigma_hat[[i]], "truncated"))
  } else {
    truncated = rep(0, length(estimate$Sigma_hat))
  }

  min_eigenvalue = sapply(estimate$Sigma_hat, function(Sigma) min(eigen(Sigma)$val))

  if (smoothed) {
    time = system.time({estimate$Sigma_hat = lapply(estimate$Sigma_hat, smooth_cov)})[3] + time
  }

  if (Sigma == "data" & !grepl("REML|cv", method) & grepl("mv", method)) {
    r2_beta_hat = r2_beta_hat(Y, D_list, estimate, estimator, Sigma_list_truth)
  } else {
    r2_beta_hat = NA
  }

  rs = intersect(c(1, 3, 5), 1:(q-1))
  max_principal_angle = sapply(rs, function(r) mapply(max_principal_angle, estimate$Sigma_hat, Sigma_list_truth, r = r))
  rownames(max_principal_angle) = paste0("Sigma_", 1:length(D_list) - 1)
  colnames(max_principal_angle) = rs
  max_principal_angle = reshape2::melt(max_principal_angle, varnames = c("estimate", "r"))

  return(list(
    time = time,
    estimate = estimate,
    Sigma_list_truth = Sigma_list_truth,
    truncated = truncated,
    min_eigenvalue = min_eigenvalue,
    spectral_error = mapply(spectral_error, estimate$Sigma_hat, Sigma_list_truth),
    squared_error = mapply(squared_error, estimate$Sigma_hat, Sigma_list_truth),
    diag_squared_error = mapply(diag_squared_error, estimate$Sigma_hat, Sigma_list_truth),
    max_principal_angle = max_principal_angle,
    h2_error = h2_error(estimate$Sigma_hat, Sigma_list_truth),
    r2_beta_hat = r2_beta_hat
  ))

}
