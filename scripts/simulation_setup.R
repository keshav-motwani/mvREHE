library(mvREHE)
library(Matrix)

source("scripts/latent_regression.R")

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
  chol = as(chol, "dgCMatrix")
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
  K_C = as(K_C, "dsCMatrix")
  K_C = (Matrix::bdiag(replicate(ceiling(n / 1000), K_C, simplify = FALSE)))
  attr(K_C, "chol") = chol_C
  K_C = as(K_C, "dsCMatrix")

  chol_G = modified_chol(K_G, n)
  K_G = as(K_G, "dsCMatrix")
  K_G = (Matrix::bdiag(replicate(ceiling(n / 1000), K_G, simplify = FALSE)))
  attr(K_G, "chol") = chol_G
  K_G = as(K_G, "dsCMatrix")

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



r2_beta_hat = function(Y, D_list, fit, estimator, Sigma_true, method, indices = NULL, fit_full = NULL) {

  if (method == "mvREHE-all") {
    outcomes     = indices[1]
    covariates   = indices[-1]
    fit_for_grid = fit_full
  } else if (ncol(Y) == 110) {
    outcomes     = setdiff(1:55, c(1, cumsum(10:1)[1:9] + 1))
    covariates   = 56:110
    fit_for_grid = fit
  } else {
    outcomes     = 1
    covariates   = 2:ncol(Y)
    fit_for_grid = fit
  }

  outcomes_local   = if (method == "mvREHE-all") 1 else outcomes
  covariates_local = if (method == "mvREHE-all") 2:ncol(Sigma_true[[1]]) else covariates

  r2_from_beta = function(beta_list) {
    do.call(rbind, lapply(seq_along(D_list), function(k)
      sapply(seq_along(outcomes_local), function(o_idx) {
        o   = outcomes_local[o_idx]
        b   = beta_list[[k]][, o_idx]
        Sig = Sigma_true[[k]]
        1 - (Sig[o,o] - 2*t(b) %*% Sig[covariates_local,o] +
               t(b) %*% Sig[covariates_local,covariates_local] %*% b) / Sig[o,o]
      })))
  }

  pinv_beta = function(Sigma_list)
    lapply(Sigma_list, function(Sig)
      sapply(outcomes_local, function(o)
        MASS::ginv(Sig[covariates_local, covariates_local]) %*% Sig[covariates_local, o]))

  r2_true = r2_from_beta(pinv_beta(Sigma_true))

  estimator_fn = function(Y, D_list) estimator(Y, D_list)
  family_ids   = family_ids_from_kinship(D_list[[2]])

  ridge_grid_fn = function(fit, outcomes, covariates)
    ridge_tuning_grid_fn(fit, outcomes, covariates, n_lambda = 20)
  lasso_grid_fn = function(fit, outcomes, covariates)
    lasso_tuning_grid_fn(fit, outcomes, covariates, n_lambda = 20)
  matrix_grid_fn = function(fit, outcomes, covariates)
    matrix_tuning_grid_fn(fit, outcomes, covariates, rank_seq = 1:3,
                          lambda_seq = 10^seq(2, -3, length.out = 5))

  cv_ridge = cv_latent_regression(Y, D_list, family_ids, fit_for_grid, outcomes, covariates,
                                   tuning_grid_fn = ridge_grid_fn, estimator_fn = estimator_fn,
                                   regression_fn = ridge_regression_from_cov, K = 5)
  cv_lasso = if (length(outcomes) > 1) {
    cv_latent_regression(Y, D_list, family_ids, fit_for_grid, outcomes, covariates,
                         tuning_grid_fn = lasso_grid_fn, estimator_fn = estimator_fn,
                         regression_fn = lasso_regression_from_cov, K = 5)
  } else NULL
  cv_matrix = if (length(outcomes) > 1) {
    cv_latent_regression(Y, D_list, family_ids, fit_for_grid, outcomes, covariates,
                         tuning_grid_fn = matrix_grid_fn, estimator_fn = estimator_fn,
                         regression_fn = low_rank_regression_from_cov, K = 5)
  } else NULL

  m = length(outcomes)
  get_beta_cv = function(cv_out, reg_fn) {
    lapply(seq_along(D_list), function(comp) {
      beta = reg_fn(fit_for_grid$Sigma_hat[[comp]], outcomes, covariates,
                    cv_out$tuning_grids[[comp]])
      mat = matrix(0, length(covariates), m)
      for (o in seq_len(m)) {
        sel = cv_out$selected[comp, o]
        if (sel > 0L) mat[, o] = beta[, o, sel]
      }
      mat
    })
  }

  result = c(
    list(ridge = r2_from_beta(get_beta_cv(cv_ridge, ridge_regression_from_cov)) / r2_true),
    if (!is.null(cv_lasso)) list(lasso = r2_from_beta(get_beta_cv(cv_lasso, lasso_regression_from_cov)) / r2_true) else NULL,
    if (!is.null(cv_matrix)) list(matrix = r2_from_beta(get_beta_cv(cv_matrix, low_rank_regression_from_cov)) / r2_true) else NULL
  )

  do.call(rbind, lapply(names(result), function(meth) {
    avg = rowMeans(result[[meth]], na.rm = TRUE)
    data.frame(estimate = paste0("Sigma_", seq_along(avg) - 1), r2 = avg, method = meth)
  }))

}

simulation = function(components, n, q, Sigma, method, commonenv, id, replicate, DATA_ANALYSIS_RESULT_PATH) {

  D_0 = as(as(Matrix::Diagonal(n), "dgCMatrix"), "dsCMatrix")
  D_1_and_2 = hcp_kinship(n)
  D_1 = D_1_and_2[[1]]
  D_2 = D_1_and_2[[2]]

  if (commonenv == "new") {
    D_2 = as(kronecker(Matrix::Diagonal(n = n / 10), tcrossprod(rep(1, 10))), "dsCMatrix")
    attr(D_2, "chol") = D_2 / 2
  }

  colnames(D_0) = colnames(D_1) = colnames(D_2) = rownames(D_0) = rownames(D_1) = rownames(D_2) = as.character(1:n)

  heritability_prop = c(0.80018312, 0.12588264, 0.07393423)

  indices = NULL

  if (Sigma != "data") {

    set.seed(123)
    Sigma_0 = heritability_prop[1] * get(paste0("generate_", Sigma, "_Sigma"))(q)
    Sigma_1 = heritability_prop[2] * get(paste0("generate_", Sigma, "_Sigma"))(q)
    Sigma_2 = heritability_prop[3] * get(paste0("generate_", Sigma, "_Sigma"))(q)
    sqrt_Sigma_0 = sqrt(heritability_prop[1]) * attr(Sigma_0, "sqrt")
    sqrt_Sigma_1 = sqrt(heritability_prop[2]) * attr(Sigma_1, "sqrt")
    sqrt_Sigma_2 = sqrt(heritability_prop[3]) * attr(Sigma_2, "sqrt")

  } else {

    fit = readRDS(file.path(DATA_ANALYSIS_RESULT_PATH, "fit.rds"))
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
  Epsilon = t(chol_D_0) %*% (matrix(rnorm(n * q), nrow = n) %*% t(sqrt_Sigma_0))
  Gamma_1 = t(chol_D_1) %*% (matrix(rnorm(nrow(chol_D_1) * q), nrow = nrow(chol_D_1)) %*% t(sqrt_Sigma_1))
  Gamma_2 = t(chol_D_2) %*% (matrix(rnorm(nrow(chol_D_2) * q), nrow = nrow(chol_D_2)) %*% t(sqrt_Sigma_2))

  Y = Epsilon + Gamma_1
  if (components == 3) {
    Y = Y + Gamma_2
  }

  Y = as.matrix(Y)

  if (Sigma == "lowdim1") {
    Y_full = Y
    Y = Y[, indices]
  }

  if (grepl("smoothed", method)) {
    smoothed = TRUE
    method = gsub("-smoothed", "", method)
  } else {
    smoothed = FALSE
  }

  if (!grepl("HE", method)) {
    D_list = lapply(D_list, as.matrix)
    # for (k in setdiff(1:length(D_list), 1)) {
    #   D_list[[k]] = D_list[[k]] + diag(1e-6, nrow(D_list[[k]]))
    # }
  }

  if (method == "mvHE") {
    estimator = mvHE
  } else if (method == "mvREHE") {
    estimator = mvREHE
  } else if (method == "mvREHE-all") {
    estimator = mvREHE
    Y = Y_full
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

  estimate_full = estimate
  if (method == "mvREHE-all") {
    estimate$Sigma_hat = lapply(estimate$Sigma_hat, function(Sigma) Sigma[indices, indices])
  }

  if (method == "mvHE") {
    truncated = 1 * sapply(1:length(estimate$Sigma_hat), function(i) attr(estimate$Sigma_hat[[i]], "truncated"))
  } else {
    truncated = rep(0, length(estimate$Sigma_hat))
  }

  min_eigenvalue = sapply(estimate$Sigma_hat, function(Sigma) if (any(!is.finite(Sigma))) NA_real_ else min(eigen(Sigma)$val))

  if (smoothed) {
    time = system.time({estimate$Sigma_hat = lapply(estimate$Sigma_hat, smooth_cov)})[3] + time
  }

  if (grepl("mv|GEMMA", method) && (Sigma == "data")) {
    r2_bh   = r2_beta_hat(Y, D_list, estimate, estimator, Sigma_list_truth, method, indices, estimate_full)
  } else {
    r2_bh    = NA
    r2_error = NA
    r2_bias  = NA
  }

  if (!grepl("smooth", Sigma)) {
    h2_error = h2_error(estimate$Sigma_hat, Sigma_list_truth)
  } else {
    h2_error = NA
  }

  return(list(
    time = time,
    estimate = estimate,
    Sigma_list_truth = Sigma_list_truth,
    truncated = truncated,
    min_eigenvalue = min_eigenvalue,
    spectral_error = mapply(spectral_error, estimate$Sigma_hat, Sigma_list_truth),
    squared_error = mapply(squared_error, estimate$Sigma_hat, Sigma_list_truth),
    diag_squared_error = mapply(diag_squared_error, estimate$Sigma_hat, Sigma_list_truth),
    h2_error = h2_error,
    r2_beta_hat = r2_bh
  ))

}
