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

extract_blocks = function(mat) {
  g = igraph::graph.adjacency(mat, weighted = TRUE)
  groups = unique(lapply(Map(sort, igraph::neighborhood(g, nrow(mat))), as.numeric))
  return(groups)
}

make_t_distributed = function(mat, blocks, df) {

  diag(rep(sqrt(df / rchisq(length(blocks), df)), lengths(blocks))) %*% mat * sqrt((df - 2) / df)

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

  h2_estimate = sapply(1:ncol(Sigma_list_estimate[[1]]), function(j) Sigma_list_estimate[[2]][j, j] / sum(sapply(Sigma_list_estimate, function(Sigma) Sigma[j, j])))
  h2_truth = sapply(1:ncol(Sigma_list_truth[[1]]), function(j) Sigma_list_truth[[2]][j, j] / sum(sapply(Sigma_list_estimate, function(Sigma) Sigma[j, j])))

  sqrt(sum((h2_estimate - h2_truth)^2))

}

max_principal_angle = function(cov_estimate, cov_truth, r) {

  if (!is.null(cov_estimate) & !is.null(cov_truth) & !any(is.na(cov_estimate))) {

    if (ncol(cov_truth) > ncol(cov_estimate)) {
      cov_estimate = expand_estimate(cov_estimate, ncol(cov_truth))
    }

    cor_estimate = cov2cor(cov_estimate)
    cor_estimate[is.na(cor_estimate)] = 0

    X = eigen(cor_estimate)$vectors[, 1:r, drop = FALSE]
    Y = eigen(cov2cor(cov_truth))$vectors[, 1:r, drop = FALSE]

    pracma::subspace(X, Y) * 180 / pi

  } else {
    NA
  }

}

pseudoinverse = function(mat) {
  eig = eigen(mat)
  r = sum(eig$val > 1e-12)
  eig$vec[, 1:r] %*% diag(1/eig$val[1:r]) %*% t(eig$vec[, 1:r])
}

beta_error = function(cov_estimate, cov_truth, Y, D_list, component, covariates, outcomes, estimator) {

  cor_estimate = cov2cor(cov_estimate)
  cor_estimate[is.na(cor_estimate)] = 0
  max_eigenvalue = max(eigen(cor_estimate)$values)
  lambda_seq = 10^seq(log10(max_eigenvalue / 1000), log10(max_eigenvalue), length.out = 100)

  lambda = cv_component_ridge_regression(Y, D_list, component, covariates, outcomes, estimator, lambda_seq = lambda_seq)
  beta_estimate = solve(cor_estimate[covariates, covariates] + diag(lambda, length(covariates), length(covariates))) %*% cor_estimate[covariates, outcomes]

  cor_truth = cov2cor(cov_truth)
  beta_truth = solve(cor_truth[covariates, covariates]) %*% cor_truth[covariates, outcomes]

  return(sqrt(sum((beta_estimate - beta_truth)^2)))

}

cv_component_ridge_regression = function(Y, D_list, component, covariates, outcomes, estimator, lambda_seq, K = 2, folds = NULL) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y)/K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }

  cv_loss = numeric(length(lambda_seq))

  for (k in 1:K) {

    fit_train = estimator(Y[-folds[[k]], ], D_list = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]]))
    cor_hat_train = cov2cor(fit_train$Sigma_hat[[component]])
    cor_hat_train[is.na(cor_hat_train)] = 0

    fit_test = estimator(Y[folds[[k]], ], D_list = lapply(D_list, function(D) D[folds[[k]], folds[[k]]]))
    cor_hat_test = cov2cor(fit_test$Sigma_hat[[component]])
    cor_hat_test[is.na(cor_hat_test)] = 0

    for (l in 1:length(lambda_seq)) {

      beta_hat = solve(cor_hat_train[covariates, covariates] + diag(lambda_seq[l], length(covariates), length(covariates))) %*% cor_hat_train[covariates, outcomes]

      cv_loss[l] = cv_loss[l] - 2 * cor_hat_test[outcomes, covariates] %*% beta_hat + t(beta_hat) %*% cor_hat_test[covariates, covariates] %*% beta_hat

    }

  }

  lambda = lambda_seq[which.min(cv_loss)]
  attr(lambda, "cv_loss") = cv_loss

  return(lambda)

}

mvREML_DR5 = function(Y, D_list) {

  R = 5

  PC_Y = prcomp(Y, center = FALSE, scale. = FALSE)
  Y = PC_Y$x[, 1:R]

  estimate = mvREML(Y, D_list)
  estimate$Sigma_hat = lapply(estimate$Sigma_hat, function(x) PC_Y$rotation[, 1:R] %*% x %*% t(PC_Y$rotation[, 1:R]))

  return(estimate)

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

    outcome = 25
    covariates = 92:182

    Sigma_hat = fit$Sigma_hat
    q = ncol(Sigma_hat[[1]])

    cond_num = as.numeric(strsplit(as.character(Sigma), "_")[[1]][2])
    for (k in 1:length(Sigma_hat)) {
      eig = eigen(Sigma_hat[[k]][covariates, covariates])
      diag(Sigma_hat[[k]])[covariates] = diag(Sigma_hat[[k]])[covariates] + eig$val[1] / (cond_num - 1)
      attr(Sigma_hat[[k]], "sqrt") = sqrt_matrix(Sigma_hat[[k]])
    }

    if (length(Sigma_hat) == 2) {
      Sigma_hat[[3]] = NA
    }

    Sigma_0 = Sigma_hat[[1]]
    Sigma_1 = Sigma_hat[[2]]
    Sigma_2 = Sigma_hat[[3]]
    sqrt_Sigma_0 = attr(Sigma_0, "sqrt")
    sqrt_Sigma_1 = attr(Sigma_1, "sqrt")
    sqrt_Sigma_2 = attr(Sigma_2, "sqrt")

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
  }

  time = system.time({estimate = estimator(Y, D_list)})[3]

  if (method == "mvHE") {
    truncated = 1 * sapply(1:length(estimate$Sigma_hat), function(i) attr(estimate$Sigma_hat[[i]], "truncated"))
  } else {
    truncated = rep(0, length(estimate$Sigma_hat))
  }

  if (method == "mvHE") {
    min_eigenvalue = sapply(1:length(estimate$Sigma_hat), function(i) attr(estimate$Sigma_hat[[i]], "min_eigenvalue"))
  } else {
    min_eigenvalue = rep(0, length(estimate$Sigma_hat))
  }

  if (smoothed) {
    time = system.time({estimate$Sigma_hat = lapply(estimate$Sigma_hat, smooth_cov)})[3] + time
  }

  if (grepl("lowdim|data", id) & !grepl("REML|cv", method) & grepl("mv", method)) {
    beta_error = sapply(1:length(estimate$Sigma_hat), function(k) {
      beta_error(estimate$Sigma_hat[[k]], Sigma_list_truth[[k]], Y, D_list, k, covariates, outcome, estimator)
    })
  } else {
    beta_error = NA
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
    beta_error = beta_error
  ))

}
