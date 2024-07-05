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

univariate = function(Y, D_list, method_fn) {

  q = ncol(Y)

  estimates = apply(Y, 2, FUN = method_fn, D_list = D_list, simplify = FALSE)

  Sigma_hat = lapply(1:length(D_list), function(i) matrix(0, q, q))

  for (j in 1:q) {
    for (k in 1:length(D_list)) {
      Sigma_hat[[k]][j, j] = estimates[[j]]$Sigma_hat[[k]][1, 1]
    }
  }

  return(list(Sigma_hat = Sigma_hat))

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

var_prop = function(Sigma_list) {

  sapply(Sigma_list, function(x) sum(diag(x))) / sum(sapply(Sigma_list, function(x) sum(diag(x))))

}

max_principal_angle = function(cov_estimate, cov_truth, r) {

  if (!is.null(cov_estimate) & !is.null(cov_truth)) {

    if (ncol(cov_truth) > ncol(cov_estimate)) {
      cov_estimate = expand_cov_estimate(cov_estimate, ncol(cov_truth))
    }

    X = eigen(cov_estimate)$vectors[, 1:r, drop = FALSE]
    Y = eigen(cov_truth)$vectors[, 1:r, drop = FALSE]

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
  lambda_seq = seq(max_eigenvalue / 1000, max_eigenvalue, length.out = 10)

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

  for (l in 1:length(lambda_seq)) {

    for (k in 1:K) {

      # fit_train = mvREHE::mvREHE_cvDR(Y[-folds[[k]], ], D_list = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]]), r_seq = 1:4 * 10, V_function = separate_svd_irlba, tolerance = tolerance, max_iter = max_iter)
      # fit_train$Sigma_hat = lapply(fit_train$Sigma_r_hat, function(Sigma) fit_train$V %*% Sigma %*% t(fit_train$V))
      fit_train = estimator(Y[-folds[[k]], ], D_list = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]]))
      cor_hat_train = cov2cor(fit_train$Sigma_hat[[component]])
      cor_hat_train[is.na(cor_hat_train)] = 0
      beta_hat = solve(cor_hat_train[covariates, covariates] + diag(lambda_seq[l], length(covariates), length(covariates))) %*% cor_hat_train[covariates, outcomes]

      # fit_test = mvREHE::mvREHE_cvDR(Y[folds[[k]], ], D_list = lapply(D_list, function(D) D[folds[[k]], folds[[k]]]),  r_seq = 1:4 * 10, V_function = separate_svd_irlba, tolerance = tolerance, max_iter = max_iter)
      # fit_test$Sigma_hat = lapply(fit_test$Sigma_r_hat, function(Sigma) fit_test$V %*% Sigma %*% t(fit_test$V))
      fit_test = estimator(Y[folds[[k]], ], D_list = lapply(D_list, function(D) D[folds[[k]], folds[[k]]]))
      cor_hat_test = cov2cor(fit_test$Sigma_hat[[component]])
      cor_hat_test[is.na(cor_hat_test)] = 0
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

  estimate = mvREML_2(Y, D_list)
  estimate$Sigma_hat = lapply(estimate$Sigma_hat, function(x) PC_Y$rotation[, 1:R] %*% x %*% t(PC_Y$rotation[, 1:R]))

  return(estimate)

}

simulation = function(n, q, Sigma, method, id, replicate) {

  D_0 = diag(1, nrow = n, ncol = n)
  D_1_and_2 = hcp_kinship(n)
  D_1 = D_1_and_2[[1]]

  colnames(D_0) = colnames(D_1) = rownames(D_0) = rownames(D_1) = as.character(1:n)

  fit = readRDS(file.path(DATA_ANALYSIS_RESULT_PATH, "fit.rds"))

  heritability_prop = sapply(fit$Sigma_hat, function(x) sum(diag(x)))/sum(sapply(fit$Sigma_hat, function(x) sum(diag(x))))

  if (Sigma != "data") {

    set.seed(123)
    Sigma_0 = heritability_prop[1] * get(paste0("generate_", Sigma, "_Sigma"))(q)
    Sigma_1 = heritability_prop[2] * get(paste0("generate_", Sigma, "_Sigma"))(q)
    sqrt_Sigma_0 = sqrt(heritability_prop[1]) * attr(Sigma_0, "sqrt")
    sqrt_Sigma_1 = sqrt(heritability_prop[2]) * attr(Sigma_1, "sqrt")

  } else if (id == 4) {

    outcome = 9
    covariates = 92:182

    Sigma_hat = fit$Sigma_hat
    q = ncol(fit$Sigma_hat)

    for (k in 1:length(Sigma_hat)) {
      eig = eigen(Sigma_hat[[k]][covariates, covariates])
      diag(Sigma_hat[[k]])[covariates] = diag(Sigma_hat[[k]])[covariates] + eig$val[1] / (cond_num - 1)
      attr(Sigma_hat[[k]], "sqrt") = sqrt_matrix(Sigma_hat[[k]])
    }

    Sigma_0 = Sigma_hat[[1]]
    Sigma_1 = Sigma_hat[[2]]
    sqrt_Sigma_0 = attr(Sigma_0, "sqrt")
    sqrt_Sigma_1 = attr(Sigma_1, "sqrt")

  }

  chol_D_0 = D_0
  chol_D_1 = attr(D_1, "chol")

  Sigma_list_truth = list(Sigma_0, Sigma_1)
  D_list = list(D_0, D_1)

  if (grepl("smooth", Sigma)) {
    Sigma_0 = Sigma_0 + diag(1, q, q)
    sqrt_Sigma_0 = sqrt_matrix(Sigma_0)
  }

  set.seed(replicate)
  Epsilon = t(chol_D_0) %*% matrix(rnorm(n * q), nrow = n) %*% t(sqrt_Sigma_0)
  Gamma_1 = t(chol_D_1) %*% matrix(rnorm(nrow(chol_D_1) * q), nrow = nrow(chol_D_1)) %*% t(sqrt_Sigma_1)
  Gamma_2 = t(chol_D_2) %*% matrix(rnorm(nrow(chol_D_2) * q), nrow = nrow(chol_D_2)) %*% t(sqrt_Sigma_2)
  Y = Epsilon + Gamma_1 + Gamma_2

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
      fit = mvREHE_cvDR(Y, D_list, K = 5, r_seq = intersect(1:q, c(5, 10, 50, 100, q)), compute_full_Sigma = TRUE)
      fit$Sigma_hat = lapply(fit$Sigma_r_hat, function(Sigma) fit$V %*% Sigma %*% t(fit$V))
      fit
    }
  } else if (method == "mvREML") {
    source("scripts/mvREML.R")
    estimator = mvREML_2
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
      univariate(Y, D_list, mvREML_2)
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

  if (Sigma == "data" & !grepl("REML", method) & grepl("mv", method)) {
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
    h2_error = (var_prop(estimate$Sigma_hat)[2] - var_prop(Sigma_list_truth)[2])^2,
    beta_error = beta_error
  ))

}

SIMULATION_ID = as.numeric(commandArgs(trailingOnly=TRUE)[1])

RESULT_PATH = paste0("simulation_2_components_", SIMULATION_ID)
DATA_ANALYSIS_RESULT_PATH = "data_analysis_2_components"
dir.create(RESULT_PATH, recursive = TRUE)

replicates = 1:50

if (SIMULATION_ID == "lowdim1") { # 5700
  methods = c("mvHE", "mvREHE", "mvREML_DR5", "HE", "REHE", "REML")
  Sigmas = "uniform"
  ns = c(250, 500, 1000, 2000, 4000, 8000)
  qs = c(10, 20)
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n")
  qs = 5
  grid = rbind(grid, expand.grid(method = c(methods, "mvREML"), replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n"))
} else if (SIMULATION_ID == "lowdim2") { # 5700
  methods = c("mvHE", "mvREHE", "mvREML_DR5", "HE", "REHE", "REML")
  Sigmas = "moderate"
  ns = c(250, 500, 1000, 2000, 4000, 8000)
  qs = c(10, 20)
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n")
  qs = 5
  grid = rbind(grid, expand.grid(method = c(methods, "mvREML"), replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n"))
} else if (SIMULATION_ID == "highdim") { # 3000
  methods = c("mvHE", "mvREHE", "mvREHE_cvDR", "mvREML_DR5")
  Sigmas = c("fast", "moderate", "slow")
  ns = c(500, 1000, 2000, 4000, 8000)
  qs = 1000
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n")
} else if (SIMULATION_ID == "smooth") { # 2000
  methods = c("mvHE", "mvREHE", "mvHE-smoothed", "mvREHE-smoothed")
  Sigmas = c("smooth_1", "smooth_2")
  qs = 100
  ns = c(500, 1000, 2000, 4000, 8000)
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n")
} else if (SIMULATION_ID == "data") { # 1000
  methods = c("mvHE", "mvREHE", "mvREHE_cvDR", "mvREML_DR5")
  Sigmas = "data"
  ns = c(500, 1000, 2000, 4000, 8000)
  qs = NA
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n")
}

PARAMETER_ID = as.numeric(commandArgs(trailingOnly=TRUE)[2])
print(grid[PARAMETER_ID, ])
replicate = grid[PARAMETER_ID, "replicate"]
n = grid[PARAMETER_ID, "n"]
q = grid[PARAMETER_ID, "q"]
Sigma = grid[PARAMETER_ID, "Sigma"]
method = as.character(grid[PARAMETER_ID, "method"])
experiment = grid[PARAMETER_ID, "experiment"]
cond_num = 100

output = simulation(n, q, Sigma, method, SIMULATION_ID, replicate)
estimate = paste0("Sigma_", 1:2 - 1)

diag_squared_error = data.frame(replicate = replicate, estimate = estimate, diag_squared_error = output$diag_squared_error, n = n, q = q, Sigma = Sigma, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
squared_error = data.frame(replicate = replicate, estimate = estimate, squared_error = output$squared_error, n = n, q = q, Sigma = Sigma, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
spectral_error = data.frame(replicate = replicate, estimate = estimate, spectral_error = output$spectral_error, n = n, q = q, Sigma = Sigma, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
h2_error = data.frame(replicate = replicate, h2_error = output$h2_error, n = n, q = q, Sigma = Sigma, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
beta_error = data.frame(replicate = replicate, estimate = estimate, beta_error = output$beta_error, n = n, q = q, Sigma = Sigma, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
max_principal_angle = cbind(output$max_principal_angle, replicate = replicate, n = n, q = q, Sigma = Sigma, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
time = data.frame(replicate = replicate, method = method, time = output$time, n = n, q = q, Sigma = Sigma, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
truncated = data.frame(estimate = estimate, truncated = output$truncated, n = n, q = q, Sigma = Sigma, method = method, replicate = replicate, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
min_eigenvalue = data.frame(estimate = estimate, min_eigenvalue = output$min_eigenvalue, n = n, q = q, Sigma = Sigma, method = method, replicate = replicate, experiment = experiment, SIMULATION_ID = SIMULATION_ID)

saveRDS(list(diag_squared_error = diag_squared_error, squared_error = squared_error, spectral_error = spectral_error, h2_error = h2_error, beta_error = beta_error, max_principal_angle = max_principal_angle, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue), file.path(RESULT_PATH, paste0("n", n, "_q", q, "_Sigma", Sigma, "_replicate", replicate, "_experiment", experiment, "_method", method, ".rds")))
print(list(diag_squared_error = diag_squared_error, squared_error = squared_error, spectral_error = spectral_error, h2_error = h2_error, beta_error = beta_error, max_principal_angle = max_principal_angle, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue))
