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
  if (q <= 100) {
    matrix = clusterGeneration::rcorrmatrix(q)
    sqrt_matrix = sqrt_matrix(matrix)
  } else {
    blocks = lapply(1:(q / 100), function(x) clusterGeneration::rcorrmatrix(100))
    sqrt_blocks = lapply(blocks, sqrt_matrix)
    matrix = as.matrix(Matrix::bdiag(blocks))
    sqrt_matrix = as.matrix(Matrix::bdiag(sqrt_blocks))
  }
  attr(matrix, "sqrt") = sqrt_matrix
  matrix
}

generate_fast_Sigma = function(q) {
  V = pracma::randortho(q)
  matrix = V %*% diag(1/(1:q)^1.25) %*% t(V)
  attr(matrix, "sqrt") = V %*% diag(1/(1:q)^(1.25/2)) %*% t(V)
  matrix
}

generate_moderate_Sigma = function(q) {
  V = pracma::randortho(q)
  matrix = V %*% diag(1/(1:q)^1) %*% t(V)
  attr(matrix, "sqrt") = V %*% diag(1/(1:q)^0.5) %*% t(V)
  matrix
}

generate_slow_Sigma = function(q) {
  V = pracma::randortho(q)
  matrix = V %*% diag(1/(1:q)^0.75) %*% t(V)
  attr(matrix, "sqrt") = V %*% diag(1/(1:q)^(0.75/2)) %*% t(V)
  matrix
}

generate_constant_Sigma = function(q) {
  matrix = diag(1, q, q)
  attr(matrix, "sqrt") = matrix
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

ar1_cor = function(n, m, rho) {

  exponent = abs(matrix(1:m - 1, nrow = m, ncol = m, byrow = TRUE) - (1:m - 1))
  mat = rho^exponent

  chol = as.matrix(Matrix::bdiag(replicate(n / m, chol(mat), simplify = FALSE)))
  mat = as.matrix(Matrix::bdiag(replicate(n / m, mat, simplify = FALSE)))
  attr(mat, "chol") = chol

  mat

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

h2_prop = function(Sigma_list, G_index) {

  sum(diag(Sigma_list[[G_index]])) / sum(sapply(Sigma_list, function(x) sum(diag(x))))

}

max_principal_angle = function(estimate, truth, r) {

  if (!is.null(estimate) & !is.null(truth)) {

    if (ncol(truth) > ncol(estimate)) {
      estimate = expand_estimate(estimate, ncol(truth))
    }

    X = eigen(estimate)$vectors[, 1:r, drop = FALSE]
    Y = eigen(truth)$vectors[, 1:r, drop = FALSE]

    pracma::subspace(X, Y) * 180 / pi

  } else {
    NA
  }

}

simulation = function(n, q, Sigma, method, replicate) {

  D_0 = diag(1, nrow = n, ncol = n)
  D_1_and_2 = hcp_kinship(n)
  D_1 = D_1_and_2[[1]]
  D_2 = D_1_and_2[[2]]

  colnames(D_0) = colnames(D_1) = colnames(D_2) = rownames(D_0) = rownames(D_1) = rownames(D_2) = as.character(1:n)

  heritability_prop = c(0.8, 0.1, 0.1)

  set.seed(123)
  Sigma_0 = heritability_prop[1] * get(paste0("generate_", Sigma, "_Sigma"))(q)
  Sigma_1 = heritability_prop[2] * get(paste0("generate_", Sigma, "_Sigma"))(q)
  Sigma_2 = heritability_prop[3] * get(paste0("generate_", Sigma, "_Sigma"))(q)
  sqrt_Sigma_0 = sqrt(heritability_prop[1]) * attr(Sigma_0, "sqrt")
  sqrt_Sigma_1 = sqrt(heritability_prop[2]) * attr(Sigma_1, "sqrt")
  sqrt_Sigma_2 = sqrt(heritability_prop[3]) * attr(Sigma_2, "sqrt")

  chol_D_0 = D_0
  chol_D_1 = attr(D_1, "chol")
  chol_D_2 = attr(D_2, "chol")

  if (grepl("smooth", Sigma)) {
   Sigma_0 = Sigma_0 + diag(1, q, q)
   sqrt_Sigma_0 = sqrt_matrix(Sigma_0)
  }

  set.seed(replicate)
  Epsilon = t(chol_D_0) %*% matrix(rnorm(n * q), nrow = n) %*% t(sqrt_Sigma_0)
  Gamma_1 = t(chol_D_1) %*% matrix(rnorm(nrow(chol_D_1) * q), nrow = nrow(chol_D_1)) %*% t(sqrt_Sigma_1)
  Gamma_2 = t(chol_D_2) %*% matrix(rnorm(nrow(chol_D_2) * q), nrow = nrow(chol_D_2)) %*% t(sqrt_Sigma_2)
  Y = Epsilon + Gamma_1 + Gamma_2

  if (grepl("-", method) && substring(method, 1, 2) == "DR") {
    PC = TRUE
    PC_Y = prcomp(Y, center = FALSE, scale. = FALSE)
    R = as.numeric(gsub("DR", "", strsplit(method, "-")[[1]][1]))
    Y = PC_Y$x[, 1:R]
    method = strsplit(method, "-")[[1]][2]
  } else {
    PC = FALSE
  }

  if (grepl("smoothed", method)) {
    smoothed = TRUE
    method = gsub("-smoothed", "", method)
  } else {
    smoothed = FALSE
  }

  D_list = list(D_0, D_1, D_2)

  if (method == "mvHE") {
    time = system.time({estimate = mvHE(Y, D_list)})[3]
  } else if (method == "mvREHE") {
    time = system.time({estimate = mvREHE(Y, D_list)})[3]
  } else if (method == "mvREHE_cvDR") {
    time = system.time({estimate = mvREHE_cvDR(Y, D_list, K = 5, V_seq = intersect(1:q, c(5, 10, 50, 100, q)))})[3]
  } else if (method == "mvREML") {
    source("scripts/mvREML.R")
    time = system.time({estimate = mvREML(Y, D_list)})[3]
  } else if (method == "HE") {
    time = system.time({estimate = univariate(Y, D_list, mvHE)})[3]
  } else if (method == "REHE") {
    time = system.time({estimate = univariate(Y, D_list, mvREHE)})[3]
  } else if (method == "REML") {
    source("scripts/mvREML.R")
    time = system.time({estimate = univariate(Y, D_list, mvREML)})[3]
  }

  if (PC) {
    estimate$Sigma_hat = lapply(estimate$Sigma_hat, function(x) PC_Y$rotation[, 1:R] %*% x %*% t(PC_Y$rotation[, 1:R]))
  }

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

  if (grepl("smooth", Sigma)) {
    set.seed(123)
    Sigma_0 = heritability_prop[1] * get(paste0("generate_", Sigma, "_Sigma"))(1000)
    Sigma_1 = heritability_prop[2] * get(paste0("generate_", Sigma, "_Sigma"))(1000)
    Sigma_2 = heritability_prop[3] * get(paste0("generate_", Sigma, "_Sigma"))(1000)
  }

  true = list(Sigma_0, Sigma_1, Sigma_2)

  rs = intersect(c(1, 3, 5), 1:(q-1))
  max_principal_angle = sapply(rs, function(r) mapply(max_principal_angle, estimate$Sigma_hat, true, r = r))
  rownames(max_principal_angle) = paste0("Sigma_", 1:length(D_list) - 1)
  colnames(max_principal_angle) = rs
  max_principal_angle = reshape2::melt(max_principal_angle, varnames = c("estimate", "r"))

  return(list(
    time = time,
    estimate = estimate,
    true = true,
    truncated = truncated,
    min_eigenvalue = min_eigenvalue,
    spectral_error = mapply(spectral_error, estimate$Sigma_hat, true),
    squared_error = mapply(squared_error, estimate$Sigma_hat, true),
    diag_squared_error = mapply(diag_squared_error, estimate$Sigma_hat, true),
    max_principal_angle = max_principal_angle,
    h2_error = (h2_prop(estimate$Sigma_hat, 2) - h2_prop(true, 2))^2
  ))

}

SIMULATION_ID = 1 # as.numeric(commandArgs(trailingOnly=TRUE)[1])

RESULT_PATH = paste0("simulation_hcp_results_", SIMULATION_ID)
dir.create(RESULT_PATH, recursive = TRUE)

replicates = 1:50

if (SIMULATION_ID == 1) { # 4800
  methods = c("mvHE", "mvREHE", "HE", "REHE", "REML")
  Sigmas = "uniform"
  ns = c(250, 500, 1000, 2000, 4000, 8000)
  qs = c(20, 100)
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n")
  qs = 5
  grid = rbind(grid, expand.grid(method = c(methods, "mvREML"), replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n"))
} else if (SIMULATION_ID == 2) { # 3000
  methods = c("mvHE", "mvREHE", "mvREHE_cvDR", paste0("DR", c(5), "-mvREML"))
  Sigmas = c("fast", "moderate", "slow")
  ns = c(500, 1000, 2000, 4000, 8000)
  qs = 1000
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n")
} else if (SIMULATION_ID == 3) { # 2400
  methods = c("mvHE", "mvREHE", "mvHE-smoothed", "mvREHE-smoothed")
  Sigmas = c("smooth_1", "smooth_2")
  qs = 100
  ns = c(125, 250, 500, 1000, 2000, 4000)
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n")
}

PARAMETER_ID = 1 # as.numeric(commandArgs(trailingOnly=TRUE)[2])
print(grid[PARAMETER_ID, ])
replicate = grid[PARAMETER_ID, "replicate"]
n = grid[PARAMETER_ID, "n"]
q = grid[PARAMETER_ID, "q"]
Sigma = grid[PARAMETER_ID, "Sigma"]
method = as.character(grid[PARAMETER_ID, "method"])
experiment = grid[PARAMETER_ID, "experiment"]

output = simulation(n, q, Sigma, method, replicate)
estimate = paste0("Sigma_", 1:3 - 1)

diag_squared_error = data.frame(replicate = replicate, estimate = estimate, diag_squared_error = output$diag_squared_error, n = n, q = q, Sigma = Sigma, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
squared_error = data.frame(replicate = replicate, estimate = estimate, squared_error = output$squared_error, n = n, q = q, Sigma = Sigma, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
spectral_error = data.frame(replicate = replicate, estimate = estimate, spectral_error = output$spectral_error, n = n, q = q, Sigma = Sigma, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
h2_error = data.frame(replicate = replicate, h2_error = output$h2_error, n = n, q = q, Sigma = Sigma, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
max_principal_angle = cbind(output$max_principal_angle, replicate = replicate, n = n, q = q, Sigma = Sigma, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
time = data.frame(replicate = replicate, method = method, time = output$time, n = n, q = q, Sigma = Sigma, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
truncated = data.frame(estimate = estimate, truncated = output$truncated, n = n, q = q, Sigma = Sigma, method = method, replicate = replicate, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
min_eigenvalue = data.frame(estimate = estimate, min_eigenvalue = output$min_eigenvalue, n = n, q = q, Sigma = Sigma, method = method, replicate = replicate, experiment = experiment, SIMULATION_ID = SIMULATION_ID)

saveRDS(list(diag_squared_error = diag_squared_error, squared_error = squared_error, spectral_error = spectral_error, h2_error = h2_error, max_principal_angle = max_principal_angle, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue), file.path(RESULT_PATH, paste0("n", n, "_q", q, "_Sigma", Sigma, "_replicate", replicate, "_experiment", experiment, "_method", method, ".rds")))
print(list(diag_squared_error = diag_squared_error, squared_error = squared_error, spectral_error = spectral_error, h2_error = h2_error, max_principal_angle = max_principal_angle, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue))
