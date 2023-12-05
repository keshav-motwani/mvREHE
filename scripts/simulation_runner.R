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

make_low_rank = function(A, r) {

  eig = eigen(A)
  eig$val[setdiff(1:nrow(A), 1:r)] = 0
  A = eig$vec %*% diag(eig$val) %*% t(eig$vec)
  attr(A, "sqrt") = eig$vec %*% diag(sqrt(eig$val)) %*% t(eig$vec)

  A

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
  matrix = V %*% diag(1/(1:q)^1) %*% t(V)
  attr(matrix, "sqrt") = V %*% diag(1/(1:q)^0.5) %*% t(V)
  matrix
}

generate_slow_Sigma = function(q) {
  V = pracma::randortho(q)
  matrix = V %*% diag(1/(1:q)^0.5) %*% t(V)
  attr(matrix, "sqrt") = V %*% diag(1/(1:q)^0.25) %*% t(V)
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

hcp_kinship = function(n) {

  kinship = R.matlab::readMat('data/kinship.mat'); # This is 2*K in the Solar-Eclipse notation
  K_G = as(kinship$K[[1]], "TsparseMatrix"); # Kinship matrix
  K_G = as.matrix(K_G)
  order = hclust(as.dist(-K_G))$order
  K_G = K_G[order, order]
  K_G = K_G[1:min(n, 1000), 1:min(n, 1000)]

  chol = chol(K_G[!duplicated(K_G), !duplicated(K_G)])
  expand = diag(1, nrow(K_G), nrow(K_G))
  expand = expand[, !duplicated(K_G)]
  expand[cbind(which(duplicated(K_G)), which(duplicated(K_G)) - 1:sum(duplicated(K_G)))] = 1
  chol = chol %*% t(expand)
  chol = as.matrix(Matrix::bdiag(replicate(ceiling(n / 1000), chol, simplify = FALSE)))

  K_G = as.matrix(Matrix::bdiag(replicate(ceiling(n / 1000), K_G, simplify = FALSE)))
  attr(K_G, "chol") = chol

  K_G

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

simulation = function(n, q, r, Sigma, method, replicate) {

  # D_1 = ar1_cor(n, 10, 0.5)
  D_1 = hcp_kinship(n)

  D_0 = diag(1, nrow = n, ncol = n)

  colnames(D_0) = colnames(D_1) = rownames(D_0) = rownames(D_1) = as.character(1:n)

  set.seed(123)
  Sigma_1 = get(paste0("generate_", Sigma, "_Sigma"))(q)
  Sigma_0 = get(paste0("generate_", Sigma, "_Sigma"))(q)

  chol_D_1 = attr(D_1, "chol")
  chol_D_0 = D_0
  sqrt_Sigma_1 = attr(Sigma_1, "sqrt")
  sqrt_Sigma_0 = attr(Sigma_0, "sqrt")

#  if (grepl("smooth", Sigma)) {
#    Sigma_0 = Sigma_0 + diag(1, q, q)
#    sqrt_Sigma_0 = sqrt_matrix(Sigma_0)
#  }

  set.seed(replicate)
  Gamma_1 = t(chol_D_1) %*% matrix(rnorm(nrow(chol_D_1) * q), nrow = nrow(chol_D_1)) %*% t(sqrt_Sigma_1)
  Epsilon = t(chol_D_0) %*% matrix(rnorm(n * q), nrow = n) %*% t(sqrt_Sigma_0)

  Y = Gamma_1 + Epsilon

  if (grepl("-", method) && substring(method, 1, 1) == "R") {
    PC = TRUE
    PC_Y = prcomp(Y, center = TRUE, scale. = FALSE)
    R = as.numeric(gsub("R", "", strsplit(method, "-")[[1]][1]))
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

  if (method == "mvHE") {
    time = system.time({estimate = mvHE(Y, list(D_0, D_1))})[3]
  } else if (method == "mvREHE") {
    time = system.time({estimate = mvREHE(Y, list(D_0, D_1))})[3]
  } else if (method == "mvREML") {
    source("scripts/mvREML.R")
    time = system.time({estimate = mvREML(Y, list(D_0, D_1))})[3]
  } else if (method == "HE") {
    time = system.time({estimate = univariate(Y, list(D_0, D_1), mvHE)})[3]
  } else if (method == "REHE") {
    time = system.time({estimate = univariate(Y, list(D_0, D_1), mvREHE)})[3]
  } else if (method == "REML") {
    source("scripts/mvREML.R")
    time = system.time({estimate = univariate(Y, list(D_0, D_1), mvREML)})[3]
  }

  if (PC) {
    estimate$Sigma_hat = lapply(estimate$Sigma_hat, function(x) PC_Y$rotation[, 1:R] %*% x %*% t(PC_Y$rotation[, 1:R]))
  }

  if (method == "mvHE") {
    truncated = 1 * c(attr(estimate$Sigma_hat[[1]], "truncated"), attr(estimate$Sigma_hat[[2]], "truncated"))
  } else {
    truncated = c(0, 0)
  }

  if (method == "mvHE") {
    min_eigenvalue = c(attr(estimate$Sigma_hat[[1]], "min_eigenvalue"), attr(estimate$Sigma_hat[[2]], "min_eigenvalue"))
  } else {
    min_eigenvalue = c(0, 0)
  }

  if (smoothed) {
    time = system.time({estimate$Sigma_hat = lapply(estimate$Sigma_hat, smooth_cov)})[3] + time
  }

  if (grepl("smooth", Sigma)) {
    Sigma_1 = get(paste0("generate_", Sigma, "_Sigma"))(1000)
    Sigma_0 = get(paste0("generate_", Sigma, "_Sigma"))(1000)
  }

  return(list(
    time = time,
    estimate = estimate,
    true = list(Sigma_0, Sigma_1),
    truncated = truncated,
    min_eigenvalue = min_eigenvalue,
    spectral_error = mapply(spectral_error, estimate$Sigma_hat, list(Sigma_0, Sigma_1)),
    squared_error = mapply(squared_error, estimate$Sigma_hat, list(Sigma_0, Sigma_1)),
    diag_squared_error = mapply(diag_squared_error, estimate$Sigma_hat, list(Sigma_0, Sigma_1))
  ))

}

SIMULATION_ID = as.numeric(commandArgs(trailingOnly=TRUE)[1])

RESULT_PATH = paste0("simulation_hcp_results_", SIMULATION_ID)
dir.create(RESULT_PATH, recursive = TRUE)

replicates = 1:50
rs = c(Inf)

if (SIMULATION_ID == 1) {
  methods = c("mvHE", "mvREHE", "HE", "REHE", "REML")
  Sigmas = "uniform"
  ns = c(500, 1000, 2000, 4000, 8000, 16000)
  qs = c(20, 100)
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, r = rs, Sigma = Sigmas, experiment = "n")
  qs = 5
  grid = rbind(grid, expand.grid(method = c(methods, "mvREML"), replicate = replicates, n = ns, q = qs, r = rs, Sigma = Sigmas, experiment = "n"))
} else if (SIMULATION_ID == 2) {
  methods = c(paste0("R", c(5, 10, 50, 100), "-mvREHE"), "mvREHE", paste0("R", c(5), "-mvREML"))
  Sigmas = c("constant", "slow", "fast")
  ns = c(500, 1000, 2000, 4000, 8000)
  qs = 1000
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, r = rs, Sigma = Sigmas, experiment = "n")
} else if (SIMULATION_ID == 3) {
  methods = c("mvHE", "mvREHE")
  Sigmas = c("smooth_1", "smooth_2")
  qs = c(25, 50, 100, 200, 500, 1000)
  ns = 500
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, r = rs, Sigma = Sigmas, experiment = "q")
  qs = 100
  ns = c(125, 250, 500, 1000, 2000, 4000)
  grid = rbind(grid, expand.grid(method = methods, replicate = replicates, n = ns, q = qs, r = rs, Sigma = Sigmas, experiment = "n"))
}

PARAMETER_ID = as.numeric(commandArgs(trailingOnly=TRUE)[2])
print(grid[PARAMETER_ID, ])
replicate = grid[PARAMETER_ID, "replicate"]
n = grid[PARAMETER_ID, "n"]
q = grid[PARAMETER_ID, "q"]
r = grid[PARAMETER_ID, "r"]
Sigma = grid[PARAMETER_ID, "Sigma"]
method = as.character(grid[PARAMETER_ID, "method"])
experiment = grid[PARAMETER_ID, "experiment"]

output = simulation(n, q, r, Sigma, method, replicate)

diag_squared_error = data.frame(replicate = replicate, estimate = c("Sigma_0", "Sigma_1"), diag_squared_error = output$diag_squared_error, n = n, q = q, Sigma = Sigma, r = r, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
squared_error = data.frame(replicate = replicate, estimate = c("Sigma_0", "Sigma_1"), squared_error = output$squared_error, n = n, q = q, Sigma = Sigma, r = r, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
spectral_error = data.frame(replicate = replicate, estimate = c("Sigma_0", "Sigma_1"), spectral_error = output$spectral_error, n = n, q = q, Sigma = Sigma, r = r, method = method, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
time = data.frame(replicate = replicate, method = method, time = output$time, n = n, q = q, Sigma = Sigma, r = r, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
truncated = data.frame(estimate = c("Sigma_0", "Sigma_1"), truncated = output$truncated, n = n, q = q, Sigma = Sigma, r = r, method = method, replicate = replicate, experiment = experiment, SIMULATION_ID = SIMULATION_ID)
min_eigenvalue = data.frame(estimate = c("Sigma_0", "Sigma_1"), min_eigenvalue = output$min_eigenvalue, n = n, q = q, Sigma = Sigma, r = r, method = method, replicate = replicate, experiment = experiment, SIMULATION_ID = SIMULATION_ID)

saveRDS(list(diag_squared_error = diag_squared_error, squared_error = squared_error, spectral_error = spectral_error, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue), file.path(RESULT_PATH, paste0("n", n, "_q", q, "_Sigma", Sigma, "_r", r, "_replicate", replicate, "_experiment", experiment, "_method", method, ".rds")))
print(list(diag_squared_error = diag_squared_error, squared_error = squared_error, spectral_error = spectral_error, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue))
