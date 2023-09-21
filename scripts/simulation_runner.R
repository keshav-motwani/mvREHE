library(mvREHE)
library(Matrix)

spectral_error = function(A, B) {
  if (!is.null(A) & !is.null(B)) {
    norm(A - B, "2")
  } else {
    NA
  }
}

squared_error = function(A, B) {
  if (!is.null(A) & !is.null(B)) {
    norm(A - B, "F")
  } else {
    NA
  }
}

diag_squared_error = function(A, B) {
  if (!is.null(A) & !is.null(B)) {
    sqrt(sum(diag(A - B)^2))
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

simulate_Sigma = function(q) {
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

sqrt_matrix = function(A) {

  eig = eigen(A)
  eig$vec %*% diag(sqrt(eig$val)) %*% t(eig$vec)

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

simulation = function(n, q, r, method) {

  # D_1 = ar1_cor(n, 10, 0.5)
  D_1 = hcp_kinship(n)

  D_0 = diag(1, nrow = n, ncol = n)

  colnames(D_0) = colnames(D_1) = rownames(D_0) = rownames(D_1) = as.character(1:n)

  Sigma_1 = simulate_Sigma(q)
  Sigma_0 = simulate_Sigma(q)

  chol_D_1 = attr(D_1, "chol")
  chol_D_0 = D_0
  sqrt_Sigma_1 = attr(Sigma_1, "sqrt")
  sqrt_Sigma_0 = attr(Sigma_0, "sqrt")

  Gamma_1 = t(chol_D_1) %*% matrix(rnorm(nrow(chol_D_1) * q), nrow = nrow(chol_D_1)) %*% t(sqrt_Sigma_1)
  Epsilon = t(chol_D_0) %*% matrix(rnorm(n * q), nrow = n) %*% t(sqrt_Sigma_0)

  Y = Gamma_1 + Epsilon

  if (grepl("-", method)) {
    PC = TRUE
    PC_Y = prcomp(Y, center = TRUE, scale. = FALSE)
    R = as.numeric(gsub("R", "", strsplit(method, "-")[[1]][1]))
    Y = PC_Y$x[, 1:R]
    method = strsplit(method, "-")[[1]][2]
  } else {
    PC = FALSE
  }

  if (method == "mvHE") {
    time = system.time({estimate = mvHE(Y, list(D_0, D_1))})[3]
  } else if (method == "mvREHE") {
    time = system.time({estimate = mvREHE(Y, list(D_0, D_1))})[3]
  } else if (method == "cv_mvREHE_L2") {
    time = system.time({estimate = cv_mvREHE_L2(Y, list(D_0, D_1))})[3]
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
    truncated = truncated = 1 * c(attr(estimate$Sigma_hat[[1]], "truncated"), attr(estimate$Sigma_hat[[2]], "truncated"))
  } else {
    truncated = c(0, 0)
  }

  if (method == "mvHE") {
    min_eigenvalue = c(attr(estimate$Sigma_hat[[1]], "min_eigenvalue"), attr(estimate$Sigma_hat[[2]], "min_eigenvalue"))
  } else {
    min_eigenvalue = c(0, 0)
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

ns = c(500, 1000, 2000, 4000, 8000)
replicates = 1:50
rs = c(Inf)

if (SIMULATION_ID == 1) {
  methods = c("mvHE", "mvREHE", "HE", "REHE", "REML")
  qs = c(20, 100)
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, r = rs, experiment = "n")
  qs = 5
  grid = rbind(grid, expand.grid(method = c(methods, "mvREML"), replicate = replicates, n = ns, q = qs, r = rs, experiment = "n"))
} else if (SIMULATION_ID == 2) {
  methods = c(paste0("R", c(5, 10, 50, 100), "-mvREHE"), "mvREHE", paste0("R", c(5, 10), "-mvREML"))
  qs = c(500, 1000, 2000)
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, r = rs, experiment = "n")
}

PARAMETER_ID = as.numeric(commandArgs(trailingOnly=TRUE)[2])
print(grid[PARAMETER_ID, ])
replicate = grid[PARAMETER_ID, "replicate"]
n = grid[PARAMETER_ID, "n"]
q = grid[PARAMETER_ID, "q"]
r = grid[PARAMETER_ID, "r"]
method = grid[PARAMETER_ID, "method"]
experiment = grid[PARAMETER_ID, "experiment"]

set.seed(replicate, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
output = simulation(n, q, r, method)

diag_squared_error = data.frame(replicate = replicate, estimate = c("Sigma_0", "Sigma_1"), diag_squared_error = output$diag_squared_error, n = n, q = q, r = r, method = method, experiment = experiment)
squared_error = data.frame(replicate = replicate, estimate = c("Sigma_0", "Sigma_1"), squared_error = output$squared_error, n = n, q = q, r = r, method = method, experiment = experiment)
spectral_error = data.frame(replicate = replicate, estimate = c("Sigma_0", "Sigma_1"), spectral_error = output$spectral_error, n = n, q = q, r = r, method = method, experiment = experiment)
time = data.frame(replicate = replicate, method = method, time = output$time, n = n, q = q, r = r, experiment = experiment)
truncated = data.frame(estimate = c("Sigma_0", "Sigma_1"), truncated = output$truncated, n = n, q = q, r = r, method = method, replicate = replicate, experiment = experiment)
min_eigenvalue = data.frame(estimate = c("Sigma_0", "Sigma_1"), min_eigenvalue = output$min_eigenvalue, n = n, q = q, r = r, method = method, replicate = replicate, experiment = experiment)

saveRDS(list(output = output, diag_squared_error = diag_squared_error, squared_error = squared_error, spectral_error = spectral_error, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue), file.path(RESULT_PATH, paste0("n", n, "_q", q, "_r", r, "_replicate", replicate, "_experiment", experiment, "_method", method, ".rds")))
print(list(diag_squared_error = diag_squared_error, squared_error = squared_error, spectral_error = spectral_error, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue))
