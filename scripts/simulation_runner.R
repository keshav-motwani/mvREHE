library(mvREHE)
source("scripts/mvREML.R")

RESULT_PATH = "results_new_rank_10/"
dir.create(RESULT_PATH, recursive = TRUE)

spectral_error = function(A, B) {
  if (!is.null(A) & !is.null(B)) {
    norm(A - B, "2")
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

sqrt_matrix = function(A) {

  eig = eigen(A)
  eig$vec %*% diag(sqrt(eig$val)) %*% t(eig$vec)

}

ar1_cor = function(n, m, rho) {

  exponent = abs(matrix(1:m - 1, nrow = m, ncol = m, byrow = TRUE) - (1:m - 1))
  mat = rho^exponent

  as.matrix(Matrix::bdiag(replicate(n / m, mat, simplify = FALSE)))

}

simulation = function(n, q, method) {

  D_1 = ar1_cor(n, 10, 0.5)

  D_0 = diag(1, nrow = n, ncol = n)

  colnames(D_0) = colnames(D_1) = rownames(D_0) = rownames(D_1) = as.character(1:n)

  Sigma_1 = make_low_rank(clusterGeneration::rcorrmatrix(q), 10)
  Sigma_0 = make_low_rank(clusterGeneration::rcorrmatrix(q), q)

  chol_D_1 = chol(D_1)
  chol_D_0 = D_0
  sqrt_Sigma_1 = attr(Sigma_1, "sqrt")
  sqrt_Sigma_0 = attr(Sigma_0, "sqrt")

  Gamma_1 = t(chol_D_1) %*% matrix(rnorm(n * q), nrow = n) %*% t(sqrt_Sigma_1)
  Epsilon = t(chol_D_0) %*% matrix(rnorm(n * q), nrow = n) %*% t(sqrt_Sigma_0)

  Y = Gamma_1 + Epsilon

  Sigma_init_list = lapply(1:2, function(i) clusterGeneration::rcorrmatrix(q))
  L_init_list = lapply(Sigma_init_list, function(i) t(chol(i)))

  if (method == "mvHE") {
    time = system.time({estimate = mvHE(Y, list(D_0, D_1))})[3]
  } else if (method == "mvREHE") {
    time = system.time({estimate = mvREHE(Y, list(D_0, D_1), Sigma_init_list = Sigma_init_list)})[3]
  } else if (method == "mvREHE_L2") {
    time = system.time({estimate = cv_mvREHE_L2(Y, list(D_0, D_1), K = 5, Sigma_init_list = Sigma_init_list)})[3]
  } else if (method == "mvREHE_rank") {
    time = system.time({estimate = cv_mvREHE_rank(Y, list(D_0, D_1), K = 5, Sigma_init_list = Sigma_init_list)})[3]
  } else if (method == "mvREML") {
    time = system.time({estimate = mvREML(Y, D_0, D_1)})[3]
  } else if (method == "naive") {
    time = NA
    estimate = list(Sigma_hat = list(clusterGeneration::rcorrmatrix(q), clusterGeneration::rcorrmatrix(q)))
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
    spectral_error = mapply(spectral_error, estimate$Sigma_hat, list(Sigma_0, Sigma_1))
  ))

}

methods = c("mvHE", "mvREHE", "mvREHE_L2", "mvREHE_rank", "mvREML", "naive")
replicates = 1:50
ns = 200 * 1:5
qs = c(3, 50)
grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, experiment = "n")
ns = c(200, 600)
qs = c(3, 5, 10, 25, 50, 100, 250)
grid = rbind(grid, expand.grid(method = methods, replicate = replicates, n = ns, q = qs, experiment = "q"))

PARAMETER_ID = as.numeric(commandArgs(trailingOnly=TRUE)[1])
replicate = grid[PARAMETER_ID, "replicate"]
n = grid[PARAMETER_ID, "n"]
q = grid[PARAMETER_ID, "q"]
method = grid[PARAMETER_ID, "method"]
experiment = grid[PARAMETER_ID, "experiment"]

set.seed(replicate, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
output = simulation(n, q, method)

spectral_error = data.frame(replicate = replicate, estimate = c("Sigma_0", "Sigma_1"), spectral_error = output$spectral_error, n = n, q = q, method = method, experiment = experiment)
time = data.frame(replicate = replicate, method = method, time = output$time, n = n, q = q, experiment = experiment)
truncated = data.frame(estimate = c("Sigma_0", "Sigma_1"), truncated = output$truncated, n = n, q = q, method = method, replicate = replicate, experiment = experiment)
min_eigenvalue = data.frame(estimate = c("Sigma_0", "Sigma_1"), min_eigenvalue = output$min_eigenvalue, n = n, q = q, method = method, replicate = replicate, experiment = experiment)

saveRDS(list(spectral_error = spectral_error, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue), file.path(RESULT_PATH, paste0("n", n, "_q", q, "_replicate", replicate, "_experiment", experiment, "_method", method, ".rds")))
print(list(spectral_error = spectral_error, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue))
