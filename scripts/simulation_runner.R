library(mvREHE)

RESULT_PATH = "results/"
dir.create(RESULT_PATH, recursive = TRUE)

squared_error = function(A, B) {
  if (!is.null(A) & !is.null(B)) {
    sqrt(sum((A - B) ^ 2))
  } else {
    NA
  }
}

spectral_error = function(A, B) {
  if (!is.null(A) & !is.null(B)) {
    norm(A - B, "2")
  } else {
    NA
  }
}

EV1_error = function(A, B) {
  if (!is.null(A) & !is.null(B)) {
    sqrt(sum((eigen(A)$vec[, 1] - eigen(B)$vec[, 1])^2))
  } else {
    NA
  }
}

make_low_rank = function(A, r) {

  eig = eigen(A)
  eig$val[(r + 1):length(eig$val)] = 0
  eig$vec %*% diag(eig$val) %*% t(eig$vec)

}

simulation = function(n, q) {

  r = ifelse(q > 10, round(q / 2), q)

  D_1 = matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      D_1[i, j] = 0.5^abs(i - j)
    }
  }

  D_0 = diag(1, nrow = n, ncol = n)

  colnames(D_0) = colnames(D_1) = rownames(D_0) = rownames(D_1) = as.character(1:n)

  Sigma_1 = make_low_rank(clusterGeneration::rcorrmatrix(q), r)
  Sigma_0 = make_low_rank(clusterGeneration::rcorrmatrix(q), r)

  chol_D_1 = chol(D_1)
  chol_D_0 = chol(D_0)
  e1 = eigen(Sigma_1)
  sqrt_Sigma_1 = e1$vec %*% diag(sqrt(e1$val)) %*% t(e1$vec)
  e0 = eigen(Sigma_0)
  sqrt_Sigma_0 = e0$vec %*% diag(sqrt(e0$val)) %*% t(e0$vec)

  Gamma_1 = t(chol_D_1) %*% matrix(rnorm(n * q), nrow = n) %*% t(sqrt_Sigma_1)
  Epsilon = t(chol_D_0) %*% matrix(rnorm(n * q), nrow = n) %*% t(sqrt_Sigma_0)

  Y = Gamma_1 + Epsilon

  Sigma_init_list = lapply(1:2, function(i) clusterGeneration::rcorrmatrix(q))
  L_init_list = lapply(Sigma_init_list, function(i) t(chol(i)))
  mvHE_time = system.time({mvHE_estimate = mvREHE::mvHE(Y, list(D_0, D_1))})[3]
  mvREHE_time = system.time({mvREHE_estimate = mvREHE::mvREHE(Y, list(D_0, D_1), lambda = 0, L_init_list = "mvHE", algorithm = "L-BFGS-B")})[3]
  cv_mvREHE_time = system.time({cv_mvREHE_estimate = mvREHE::cv_mvREHE(Y, list(D_0, D_1), K = 5, L_init_list = "mvHE", algorithm = "L-BFGS-B")})[3]
  mvLRHE_time = system.time({mvLRHE_estimate = mvLRHE(Y, list(D_0, D_1), r = q, Sigma_init_list = "mvHE")})[3]
  cv_mvLRHE_time = system.time({cv_mvLRHE_estimate = cv_mvLRHE(Y, list(D_0, D_1), K = 5, Sigma_init_list = "mvHE")})[3]
  orc_mvLRHE_time = system.time({orc_mvLRHE_estimate = mvLRHE(Y, list(D_0, D_1), r = r, Sigma_init_list = "mvHE")})[3]
  if (FALSE) {# q <= 5 & n <= 400) {
    mvREML_time = system.time({mvREML_estimate = mvREML(Y, D_0, D_1)})[3]
  } else {
    mvREML_time = NA
    mvREML_estimate = NULL
  }
  naive = list(Sigma_hat = list(clusterGeneration::rcorrmatrix(q), clusterGeneration::rcorrmatrix(q)))

  return(list(time = c(mvHE = mvHE_time, mvREHE = mvREHE_time, cv_mvREHE = cv_mvREHE_time, mvLRHE = mvLRHE_time, cv_mvLRHE = cv_mvLRHE_time, orc_mvLRHE = orc_mvLRHE_time, mvREML = mvREML_time),
              truncated = 1 * c(attr(mvHE_estimate$Sigma_hat[[1]], "truncated"), attr(mvHE_estimate$Sigma_hat[[2]], "truncated")),
              min_eigenvalue = c(attr(mvHE_estimate$Sigma_hat[[1]], "min_eigenvalue"), attr(mvHE_estimate$Sigma_hat[[2]], "min_eigenvalue")),
              squared_error = lapply(list(mvHE = mvHE_estimate, mvREHE = mvREHE_estimate, cv_mvREHE = cv_mvREHE_estimate, mvLRHE = mvLRHE_estimate, cv_mvLRHE = cv_mvLRHE_estimate, orc_mvLRHE = orc_mvLRHE_estimate, mvREML = mvREML_estimate, naive = naive), function(estimate) c(squared_error(Sigma_0, estimate$Sigma_hat[[1]]), squared_error(Sigma_1, estimate$Sigma_hat[[2]]))),
              spectral_error = lapply(list(mvHE = mvHE_estimate, mvREHE = mvREHE_estimate, cv_mvREHE = cv_mvREHE_estimate, mvLRHE = mvLRHE_estimate, cv_mvLRHE = cv_mvLRHE_estimate, orc_mvLRHE = orc_mvLRHE_estimate, mvREML = mvREML_estimate, naive = naive), function(estimate) c(spectral_error(Sigma_0, estimate$Sigma_hat[[1]]), spectral_error(Sigma_1, estimate$Sigma_hat[[2]]))),
              EV1_error = lapply(list(mvHE = mvHE_estimate, mvREHE = mvREHE_estimate, cv_mvREHE = cv_mvREHE_estimate, mvLRHE = mvLRHE_estimate, cv_mvLRHE = cv_mvLRHE_estimate, orc_mvLRHE = orc_mvLRHE_estimate, mvREML = mvREML_estimate, naive = naive), function(estimate) c(EV1_error(Sigma_0, estimate$Sigma_hat[[1]]), EV1_error(Sigma_1, estimate$Sigma_hat[[2]])))))

}

replicates = 1:50
ns = 200 * 1:5
qs = c(3, 10)
grid = expand.grid(replicate = replicates, n = ns, q = qs, experiment = "n")
ns = c(200)
qs = c(3, 5, 10, 25, 50)
grid = rbind(grid, expand.grid(replicate = replicates, n = ns, q = qs, experiment = "q"))
ns = c(600)
qs = c(3, 5, 10)
grid = rbind(grid, expand.grid(replicate = replicates, n = ns, q = qs, experiment = "q"))
ns = c(600)
qs = c(25, 50)
grid = rbind(grid, expand.grid(replicate = replicates, n = ns, q = qs, experiment = "q"))

PARAMETER_ID = as.numeric(commandArgs(trailingOnly=TRUE)[1])
replicate = grid[PARAMETER_ID, "replicate"]
n = grid[PARAMETER_ID, "n"]
q = grid[PARAMETER_ID, "q"]
experiment = grid[PARAMETER_ID, "experiment"]

set.seed(replicate, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
output = simulation(n, q)

spectral_error = list()
for (method in names(output$spectral_error)) {
  spectral_error = c(spectral_error, list(data.frame(replicate = replicate, estimate = c("Sigma_0", "Sigma_1"), spectral_error = output$spectral_error[[method]], n = n, q = q, method = method, experiment = experiment)))
}
spectral_error = do.call(rbind, spectral_error)

EV1_error = list()
for (method in names(output$EV1_error)) {
  EV1_error = c(EV1_error, list(data.frame(replicate = replicate, estimate = c("Sigma_0", "Sigma_1"), EV1_error = output$EV1_error[[method]], n = n, q = q, method = method, experiment = experiment)))
}
EV1_error = do.call(rbind, EV1_error)

squared_error = list()
for (method in names(output$squared_error)) {
  squared_error = c(squared_error, list(data.frame(replicate = replicate, estimate = c("Sigma_0", "Sigma_1"), squared_error = output$squared_error[[method]], n = n, q = q, method = method, experiment = experiment)))
}
squared_error = do.call(rbind, squared_error)

time = data.frame(replicate = replicate, method = names(output$time), time = output$time, n = n, q = q, experiment = experiment)

truncated = data.frame(estimate = c("Sigma_0", "Sigma_1"), truncated = output$truncated, n = n, q = q, method = "mvHE", replicate = replicate, experiment = experiment)

min_eigenvalue = data.frame(estimate = c("Sigma_0", "Sigma_1"), min_eigenvalue = output$min_eigenvalue, n = n, q = q, method = "mvHE", replicate = replicate, experiment = experiment)

saveRDS(list(EV1_error = EV1_error, spectral_error = spectral_error, squared_error = squared_error, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue), file.path(RESULT_PATH, paste0("n", n, "_q", q, "_replicate", replicate, "_experiment", experiment, ".rds")))
print(list(EV1_error = EV1_error, spectral_error = spectral_error, squared_error = squared_error, time = time, truncated = truncated, min_eigenvalue = min_eigenvalue))
