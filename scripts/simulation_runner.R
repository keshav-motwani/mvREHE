source("scripts/simulation_setup.R")

SIMULATION_ID = commandArgs(trailingOnly=TRUE)[1]
COMPONENTS = commandArgs(trailingOnly=TRUE)[2]

RESULT_PATH = paste0("simulation_", COMPONENTS, "_components_", SIMULATION_ID)
DATA_ANALYSIS_RESULT_PATH = paste0("data_analysis_", COMPONENTS, "_components")
dir.create(RESULT_PATH, recursive = TRUE)

replicates = 1:50

if (grepl("lowdim", SIMULATION_ID)) { # 5000
  methods = c("mvHE", "mvREHE", "mvREML", "HE", "REHE", "REML")
  Sigmas = "uniform"
  ns = c(250, 500, 1000, 2000, 4000, 8000)
  qs = c(5, 10, 20)
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n")
  grid = grid[grid$method != "mvREML" | (grid$method == "mvREML" & grid$q == 5) | (grid$method == "mvREML" & grid$q == 10 & grid$n <= 2000), ]
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
} else if (SIMULATION_ID == "data") { # 1250
  methods = c("mvHE", "mvREHE", "HE", "REHE", "REML")
  Sigmas = "data_100"
  ns = c(500, 1000, 2000, 4000, 8000)
  qs = NA
  grid = expand.grid(method = methods, replicate = replicates, n = ns, q = qs, Sigma = Sigmas, experiment = "n")
}

PARAMETER_ID = as.numeric(commandArgs(trailingOnly=TRUE)[3])
print(grid[PARAMETER_ID, ])
replicate = grid[PARAMETER_ID, "replicate"]
n = grid[PARAMETER_ID, "n"]
q = grid[PARAMETER_ID, "q"]
Sigma = grid[PARAMETER_ID, "Sigma"]
method = as.character(grid[PARAMETER_ID, "method"])
experiment = grid[PARAMETER_ID, "experiment"]

output = simulation(COMPONENTS, n, q, Sigma, method, SIMULATION_ID, replicate, DATA_ANALYSIS_RESULT_PATH)
estimate = paste0("Sigma_", 1:COMPONENTS - 1)

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
