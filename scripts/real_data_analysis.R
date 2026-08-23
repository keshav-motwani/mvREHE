library(tidyverse)
library(mvREHE)
library(Matrix)

source("scripts/latent_regression.R")

#################################
##########  Functions ###########
#################################

from_conn_to_vec = function(conn)
{
  as.vector(conn[lower.tri(conn, diag = TRUE)])
}

sqrt_matrix = function(A) {
  eig = eigen(A)
  eig$vec %*% diag(sqrt(pmax(eig$val, 0))) %*% t(eig$vec)
}

#################################
##########  Load data ###########
#################################

DATA_PATH = "data/"
RESULT_PATH = "real_data_analysis_ica_260804"
dir.create(RESULT_PATH, recursive = TRUE)

data = readRDS(file.path(DATA_PATH, "ica_clean_data.rds"))

fun_connectomes = data$fun_connectomes
str_connectomes = data$str_connectomes
K_G = as(data$K_G, "dsCMatrix")
X = data$X
groups = data$groups

Y_fun = sapply(fun_connectomes, from_conn_to_vec) %>% t
Y_str = sapply(str_connectomes, from_conn_to_vec) %>% t
colnames(Y_fun) = paste0("fun", 1:ncol(Y_fun))
colnames(Y_str) = paste0("str", 1:ncol(Y_str))

Y = cbind(Y_fun, Y_str)

fun_indices = grep("fun", colnames(Y))
str_indices = grep("str", colnames(Y))

connection_names = outer(groups, groups, "paste")
connection_names = connection_names[lower.tri(connection_names, diag = TRUE)]

fixed_effects = TRUE
residuals = lsfit(X, Y, intercept = FALSE)$residuals
colnames(residuals) = colnames(Y)
Y = residuals

scale_type = "column"
Y = scale(Y)
Y[, attr(Y, "scaled:scale") < 1e-14] = 0

D_list = list(as(as(Matrix::Diagonal(nrow(Y)), "dgCMatrix"), "dsCMatrix"), K_G, (K_G > 0) * 1)


print("getting family ids from kinship matrix")
fam = family_ids_from_kinship(K_G)

print("fitting mvREHE model")
if (!file.exists(file.path(RESULT_PATH, "fit.rds"))) {
  fit = mvREHE(Y, D_list)
  saveRDS(fit, file.path(RESULT_PATH, "fit.rds"))
  raw_cor = cor(Y)
  saveRDS(raw_cor, file.path(RESULT_PATH, "raw_cor.rds"))
} else {
  fit = readRDS(file.path(RESULT_PATH, "fit.rds"))
}

CORES = 180

ARRAY_ID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
JOB_ID   = Sys.getenv("SLURM_JOB_ID")
CHECKPOINT_PATH = if (nchar(JOB_ID) > 0) {
  file.path("/gscratch/stf/kmotwani/app_ckpts", paste0(JOB_ID, "_", ARRAY_ID))
} else {
  NULL
}

if (ARRAY_ID == 1) {

  estimator_fn = function(Y, D_list) {
    vars = matrixStats::colVars(Y)
    w = 1 / vars
    w[vars < 1e-8] = 1
    mvREHE(Y, D_list, return_full = FALSE, w_columns = w)
  }

  RANK_SEQ = c(5, 10, 20, 30)

  matrix_grid_fn = function(fit, outcomes, covariates) matrix_tuning_grid_fn(fit, outcomes, covariates, rank_seq = RANK_SEQ)
  print("running nested cross-validation for latent matrix regression")
  if (!file.exists(file.path(RESULT_PATH, "latent_matrix_regression_nested.rds"))) {
    set.seed(2026)
    latent_matrix_regression_nested = nested_cv_latent_regression(
      Y, D_list, fam, fit, fun_indices, str_indices,
      tuning_grid_fn = matrix_grid_fn, estimator_fn = estimator_fn,
      regression_fn = low_rank_regression_from_cov,
      K_outer = 5, K_inner = 5, cores = CORES,
      checkpoint_path = CHECKPOINT_PATH)
    saveRDS(latent_matrix_regression_nested,
            file.path(RESULT_PATH, "latent_matrix_regression_nested.rds"))
  }

}

if (ARRAY_ID == 2) {

  N_LAMBDA = 10

  estimator_fn = function(Y, D_list) {
    vars = matrixStats::colVars(Y)
    w = 1 / vars
    w[vars < 1e-8] = 1
    mvREHE(Y, D_list, return_full = TRUE, w_columns = w)
  }

  lasso_grid_fn = function(fit, outcomes, covariates) lasso_tuning_grid_fn(fit, outcomes, covariates, n_lambda = N_LAMBDA)
  print("running nested cross-validation for latent lasso regression")
  if (!file.exists(file.path(RESULT_PATH, "latent_lasso_regression_nested.rds"))) {
    set.seed(2026)
    latent_lasso_regression_nested = nested_cv_latent_regression(
      Y, D_list, fam, fit, fun_indices, str_indices,
      tuning_grid_fn = lasso_grid_fn, estimator_fn = estimator_fn,
      regression_fn = lasso_regression_from_cov,
      K_outer = 5, K_inner = 5, cores = CORES,
      checkpoint_path = CHECKPOINT_PATH)
    saveRDS(latent_lasso_regression_nested,
            file.path(RESULT_PATH, "latent_lasso_regression_nested.rds"))
  }

}

if (ARRAY_ID == 3) {

  N_LAMBDA = 10

  estimator_fn = function(Y, D_list) {
    vars = matrixStats::colVars(Y)
    w = 1 / vars
    w[vars < 1e-8] = 1
    mvREHE(Y, D_list, return_full = TRUE, w_columns = w)
  }

  ridge_grid_fn = function(fit, outcomes, covariates) ridge_tuning_grid_fn(fit, outcomes, covariates, n_lambda = N_LAMBDA)
  print("running nested cross-validation for latent ridge regression")
  if (!file.exists(file.path(RESULT_PATH, "latent_ridge_regression_nested.rds"))) {
    set.seed(2026)
    latent_ridge_regression_nested = nested_cv_latent_regression(
      Y, D_list, fam, fit, fun_indices, str_indices,
      tuning_grid_fn = ridge_grid_fn, estimator_fn = estimator_fn,
      regression_fn = ridge_regression_from_cov,
      K_outer = 5, K_inner = 5, cores = CORES,
      checkpoint_path = CHECKPOINT_PATH)
    saveRDS(latent_ridge_regression_nested,
            file.path(RESULT_PATH, "latent_ridge_regression_nested.rds"))
  }

}

if (ARRAY_ID == 4) {

  RANK_SEQ = c(5, 10, 20, 30)

  raw_matrix_grid_fn = function(Sigma, outcomes, covariates) {
    lapply(seq_along(outcomes), function(o) data.frame(rank = RANK_SEQ))
  }
  print("running nested cross-validation for raw matrix regression")
  if (!file.exists(file.path(RESULT_PATH, "raw_matrix_regression_nested.rds"))) {
    set.seed(2026)
    raw_matrix_regression_nested = nested_cv_raw_regression(
      Y, fam, fun_indices, str_indices,
      tuning_grid_fn = raw_matrix_grid_fn,
      regression_fn = low_rank_regression_from_cov,
      K_outer = 5, K_inner = 5, cores = CORES, use_svd = TRUE)
    saveRDS(raw_matrix_regression_nested,
            file.path(RESULT_PATH, "raw_matrix_regression_nested.rds"))
  }

}

if (ARRAY_ID == 5) {

  N_LAMBDA = 10

  raw_lasso_grid_fn = function(Sigma, outcomes, covariates) {
    vars = c(covariates, outcomes)
    weights = weights_cov_to_cor(Sigma[vars, vars])
    s_cov = weights[seq_along(covariates)]
    s_out = weights[seq_along(outcomes) + length(covariates)]
    lapply(seq_along(outcomes), function(o) {
      outcome = outcomes[o]
      lambda_max = if (Sigma[outcome, outcome] > 1e-10)
        max(abs(Sigma[covariates, outcome] * s_cov * s_out[o]))
      else 0
      if (lambda_max > 0)
        data.frame(lambda = exp(seq(log(lambda_max), log(lambda_max * 1e-4),
                                   length.out = N_LAMBDA)))
      else
        data.frame(lambda = rep(NA_real_, N_LAMBDA))
    })
  }
  print("running nested cross-validation for raw lasso regression")
  if (!file.exists(file.path(RESULT_PATH, "raw_lasso_regression_nested.rds"))) {
    set.seed(2026)
    raw_lasso_regression_nested = nested_cv_raw_regression(
      Y, fam, fun_indices, str_indices,
      tuning_grid_fn = raw_lasso_grid_fn,
      regression_fn = lasso_regression_from_cov,
      K_outer = 5, K_inner = 5, cores = CORES)
    saveRDS(raw_lasso_regression_nested,
            file.path(RESULT_PATH, "raw_lasso_regression_nested.rds"))
  }

}

if (ARRAY_ID == 6) {

  N_LAMBDA = 10

  raw_ridge_grid_fn = function(Sigma, outcomes, covariates) raw_tuning_grid_fn(Sigma, outcomes, covariates, n_lambda = N_LAMBDA)
  print("running nested cross-validation for raw ridge regression")
  if (!file.exists(file.path(RESULT_PATH, "raw_ridge_regression_nested.rds"))) {
    set.seed(2026)
    raw_ridge_regression_nested = nested_cv_raw_regression(
      Y, fam, fun_indices, str_indices,
      tuning_grid_fn = raw_ridge_grid_fn,
      regression_fn = ridge_regression_from_cov,
      K_outer = 5, K_inner = 5, cores = CORES)
    saveRDS(raw_ridge_regression_nested,
            file.path(RESULT_PATH, "raw_ridge_regression_nested.rds"))
  }

}

if (ARRAY_ID == 7) {

  estimator_fn = function(Y, D_list) {
    vars = matrixStats::colVars(Y)
    w = 1 / vars
    w[vars < 1e-8] = 1
    mvREHE(Y, D_list, return_full = FALSE, w_columns = w)
  }
  RANK_SEQ = c(5, 10, 20, 30)
  matrix_grid_fn = function(fit, outcomes, covariates) matrix_tuning_grid_fn(fit, outcomes, covariates, rank_seq = RANK_SEQ)

  print("running cv for latent matrix regression tuning selection")
  if (!file.exists(file.path(RESULT_PATH, "latent_matrix_cv.rds"))) {
    set.seed(2026)
    latent_matrix_cv = cv_latent_regression(
      Y, D_list, fam, fit, fun_indices, str_indices,
      tuning_grid_fn = matrix_grid_fn, estimator_fn = estimator_fn,
      regression_fn = low_rank_regression_from_cov,
      K = 5, cores = CORES)
    saveRDS(latent_matrix_cv, file.path(RESULT_PATH, "latent_matrix_cv.rds"))
  } else {
    latent_matrix_cv = readRDS(file.path(RESULT_PATH, "latent_matrix_cv.rds"))
  }

  p_ltr = length(str_indices)
  p_str = round((sqrt(8 * p_ltr + 1) - 1) / 2)
  m = length(fun_indices)

  if (!file.exists(file.path(RESULT_PATH, "latent_matrix_coef.rds"))) {
    fit_latent = estimator_fn(Y, D_list)
    coef_list = lapply(seq_along(D_list), function(component) {
      result = low_rank_regression_from_cov(
        fit_latent$Sigma_hat[[component]], fun_indices, str_indices,
        latent_matrix_cv$tuning_grids[[component]], cores = CORES, V = fit_latent$V,
        return_factors = TRUE)
      coef_array = array(0, dim = c(p_str, p_str, m))
      for (o in seq_len(m)) {
        sel = latent_matrix_cv$selected[component, o]
        if (sel > 0L) {
          f = result$factors[[o]][[sel]]
          coef_array[, , o] = f$beta1 %*% t(f$beta2)
        }
      }
      coef_array
    })
    saveRDS(coef_list, file.path(RESULT_PATH, "latent_matrix_coef.rds"))
  }

}

if (ARRAY_ID == 8) {

  estimator_fn = function(Y, D_list) {
    vars = matrixStats::colVars(Y)
    w = 1 / vars
    w[vars < 1e-8] = 1
    mvREHE(Y, D_list, return_full = TRUE, w_columns = w)
  }
  N_LAMBDA = 10
  lasso_grid_fn = function(fit, outcomes, covariates) lasso_tuning_grid_fn(fit, outcomes, covariates, n_lambda = N_LAMBDA)

  print("running cv for latent lasso regression tuning selection")
  if (!file.exists(file.path(RESULT_PATH, "latent_lasso_cv.rds"))) {
    set.seed(2026)
    latent_lasso_cv = cv_latent_regression(
      Y, D_list, fam, fit, fun_indices, str_indices,
      tuning_grid_fn = lasso_grid_fn, estimator_fn = estimator_fn,
      regression_fn = lasso_regression_from_cov,
      K = 5, cores = CORES)
    saveRDS(latent_lasso_cv, file.path(RESULT_PATH, "latent_lasso_cv.rds"))
  } else {
    latent_lasso_cv = readRDS(file.path(RESULT_PATH, "latent_lasso_cv.rds"))
  }

  m = length(fun_indices)
  p = length(str_indices)

  if (!file.exists(file.path(RESULT_PATH, "latent_lasso_coef.rds"))) {
    coef_list = lapply(seq_along(D_list), function(component) {
      beta = lasso_regression_from_cov(
        fit$Sigma_hat[[component]], fun_indices, str_indices,
        latent_lasso_cv$tuning_grids[[component]], cores = CORES)
      coef_mat = matrix(0, p, m)
      for (o in seq_len(m)) {
        sel = latent_lasso_cv$selected[component, o]
        if (sel > 0L) coef_mat[, o] = beta[, o, sel]
      }
      coef_mat
    })
    saveRDS(coef_list, file.path(RESULT_PATH, "latent_lasso_coef.rds"))
  }

}

if (ARRAY_ID == 9) {

  estimator_fn = function(Y, D_list) {
    vars = matrixStats::colVars(Y)
    w = 1 / vars
    w[vars < 1e-8] = 1
    mvREHE(Y, D_list, return_full = TRUE, w_columns = w)
  }
  N_LAMBDA = 10
  ridge_grid_fn = function(fit, outcomes, covariates) ridge_tuning_grid_fn(fit, outcomes, covariates, n_lambda = N_LAMBDA)

  print("running cv for latent ridge regression tuning selection")
  if (!file.exists(file.path(RESULT_PATH, "latent_ridge_cv.rds"))) {
    set.seed(2026)
    latent_ridge_cv = cv_latent_regression(
      Y, D_list, fam, fit, fun_indices, str_indices,
      tuning_grid_fn = ridge_grid_fn, estimator_fn = estimator_fn,
      regression_fn = ridge_regression_from_cov,
      K = 5, cores = CORES)
    saveRDS(latent_ridge_cv, file.path(RESULT_PATH, "latent_ridge_cv.rds"))
  } else {
    latent_ridge_cv = readRDS(file.path(RESULT_PATH, "latent_ridge_cv.rds"))
  }

  m = length(fun_indices)
  p = length(str_indices)

  if (!file.exists(file.path(RESULT_PATH, "latent_ridge_coef.rds"))) {
    coef_list = lapply(seq_along(D_list), function(component) {
      beta = ridge_regression_from_cov(
        fit$Sigma_hat[[component]], fun_indices, str_indices,
        latent_ridge_cv$tuning_grids[[component]], cores = CORES)
      coef_mat = matrix(0, p, m)
      for (o in seq_len(m)) {
        sel = latent_ridge_cv$selected[component, o]
        if (sel > 0L) coef_mat[, o] = beta[, o, sel]
      }
      coef_mat
    })
    saveRDS(coef_list, file.path(RESULT_PATH, "latent_ridge_coef.rds"))
  }

}