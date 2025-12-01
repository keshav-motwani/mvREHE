library(tidyverse)
library(mvREHE)

source("scripts/latent_matrix_regression.R")
source("scripts/latent_ridge_regression.R")
source("scripts/latent_lasso_regression.R")

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

DATA_PATH = "Data_YL"
RESULT_PATH = "real_data_analysis"
dir.create(RESULT_PATH, recursive = TRUE)

data = readRDS(file.path(DATA_PATH, "clean_data.rds"))

fun_connectomes = data$fun_connectomes
str_connectomes = data$str_connectomes
K_G = data$K_G
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

D_list = list(diag(nrow = nrow(Y), ncol = nrow(Y)), as.matrix(K_G), as.matrix(K_G > 0) * 1)
genetic_component = 2
common_env_component = 3
unique_env_component = 1


OUTCOME = as.numeric(commandArgs(trailingOnly=TRUE)[1])
print(OUTCOME)

if (OUTCOME > 0) {
  
  if (!file.exists(file.path(RESULT_PATH, paste0("latent_matrix_regression_fit_", OUTCOME,".rds")))) {
    latent_matrix_regression_fit = cv_latent_matrix_regression(Y, D_list, fun_indices[OUTCOME], str_indices, 10:20, 10^seq(2, -3, length.out = 10), function(Y, D_list) mvREHE(Y, D_list, return_full = F), 2, n_init = 5, cores = 50, tolerance = 1e-3)
    saveRDS(latent_matrix_regression_fit, file.path(RESULT_PATH, paste0("latent_matrix_regression_fit_", OUTCOME,".rds")))
  }
  if (!file.exists(file.path(RESULT_PATH, paste0("raw_matrix_regression_fit_", OUTCOME,".rds")))) {
    raw_matrix_regression_fit = cv_raw_matrix_regression(Y, fun_indices[OUTCOME], str_indices, 10:20, 10^seq(2, -3, length.out = 10), 2, n_init = 5, cores = 50, tolerance = 1e-3)
    saveRDS(raw_matrix_regression_fit, file.path(RESULT_PATH, paste0("raw_matrix_regression_fit_", OUTCOME,".rds")))
  }
  
} else {
  
  if (!file.exists(file.path(RESULT_PATH, "fit.rds"))) {
    fit = mvREHE::mvREHE(Y, D_list = D_list, return_full = F)
    saveRDS(fit, file.path(RESULT_PATH, "fit_reduced.rds"))
    fit$Sigma_hat = lapply(fit$Sigma_hat, function(Sigma) fit$V %*% Sigma %*% t(fit$V))
    fit$V = NULL
    for (k in 1:length(fit$Sigma_hat)) {
      attr(fit$Sigma_hat[[k]], "sqrt") = sqrt_matrix(fit$Sigma_hat[[k]])
    }
    saveRDS(fit, file.path(RESULT_PATH, "fit.rds"))
  } else {
    fit = readRDS(file.path(RESULT_PATH, "fit.rds"))
  }
  
  if (!file.exists(file.path(RESULT_PATH, "REML_fit.rds"))) {
    source("scripts/mvREML.R")
    REML_fit = univariate(Y, D_list, mvREML)
    saveRDS(REML_fit, file.path(RESULT_PATH, "REML_fit.rds"))
  }
  
  if (!file.exists(file.path(RESULT_PATH, paste0("latent_ridge_regression_fit",".rds")))) {
    latent_ridge_regression_fit = cv_latent_ridge_regression(Y, D_list, fit, fun_indices, str_indices, 100, mvREHE, 2, cores = 108)
    saveRDS(latent_ridge_regression_fit, file.path(RESULT_PATH, paste0("latent_ridge_regression_fit",".rds")))
  }
  if (!file.exists(file.path(RESULT_PATH, paste0("raw_ridge_regression_fit",".rds")))) {
    raw_ridge_regression_fit = cv_raw_ridge_regression(Y, fun_indices, str_indices, 100, 2, cores = 108)
    saveRDS(raw_ridge_regression_fit, file.path(RESULT_PATH, paste0("raw_ridge_regression_fit",".rds")))
  }
  
  if (!file.exists(file.path(RESULT_PATH, paste0("latent_lasso_regression_fit",".rds")))) {
    latent_lasso_regression_fit = cv_latent_lasso_regression(Y, D_list, fit, fun_indices, str_indices, 100, mvREHE, 2, cores = 108)
    saveRDS(latent_lasso_regression_fit, file.path(RESULT_PATH, paste0("latent_lasso_regression_fit",".rds")))
  }
  if (!file.exists(file.path(RESULT_PATH, paste0("raw_lasso_regression_fit",".rds")))) {
    raw_lasso_regression_fit = cv_raw_lasso_regression(Y, fun_indices, str_indices, 100, 2, cores = 108)
    saveRDS(raw_lasso_regression_fit, file.path(RESULT_PATH, paste0("raw_lasso_regression_fit",".rds")))
  }
  
}
