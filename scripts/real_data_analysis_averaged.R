library(tidyverse)
library(mvREHE)

source("scripts/latent_matrix_regression.R")
source("scripts/latent_lasso_regression.R")
source("scripts/latent_ridge_regression.R")

#################################
##########  Functions ###########
#################################

from_conn_to_vec = function(conn)
{
  as.vector(conn[lower.tri(conn, diag = TRUE)])
}

from_vec_to_conn = function(vec_conn)
{
  p = round(uniroot(function(x) x^2 + x - 2 * length(vec_conn), interval = c(0, 500))$root)
  conn = matrix(NA, nrow = p, ncol = p)
  conn[lower.tri(conn, diag = TRUE)] = vec_conn # Lower triangle
  t_conn = t(conn)
  t_conn[lower.tri(t_conn, diag = TRUE)] = vec_conn # Upper triangle
  t(t_conn)
}

plot_connectome_vec = function(connectome_vec, title, groups, community = FALSE, breaks = NULL, colors = NULL, legend = FALSE, upper_triangle = TRUE, is_full = TRUE, cluster = FALSE) {
  if (length(connectome_vec) == 0) {
    return(NA)
  }
  if (!is_full) {
    connectome = from_vec_to_conn(connectome_vec)
  } else {
    connectome = matrix(connectome_vec, sqrt(length(connectome_vec)))
  }
  groups = factor(groups)
  order = order(groups)
  connectome = connectome[order, order]
  groups = groups[order]
  if (!upper_triangle) {
    connectome[upper.tri(connectome)] = NA
  }
  if (is.null(breaks)) {
    breaks = sort(c(-quantile(abs(connectome), 0.99, na.rm = TRUE), 0, quantile(abs(connectome), 0.99, na.rm = TRUE)))
    colors = c("blue", "white", "red")
  }
  if (community) {
    P_community = t(sapply(unique(groups), FUN = function(x) as.numeric(groups == x)/sum(groups == x)))
    connectome = P_community %*% connectome %*% t(P_community)
    groups = unique(groups)
  }
  heatmap = ComplexHeatmap::Heatmap(
    connectome,
    col = circlize::colorRamp2(breaks, colors),
    cluster_rows = cluster,
    cluster_columns = cluster,
    row_split = groups,
    column_split = groups,
    row_title_gp = grid::gpar(fontsize = 5),
    column_title_gp = grid::gpar(fontsize = 5),
    column_title_rot = 90,
    row_title_rot = 0,
    show_heatmap_legend = legend,
    show_row_names = FALSE,
    show_column_names = FALSE,
    na_col = "gray99",
    heatmap_legend_param = list(direction = "horizontal", title = "", legend_height = unit(2, "cm"), labels_gp = grid::gpar(fontsize = 5), title_gp = grid::gpar(fontsize = 5))
  )
  grid::grid.grabExpr(ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom", column_title = title, column_title_gp = grid::gpar(fontsize = 8)))
}

sqrt_matrix = function(A) {
  eig = eigen(A)
  eig$vec %*% diag(sqrt(pmax(eig$val, 0))) %*% t(eig$vec)
}

#################################
##########  Load data ###########
#################################

DATA_PATH = "Data_YL"
RESULT_PATH = "data_analysis_3_components_10ROI"
dir.create(RESULT_PATH, recursive = TRUE)

data = readRDS(file.path(DATA_PATH, "clean_data.rds"))

fun_connectomes = data$fun_connectomes_averaged
str_connectomes = data$str_connectomes_averaged
K_G = data$K_G
X = data$X
groups = colnames(fun_connectomes[[1]])

Y_fun = sapply(fun_connectomes, c) %>% t
Y_str = sapply(str_connectomes, c) %>% t
colnames(Y_fun) = paste0("fun", 1:ncol(Y_fun))
colnames(Y_str) = paste0("str", 1:ncol(Y_str))

Y = cbind(Y_fun, Y_str)

fun_indices = grep("fun", colnames(Y))
str_indices = grep("str", colnames(Y))

connection_names = outer(groups, groups, "paste")
# connection_names = connection_names[lower.tri(connection_names, diag = TRUE)]

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

fit = mvREHE::mvREHE(Y, D_list = D_list, return_full = T)
for (k in 1:length(fit$Sigma_hat)) {
  attr(fit$Sigma_hat[[k]], "sqrt") = sqrt_matrix(fit$Sigma_hat[[k]])
}

saveRDS(fit, file.path(RESULT_PATH, "fit.rds"))

latent_matrix_regression_fit = cv_latent_matrix_regression(Y, D_list, fun_indices, str_indices, 1:3, 10^seq(2, -3, length.out = 25), function(Y, D_list) mvREHE(Y, D_list, return_full = TRUE), 2, n_init = 25, cores = 40)
latent_ridge_regression_fit = cv_latent_ridge_regression(Y, D_list, fit, fun_indices, str_indices, 100, mvREHE, 2, cores = 40)
latent_lasso_regression_fit = cv_latent_lasso_regression(Y, D_list, fit, fun_indices, str_indices, 100, mvREHE, 2, cores = 40)

raw_matrix_regression_fit = cv_raw_matrix_regression(Y, fun_indices, str_indices, 1:3, 10^seq(2, -3, length.out = 25), 2, n_init = 25, cores = 40)
raw_ridge_regression_fit = cv_raw_ridge_regression(Y, fun_indices, str_indices, 100, 2, cores = 40)
raw_lasso_regression_fit = cv_raw_lasso_regression(Y, fun_indices, str_indices, 100, 2, cores = 40)

r2_components = array(c(latent_matrix_regression_fit$cv_r2, latent_lasso_regression_fit$cv_r2, latent_lasso_regression_fit$cv_r2), dim = c(3, ncol(latent_lasso_regression_fit$cv_r2), 3))
r2_components = apply(r2_components, c(1, 2), max)
r2_components = lapply(1:3, function(i) r2_components[i, ])

r2_raw = array(c(raw_matrix_regression_fit$cv_r2, raw_lasso_regression_fit$cv_r2, raw_lasso_regression_fit$cv_r2), dim = c(ncol(latent_lasso_regression_fit$cv_r2), 3))
r2_raw = apply(r2_raw, 1, max)

community = FALSE
figure = list()
figure[[1]] = plot_connectome_vec(r2_raw, bquote("Observed"~R^2~"\n(Min" == .(round(min(r2_raw, na.rm = TRUE), 2))~", Max" == .(round(max(r2_raw, na.rm = TRUE), 2))~")"), groups = groups, community = community, breaks = c(0, 0.2), colors = c("white", "red"), upper_triangle = TRUE)
figure[[2]] = plot_connectome_vec(r2_components[[genetic_component]], bquote("Genetic"~R^2~"\n(Min" == .(round(min(r2_components[[genetic_component]], na.rm = TRUE), 2))~", Max" == .(round(max(r2_components[[genetic_component]], na.rm = TRUE), 2))~")"), groups = groups, community = community, breaks = c(0, 0.2), colors = c("white", "red"), upper_triangle = TRUE)
figure[[3]] = plot_connectome_vec(r2_components[[common_env_component]], bquote("Common Env"~R^2~"\n(Min" == .(round(min(r2_components[[common_env_component]], na.rm = TRUE), 2))~", Max" == .(round(max(r2_components[[common_env_component]], na.rm = TRUE), 2))~")"), groups = groups, community = community, breaks = c(0, 0.2), colors = c("white", "red"), upper_triangle = TRUE)
figure[[4]] = plot_connectome_vec(r2_components[[unique_env_component]], bquote("Unique Env"~R^2~"\n(Min" == .(round(min(r2_components[[unique_env_component]], na.rm = TRUE), 2))~", Max" == .(round(max(r2_components[[unique_env_component]], na.rm = TRUE), 2))~")"), groups = groups, community = community, breaks = c(0, 0.2), colors = c("white", "red"), upper_triangle = TRUE)

pdf(file.path(RESULT_PATH, "r2.pdf"), height = 2.9, width = 10)
patchwork::wrap_plots(figure, ncol = 4, byrow = TRUE) +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c", "d")))
dev.off()

max_index = which.max(r2_components[[2]])
connection_names[max_index]

latent_matrix_fit = latent_matrix_regression(fit, 2, max_index, 101:200, latent_matrix_regression_fit$rank[2, max_index], latent_matrix_regression_fit$lambda[2, max_index])
latent_ridge_fit = latent_ridge_regression(fit, 2, max_index, 101:200, latent_ridge_regression_fit$lambda[2, max_index])
latent_lasso_fit = latent_lasso_regression(fit, 2, max_index, 101:200, latent_lasso_regression_fit$lambda[2, max_index])

groups = factor(groups, levels = groups)

community = FALSE
figure = list()
figure[[1]] = plot_connectome_vec(c(tcrossprod(latent_matrix_fit$beta1, latent_matrix_fit$beta2)), "", groups = groups, community = community, breaks = c(-1, 0, 1), colors = c("blue", "white", "red"), is_full = TRUE)
figure[[2]] = plot_connectome_vec(get_lower_tri_beta(latent_lasso_fit), "", groups = groups, community = community, breaks = c(-1, 0, 1), colors = c("blue", "white", "red"), upper_triangle = FALSE, is_full = FALSE)
figure[[3]] = plot_connectome_vec(get_lower_tri_beta(latent_ridge_fit), "", groups = groups, community = community, breaks = c(-1, 0, 1), colors = c("blue", "white", "red"), upper_triangle = FALSE, is_full = FALSE)

pdf(file.path(RESULT_PATH, "beta_G.pdf"), height = 2.9, width = 7.5)
patchwork::wrap_plots(figure, ncol = 3, byrow = TRUE) +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c")))
dev.off()

get_lower_tri_beta = function(beta) {
  
  beta = matrix(beta, sqrt(length(beta)))
  beta = beta + t(beta) - diag(beta)
  beta[lower.tri(beta, diag = TRUE)]
  
}
