library(ggplot2)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(grid)

from_vec_to_conn = function(vec_conn)
{
  p = round(uniroot(function(x) x^2 + x - 2 * length(vec_conn), interval = c(0, 500))$root)
  conn = matrix(NA, nrow = p, ncol = p)
  conn[lower.tri(conn, diag = TRUE)] = vec_conn # Lower triangle
  t_conn = t(conn)
  t_conn[lower.tri(t_conn, diag = TRUE)] = vec_conn # Upper triangle
  t(t_conn)
}

plot_connectome_vec = function(connectome_vec, title, groups, community = FALSE, breaks = NULL, colors = NULL, legend = FALSE, upper_triangle = TRUE, is_full = FALSE, cluster = FALSE) {
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

RESULT_PATH = "real_data_analysis"

latent_ridge_regression_fit = readRDS(file.path(RESULT_PATH, paste0("latent_ridge_regression_fit",".rds")))
raw_ridge_regression_fit = readRDS(file.path(RESULT_PATH, paste0("raw_ridge_regression_fit",".rds")))

latent_lasso_regression_fit = readRDS(file.path(RESULT_PATH, paste0("latent_lasso_regression_fit",".rds")))
raw_lasso_regression_fit = readRDS(file.path(RESULT_PATH, paste0("raw_lasso_regression_fit",".rds")))

RESULT_PATH = "data_analysis_3_components_68ROI/"

latent_matrix_regression_fit = list(cv_r2 = matrix(0, nrow = 3, ncol = ncol(latent_ridge_regression_fit$cv_r2)), lambda = matrix(0, nrow = 3, ncol = ncol(latent_ridge_regression_fit$cv_r2)), rank = matrix(0, nrow = 3, ncol = ncol(latent_ridge_regression_fit$cv_r2)))
raw_matrix_regression_fit = list(cv_r2 = numeric(ncol(latent_ridge_regression_fit$cv_r2)), lambda = numeric(ncol(latent_ridge_regression_fit$cv_r2)), rank = numeric(ncol(latent_ridge_regression_fit$cv_r2)))
for (i in 1:ncol(latent_matrix_regression_fit$cv_r2)) {
  if (file.exists(file.path(RESULT_PATH, paste0("latent_matrix_regression_fit_", i, ".rds")))) {
    print(file.path(RESULT_PATH, paste0("latent_matrix_regression_fit_", i, ".rds")))
    fit = readRDS(file.path(RESULT_PATH, paste0("latent_matrix_regression_fit_", i, ".rds")))
    latent_matrix_regression_fit$cv_r2[, i] = fit$cv_r2
    latent_matrix_regression_fit$lambda[, i] = fit$lambda
    latent_matrix_regression_fit$rank[, i] = fit$rank
  } else {
    cat(i, ",")
  }
}
for (i in 1:length(raw_matrix_regression_fit$cv_r2)) {
  if (file.exists(file.path(RESULT_PATH, paste0("raw_matrix_regression_fit_", i, ".rds")))) {
    fit = readRDS(file.path(RESULT_PATH, paste0("raw_matrix_regression_fit_", i, ".rds")))
    raw_matrix_regression_fit$cv_r2[i] = fit$cv_r2
    raw_matrix_regression_fit$lambda[i] = fit$lambda
    raw_matrix_regression_fit$rank[i] = fit$rank
  } else {
    print(i)
  }
}

r2_components = array(c(latent_matrix_regression_fit$cv_r2, latent_ridge_regression_fit$cv_r2, latent_lasso_regression_fit$cv_r2), dim = c(3, ncol(latent_ridge_regression_fit$cv_r2), 3))
r2_components = apply(r2_components, c(1, 2), max)
r2_components = lapply(1:3, function(i) r2_components[i, ])

r2_raw = array(c(raw_matrix_regression_fit$cv_r2, raw_ridge_regression_fit$cv_r2, raw_lasso_regression_fit$cv_r2), dim = c(ncol(latent_ridge_regression_fit$cv_r2), 3))
r2_raw = apply(r2_raw, c(1), max)

genetic_component = 2
common_env_component = 3
unique_env_component = 1

data = readRDS("Data_YL/clean_data.rds")
fun_groups = str_groups = data$groups

min1 = function(x) min(x, na.rm = TRUE)
max1 = function(x) max(x, na.rm = TRUE)

community = FALSE
figure = list()
figure[[1]] = plot_connectome_vec(r2_raw, bquote("Observed"~R^2~"\n(min" == .(round(min1(r2_raw), 2))~", max" == .(round(max1(r2_raw), 2))~")"), groups = fun_groups, community = community, breaks = c(0, 0.35), colors = c("white", "red"))
figure[[2]] = plot_connectome_vec(r2_components[[genetic_component]], bquote("Genetic"~R^2~"\n(min" == .(round(min1(r2_components[[genetic_component]]), 2))~", max" == .(round(max1(r2_components[[genetic_component]]), 2))~")"), groups = fun_groups, community = community, breaks = c(0, 0.35), colors = c("white", "red"))
figure[[3]] = plot_connectome_vec(r2_components[[common_env_component]], bquote("Common Env"~R^2~"\n(min" == .(round(min1(r2_components[[common_env_component]]), 2))~", max" == .(round(max1(r2_components[[common_env_component]]), 2))~")"), groups = fun_groups, community = community, breaks = c(0, 0.35), colors = c("white", "red"))
figure[[4]] = plot_connectome_vec(r2_components[[unique_env_component]], bquote("Unique Env"~R^2~"\n(min" == .(round(min1(r2_components[[unique_env_component]]), 2))~", max" == .(round(max1(r2_components[[unique_env_component]]), 2))~")"), groups = fun_groups, community = community, breaks = c(0, 0.35), colors = c("white", "red"))

pdf(file.path(RESULT_PATH, "r2.pdf"), height = 2.9, width = 10)
patchwork::wrap_plots(figure, ncol = 4, byrow = TRUE) +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c", "d")))
dev.off()

source("scripts/latent_lasso_regression.R")
source("scripts/latent_ridge_regression.R")
source("scripts/latent_matrix_regression.R")

RESULT_PATH = "real_data_analysis"
fit_reduced = readRDS(file.path(RESULT_PATH, "fit_reduced.rds"))
fit = readRDS(file.path(RESULT_PATH, "fit.rds"))

groups = fun_groups

connection_names = outer(groups, groups, "paste")
connection_names = connection_names[lower.tri(connection_names, diag = TRUE)]

max_index = which.max(r2_components[[2]])
connection_names[max_index]
latent_matrix_fit = latent_matrix_regression(fit, 2, max_index, 2347:4692, latent_matrix_regression_fit$rank[2, max_index], latent_matrix_regression_fit$lambda[2, max_index], tolerance = 1e-3)
latent_matrix_fit = c(tcrossprod(latent_matrix_fit$beta1, latent_matrix_fit$beta2))
latent_ridge_fit = latent_ridge_regression(fit, 2, max_index, 2347:4692, latent_ridge_regression_fit$lambda[2, max_index])
latent_lasso_fit = latent_lasso_regression(fit, 2, max_index, 2347:4692, latent_lasso_regression_fit$lambda[2, max_index])
plot(latent_ridge_fit, latent_lasso_fit)


community = FALSE
figure = list()
figure[[1]] = plot_connectome_vec(latent_matrix_fit, "Genetic Matrix Regression Coefficients", groups = groups, community = community, is_full = TRUE, legend = TRUE)
figure[[2]] = plot_connectome_vec(latent_lasso_fit, "Genetic Lasso Regression Coefficients", groups = groups, community = community, upper_triangle = FALSE, legend = TRUE)
figure[[3]] = plot_connectome_vec(latent_ridge_fit, "Genetic Ridge Regression Coefficients", groups = groups, community = community, upper_triangle = FALSE, legend = TRUE)

pdf(file.path(RESULT_PATH, "beta_G.pdf"), height = 2.9, width = 7)
patchwork::wrap_plots(figure, ncol = 3, byrow = TRUE)  +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c")))
dev.off()
