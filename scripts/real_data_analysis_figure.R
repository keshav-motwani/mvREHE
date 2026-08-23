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

RESULT_PATH = "real_data_analysis_ica_260804"

latent_matrix_regression_nested = readRDS(file.path(RESULT_PATH, paste0("latent_matrix_regression_nested",".rds")))
raw_matrix_regression_nested = readRDS(file.path(RESULT_PATH, paste0("raw_matrix_regression_nested",".rds")))


latent_ridge_regression_nested = readRDS(file.path(RESULT_PATH, paste0("latent_ridge_regression_nested",".rds")))
raw_ridge_regression_nested = readRDS(file.path(RESULT_PATH, paste0("raw_ridge_regression_nested",".rds")))

latent_lasso_regression_nested = readRDS(file.path(RESULT_PATH, paste0("latent_lasso_regression_nested",".rds")))
raw_lasso_regression_nested = readRDS(file.path(RESULT_PATH, paste0("raw_lasso_regression_nested",".rds")))

r2_components = array(c(latent_matrix_regression_nested$cv_r2, latent_ridge_regression_nested$cv_r2, latent_lasso_regression_nested$cv_r2), dim = c(3, ncol(latent_ridge_regression_nested$cv_r2), 3))
r2_components = apply(r2_components, c(1, 2), max)
r2_components = lapply(1:3, function(i) r2_components[i, ])

r2_raw = array(c(raw_matrix_regression_nested$cv_r2, raw_ridge_regression_nested$cv_r2, raw_lasso_regression_nested$cv_r2), dim = c(ncol(latent_ridge_regression_nested$cv_r2), 3))
r2_raw = apply(r2_raw, c(1), max)

genetic_component = 2
common_env_component = 3
unique_env_component = 1

data = readRDS("Data_YL/clean_data.rds")
groups = fun_groups = str_groups = data$groups

pct_label = function(x) {
  q = round(quantile(x, c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE), 2)
  substitute(min==q1~"/"~Q[1]==q2~"/"~med==q3~"/"~Q[3]==q4~"/"~max==q5,
             list(q1 = q[[1]], q2 = q[[2]], q3 = q[[3]], q4 = q[[4]], q5 = q[[5]]))
}

r2_max = max(c(r2_raw, unlist(r2_components)), na.rm = TRUE)
r2_breaks = c(0, r2_max)
r2_colors = c("white", "red")

community = FALSE
figure = list()
figure[[1]] = plot_connectome_vec(r2_raw, bquote(atop("Observed"~R^2, scriptstyle(.(pct_label(r2_raw))))), groups = fun_groups, community = community, breaks = r2_breaks, colors = r2_colors, legend = TRUE)
figure[[2]] = plot_connectome_vec(r2_components[[genetic_component]], bquote(atop("Genetic"~R^2, scriptstyle(.(pct_label(r2_components[[genetic_component]]))))), groups = fun_groups, community = community, breaks = r2_breaks, colors = r2_colors, legend = TRUE)
figure[[3]] = plot_connectome_vec(r2_components[[common_env_component]], bquote(atop("Common Env"~R^2, scriptstyle(.(pct_label(r2_components[[common_env_component]]))))), groups = fun_groups, community = community, breaks = r2_breaks, colors = r2_colors, legend = TRUE)
figure[[4]] = plot_connectome_vec(r2_components[[unique_env_component]], bquote(atop("Unique Env"~R^2, scriptstyle(.(pct_label(r2_components[[unique_env_component]]))))), groups = fun_groups, community = community, breaks = r2_breaks, colors = r2_colors, legend = TRUE)

pdf(file.path(RESULT_PATH, "r2.pdf"), height = 3.4, width = 10)
patchwork::wrap_plots(figure, ncol = 4, byrow = TRUE) +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c", "d")))
dev.off()

print(file.path(RESULT_PATH, "r2.pdf"))

for (method in c("ridge", "lasso", "matrix")) {
  r2_raw_m = get(paste0("raw_", method, "_regression_nested"))$cv_r2
  cv_r2_m  = get(paste0("latent_", method, "_regression_nested"))$cv_r2
  r2_comp_m = lapply(1:3, function(i) cv_r2_m[i, ])

  r2_max_m = max(c(r2_raw_m, unlist(r2_comp_m)), na.rm = TRUE)
  r2_breaks_m = c(0, r2_max_m)

  fig_m = list()
  fig_m[[1]] = plot_connectome_vec(r2_raw_m, bquote(atop("Observed"~R^2, scriptstyle(.(pct_label(r2_raw_m))))), groups = fun_groups, community = community, breaks = r2_breaks_m, colors = r2_colors, legend = TRUE)
  fig_m[[2]] = plot_connectome_vec(r2_comp_m[[genetic_component]], bquote(atop("Genetic"~R^2, scriptstyle(.(pct_label(r2_comp_m[[genetic_component]]))))), groups = fun_groups, community = community, breaks = r2_breaks_m, colors = r2_colors, legend = TRUE)
  fig_m[[3]] = plot_connectome_vec(r2_comp_m[[common_env_component]], bquote(atop("Common Env"~R^2, scriptstyle(.(pct_label(r2_comp_m[[common_env_component]]))))), groups = fun_groups, community = community, breaks = r2_breaks_m, colors = r2_colors, legend = TRUE)
  fig_m[[4]] = plot_connectome_vec(r2_comp_m[[unique_env_component]], bquote(atop("Unique Env"~R^2, scriptstyle(.(pct_label(r2_comp_m[[unique_env_component]]))))), groups = fun_groups, community = community, breaks = r2_breaks_m, colors = r2_colors, legend = TRUE)

  pdf(file.path(RESULT_PATH, paste0("r2_", method, ".pdf")), height = 3.4, width = 10)
  print(patchwork::wrap_plots(fig_m, ncol = 4, byrow = TRUE) +
    patchwork::plot_annotation(tag_levels = list(c("a", "b", "c", "d"))))
  dev.off()
  print(file.path(RESULT_PATH, paste0("r2_", method, ".pdf")))
}

connection_names = outer(groups, groups, "paste")
connection_names = connection_names[lower.tri(connection_names, diag = TRUE)]

max_index = which.max(r2_components[[2]])
print(max_index)
print(connection_names[max_index])
print(latent_matrix_regression_nested$cv_r2[2, max_index])
print(latent_lasso_regression_nested$cv_r2[2, max_index])
print(latent_ridge_regression_nested$cv_r2[2, max_index])

cv_r2_methods = rbind(
  ridge = latent_ridge_regression_nested$cv_r2[2, ],
  lasso = latent_lasso_regression_nested$cv_r2[2, ],
  matrix = latent_matrix_regression_nested$cv_r2[2, ]
)

best_method = apply(cv_r2_methods, 2, function(x) names(which.max(x)))
best_method_pct = prop.table(table(factor(best_method, levels = c("ridge", "lasso", "matrix")))) * 100
print(best_method_pct)

latent_matrix_coef = readRDS(file.path(RESULT_PATH, "latent_matrix_coef.rds"))
latent_ridge_coef = readRDS(file.path(RESULT_PATH, "latent_ridge_coef.rds"))
latent_lasso_coef = readRDS(file.path(RESULT_PATH, "latent_lasso_coef.rds"))

community = FALSE
figure = list()
figure[[1]] = plot_connectome_vec(c(latent_matrix_coef[[2]][, , max_index]), "Genetic Tensor Regression Coefficients", groups = groups, community = community, is_full = TRUE, legend = TRUE)
figure[[2]] = plot_connectome_vec(latent_lasso_coef[[2]][, max_index], "Genetic Lasso Regression Coefficients", groups = groups, community = community, upper_triangle = FALSE, legend = TRUE)
figure[[3]] = plot_connectome_vec(latent_ridge_coef[[2]][, max_index], "Genetic Ridge Regression Coefficients", groups = groups, community = community, upper_triangle = FALSE, legend = TRUE)

pdf(file.path(RESULT_PATH, "beta_G.pdf"), height = 3.1, width = 7.5)
patchwork::wrap_plots(figure, ncol = 3, byrow = TRUE)  +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c")))
dev.off()


raw_cor = readRDS(file.path(RESULT_PATH, "raw_cor.rds"))
fit = readRDS(file.path(RESULT_PATH, "fit.rds"))

latent_spatial_coupling = lapply(fit$Sigma_hat, function(Sigma) diag(cov2cor(Sigma)[1:2346, 2347:4692]))
raw_spatial_coupling = diag(raw_cor[1:2346, 2347:4692])

all_sc = c(raw_spatial_coupling, unlist(latent_spatial_coupling))
sc_lim = quantile(abs(all_sc), 0.99, na.rm = TRUE)
sc_breaks = c(-sc_lim, 0, sc_lim)
sc_colors = c("blue", "white", "red")
sc_xlim = range(all_sc, na.rm = TRUE)

community = FALSE
figure = list()
figure[[1]] = plot_connectome_vec(raw_spatial_coupling, bquote(atop("Observed Correlation", scriptstyle(.(pct_label(raw_spatial_coupling))))), groups = groups, community = community, legend = TRUE, breaks = sc_breaks, colors = sc_colors)
figure[[2]] = plot_connectome_vec(latent_spatial_coupling[[2]], bquote(atop("Genetic Correlation", scriptstyle(.(pct_label(latent_spatial_coupling[[2]]))))), groups = groups, community = community, legend = TRUE, breaks = sc_breaks, colors = sc_colors)
figure[[3]] = plot_connectome_vec(latent_spatial_coupling[[3]], bquote(atop("Common Env Correlation", scriptstyle(.(pct_label(latent_spatial_coupling[[3]]))))), groups = groups, community = community, legend = TRUE, breaks = sc_breaks, colors = sc_colors)
figure[[4]] = plot_connectome_vec(latent_spatial_coupling[[1]], bquote(atop("Unique Env Correlation", scriptstyle(.(pct_label(latent_spatial_coupling[[1]]))))), groups = groups, community = community, legend = TRUE, breaks = sc_breaks, colors = sc_colors)

pdf(file.path(RESULT_PATH, "structure_function_correlation.pdf"), height = 3.4, width = 10)
patchwork::wrap_plots(figure, ncol = 4, byrow = TRUE)  +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c", "d")))
dev.off()
