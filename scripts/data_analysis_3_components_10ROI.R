rm(list = ls())

library(ggplot2)
library(gridExtra)

library(R.matlab)
library(readr)
library(Matrix)

library(dplyr)
library(tidyr)
library(devtools)

if (!require("mvREHE")) install_local("Packages/mvREHE")

source("scripts/latent_regression.R")

#################################
##########  Functions ###########
#################################

from_conn_to_vec = function(conn)
{
  as.vector(conn[lower.tri(conn, diag = TRUE)])
}

from_vec_to_conn = function(vec_conn, upper_triangle = FALSE)
{
  p = round(uniroot(function(x) x^2 + x - 2 * length(vec_conn), interval = c(0, 500))$root)
  conn = matrix(NA, nrow = p, ncol = p)
  conn[lower.tri(conn, diag = TRUE)] = vec_conn # Lower triangle
  if (upper_triangle) {
    t_conn = t(conn)
    t_conn[lower.tri(t_conn, diag = TRUE)] = vec_conn # Upper triangle
    t(t_conn)
  } else {
    conn
  }
}

plot_connectome_vec = function(connectome_vec, title, groups, community = FALSE, breaks = NULL, colors = NULL, legend = FALSE, upper_triangle = FALSE, cluster = FALSE) {
  if (length(connectome_vec) == 0) {
    return(NA)
  }
  connectome = from_vec_to_conn(connectome_vec, upper_triangle)
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
    heatmap_legend_param = list(title = expression("mvREHE"~hat(h)[j]^2), legend_height = unit(2, "cm"), labels_gp = grid::gpar(fontsize = 5), title_gp = grid::gpar(fontsize = 5))
  )
  grid::grid.grabExpr(ComplexHeatmap::draw(heatmap, column_title = title, column_title_gp = grid::gpar(fontsize = 8)))
}

#################################
##########  Load data ###########
#################################

# Load Kinship matrix
kinship = readMat('data/kinship.mat'); # This is 2*K in the Solar-Eclipse notation
K_G = as(kinship$K[[1]], "TsparseMatrix"); # Kinship matrix

id_subjects = as.character(kinship$K[[2]]);   # ids of the subjects
rownames(K_G) = id_subjects;
colnames(K_G) = id_subjects;

session = 1 # There are four fmri sessions for each subject: we pick the first one

# Load functional connectome data (how the different parts of the brain are functionally connected)

if (file.exists("data/Data_YL/fun_connectomes.rds")) {

  fun_connectomes = readRDS("data/Data_YL/fun_connectomes.rds")

} else {

  regions_group_fun = read_csv("data/name_regions_alternating.csv", col_names = T, col_types = cols())
  regions_group_fun$sorted_idx_gordon #indices ROIs re-ordered to cluster ROIs by macro-regions
  regions_group_fun$Var3 # name macro-region each ROI belongs to

  indices = regions_group_fun$sorted_idx_gordon
  groups = regions_group_fun$Var3

  fun_connectomes = list()

  # It can take long so you may want to load only the first few
  for (i in seq(id_subjects))
  {
    tryCatch(
      {
        connectome = as.matrix(read_csv(paste0('data/fun_glasser/3T_HCP1200_MSMAll_glasser_et_al_conn/',
                                               id_subjects[i], '_', session, '_cov.csv'), col_names = F, col_types = cols()))
        connectome = connectome[indices, indices]
        rownames(connectome) = colnames(connectome) = groups
        fun_connectomes[[i]] = connectome
      },
      warning = function(cond){
        fun_connectomes[[i]] = NULL
      },
      error = function(e){
        fun_connectomes[[i]] = NULL
      })
  }

  names(fun_connectomes) = id_subjects

  saveRDS(fun_connectomes, "data/fun_connectomes.rds")

}

# Load structural connectome data (how the different parts of the brain are structurally connected)

if (file.exists("data/Data_YL/str_connectomes.rds")) {

  str_connectomes = readRDS("data/Data_YL/str_connectomes.rds")

} else {

  regions_group_str = regions_group_fun[regions_group_fun$sorted_idx_gordon<=180, ]

  indices = regions_group_str$sorted_idx_gordon # indices structural ROIs re-ordered to cluster ROIs by macro-regions
  groups = regions_group_str$Var3 # name macro-region each structural ROI belongs to

  str_connectomes = list()

  # It can take long so you may want to load only the first few
  for (i in seq(id_subjects))
  {

    tryCatch(
      {
        str_raw = read_delim(paste0('data/str_glasser/sub-',
                                    id_subjects[i], '_ses-', session, '_run-1_dwi_Glasser_space-MNI152NLin6_res-1x1x1_connectome.csv'), col_names = F, col_types = cols(), delim = " ")
        str_raw_spmat = sparseMatrix(str_raw$X1, str_raw$X2, x = str_raw$X3, symmetric = TRUE)
        connectome = as.matrix(str_raw_spmat)
        connectome = connectome[indices, indices]
        rownames(connectome) = colnames(connectome) = groups
        str_connectomes[[i]] = connectome
      },
      warning = function(cond){
        str_connectomes[[i]] = NULL
      },
      error = function(e){
        str_connectomes[[i]] = NULL
      })
  }

  names(str_connectomes) = id_subjects

  saveRDS(str_connectomes, "data/str_connectomes.rds")

}

#######################
###### Join data ######
#######################

## This is to avoid using inner_join(by = "subject") https://stackoverflow.com/questions/64692005/how-do-i-split-a-data-frame-then-apply-inner-join
common_subjects = intersect(intersect(id_subjects, names(fun_connectomes)[!sapply(fun_connectomes, is.null)]), names(str_connectomes)[!sapply(str_connectomes, is.null)])
common_subjects = setdiff(common_subjects, id_subjects[366])
# keep only rows with common_subjects in each data frame

K_G = K_G[common_subjects, common_subjects] # Kinship matrix
rownames(K_G) = as.character(common_subjects);
colnames(K_G) = as.character(common_subjects);

K_G = as.matrix(K_G)
order = hclust(as.dist(-K_G))$order

fun_connectomes = fun_connectomes[common_subjects[order]]
str_connectomes = str_connectomes[common_subjects[order]]

K_G = K_G[common_subjects[order], common_subjects[order]]

#################################
##########  Analysis ############
#################################

X = read.csv("data/conf.csv")
rownames(X) = X$subject
X = X[common_subjects[order], c("Age", "Age.2", "Sex", "FS_IntraCranial_Vol..1.3.", "FS_BrainSeg_Vol..1.3.")]
X = cbind(rep(1, nrow(X)), X)

RESULT_PATH = "data_analysis_3_components_68ROI"
dir.create(RESULT_PATH, recursive = TRUE)

groups = read.csv("data/communities_aparc.csv")
groups = groups[order(groups$id_ROI), ][, 3]

# groups = colnames(fun_connectomes[[1]])
P_community = t(sapply(unique(groups), FUN = function(x) as.numeric(groups == x)/sum(groups == x)))
fun_connectomes = lapply(fun_connectomes, function(x) {
  connectome = P_community %*% x %*% t(P_community)
  rownames(connectome) = colnames(connectome) = unique(groups)
  connectome[sort(rownames(connectome)), sort(rownames(connectome))]
})

# groups = colnames(str_connectomes[[1]])
P_community = t(sapply(unique(groups), FUN = function(x) as.numeric(groups == x)/sum(groups == x)))
str_connectomes = lapply(str_connectomes, function(x) {
  connectome = P_community %*% x %*% t(P_community)
  rownames(connectome) = colnames(connectome) = unique(groups)
  connectome
  connectome[sort(rownames(connectome)), sort(rownames(connectome))]
})

Y_fun = sapply(fun_connectomes, from_conn_to_vec) %>% t
Y_str = sapply(str_connectomes, from_conn_to_vec) %>% t
colnames(Y_fun) = paste0("fun", 1:ncol(Y_fun))
colnames(Y_str) = paste0("str", 1:ncol(Y_str))

Y = cbind(Y_fun, Y_str)
Y_mean = colMeans(Y)

fun_indices = grep("fun", colnames(Y))
str_indices = grep("str", colnames(Y))

fun_groups = str_groups = unique(groups)

# fun_groups = colnames(fun_connectomes[[1]])
# str_groups = colnames(str_connectomes[[1]])

connection_names = outer(fun_groups, fun_groups, "paste")
connection_names = connection_names[lower.tri(connection_names, diag = TRUE)]

fixed_effects = TRUE
residuals = lsfit(X, Y, intercept = FALSE)$residuals
colnames(residuals) = colnames(Y)
Y = residuals

scale_type = "column"
Y = scale(Y)
Y[, attr(Y, "scaled:scale") == 0] = 0

D_list = list(diag(nrow = nrow(Y), ncol = nrow(Y)), as.matrix(K_G), as.matrix(K_G > 0) * 1)
genetic_component = 2
common_env_component = 3
unique_env_component = 1

# Variance component model estimation (replace code)
start.time = Sys.time()
fit = mvREHE::mvREHE(Y, D_list = D_list, return_full = T)
end.time = Sys.time()
time.taken = end.time - start.time
print(time.taken)

saveRDS(fit, file.path(RESULT_PATH, "fit.rds"))

raw_matrix_regression_fit = cv_raw_matrix_regression(Y, fun_indices, str_indices, 1:5, 10^seq(2, -3, length.out = 100), 2, n_init = 10, cores = 150)
latent_matrix_regression_fit = cv_latent_matrix_regression(Y, D_list, fun_indices, str_indices, 1:5, 10^seq(2, -3, length.out = 100), mvREHE, 2, n_init = 10, cores = 150)

saveRDS(list(latent = latent_matrix_regression_fit, raw = raw_matrix_regression_fit), file.path(RESULT_PATH, "regression_fit.rds"))

