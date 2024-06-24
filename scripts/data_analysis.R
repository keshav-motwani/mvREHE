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

if (file.exists("data/fun_connectomes.rds")) {

  fun_connectomes = readRDS("data/fun_connectomes.rds")

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

  saveRDS(fun_connectomes, "data/fun_connectomes.rds")

}

# Load structural connectome data (how the different parts of the brain are structurally connected)

if (file.exists("data/str_connectomes.rds")) {

  str_connectomes = readRDS("data/str_connectomes.rds")

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

  saveRDS(str_connectomes, "data/str_connectomes.rds")

}

#######################
###### Join data ######
#######################

## This is to avoid using inner_join(by = "subject") https://stackoverflow.com/questions/64692005/how-do-i-split-a-data-frame-then-apply-inner-join
common_subjects = intersect(id_subjects, id_subjects[!sapply(fun_connectomes,is.null)]) %>%
  intersect(id_subjects[!sapply(str_connectomes,is.null)])

# keep only rows with common_subjects in each data frame

K_G = K_G[id_subjects %in% common_subjects,id_subjects %in% common_subjects] # Kinship matrix
rownames(K_G) = as.character(common_subjects);
colnames(K_G) = as.character(common_subjects);

fun_connectomes = fun_connectomes[id_subjects %in% common_subjects]
str_connectomes = str_connectomes[id_subjects %in% common_subjects]

#################################
##########  Analysis ############
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


# Plot rotated PC modes of variation
plot_connectome_vec = function(connectome_vec, title, groups, community = FALSE, breaks = NULL, colors = NULL) {
  if (length(connectome_vec) == 0) {
    return(NA)
  }
  connectome = from_vec_to_conn(connectome_vec)
  if (is.null(breaks)) {
    breaks = sort(c(-quantile(abs(connectome), 0.99), 0, quantile(abs(connectome), 0.99)))
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
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = groups,
    column_split = groups,
    row_title_gp = grid::gpar(fontsize = 5),
    column_title_gp = grid::gpar(fontsize = 5),
    column_title_rot = 90,
    row_title_rot = 0,
    show_heatmap_legend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE
  )
  grid::grid.grabExpr(ComplexHeatmap::draw(heatmap, column_title = title, column_title_gp = grid::gpar(fontsize = 8)))
}

separate_svd_irlba = function(Y, r) {

  V1 = mvREHE:::svd_irlba(Y[, grep("fun", colnames(Y))], r)
  V2 = mvREHE:::svd_irlba(Y[, grep("str", colnames(Y))], r)

  V = as.matrix(Matrix::bdiag(V1, V2))
  attr(V, "groups") = c(rep(1, r), rep(2, r))

  V

}

X = read.csv("data/conf.csv")
rownames(X) = X$subject
X = X[common_subjects, c("Age", "Age.2", "Sex", "FS_IntraCranial_Vol..1.3.", "FS_BrainSeg_Vol..1.3.")]
X = cbind(rep(1, nrow(X)), X)

RESULT_PATH = "data_analysis_results"
dir.create(RESULT_PATH, recursive = TRUE)

for (average in c(F, T)) {

  if (average) {

    groups = colnames(fun_connectomes[[1]])
    P_community = t(sapply(unique(groups), FUN = function(x) as.numeric(groups == x)/sum(groups == x)))
    fun_connectomes = lapply(fun_connectomes, function(x) {
      connectome = P_community %*% x %*% t(P_community)
      rownames(connectome) = colnames(connectome) = unique(groups)
      connectome
    })

    groups = colnames(str_connectomes[[1]])
    P_community = t(sapply(unique(groups), FUN = function(x) as.numeric(groups == x)/sum(groups == x)))
    str_connectomes = lapply(str_connectomes, function(x) {
      connectome = P_community %*% x %*% t(P_community)
      rownames(connectome) = colnames(connectome) = unique(groups)
      connectome
    })

  }

  Y_fun = sapply(fun_connectomes, from_conn_to_vec) %>% t
  Y_str = sapply(str_connectomes, from_conn_to_vec) %>% t
  colnames(Y_fun) = paste0("fun", 1:ncol(Y_fun))
  colnames(Y_str) = paste0("str", 1:ncol(Y_str))
  Y = cbind(Y_fun, Y_str)

  Ys = list(combined = Y) # , fun = Y_fun, str = Y_str)

  for (scale_type in c("matrix")) { # , "column")) {

    for (fixed_effects in c(T)) { # , F)) {

      for (type in names(Ys)) {

        Y = Ys[[type]]
        Y_mean = colMeans(Y)

        if (fixed_effects) {
          cols = 1:ncol(X)
        } else {
          cols = 1
        }

        residuals = lsfit(X[, cols, drop = FALSE], Y, intercept = FALSE)$residuals
        colnames(residuals) = colnames(Y)

        Y = residuals

        if (type == "combined") {
          Y[, grep("fun", colnames(Y))] = Y[, grep("fun", colnames(Y))] / norm(Y[, grep("fun", colnames(Y))], "F")
          Y[, grep("str", colnames(Y))] = Y[, grep("str", colnames(Y))] / norm(Y[, grep("str", colnames(Y))], "F")
        } else {
          Y = Y / norm(Y, "F")
        }

        if (scale_type == "column") {

          Y = scale(Y)
          Y[, attr(Y, "scaled:scale") == 0] = 0

        }

        D_list = list(diag(nrow = nrow(Y), ncol = nrow(Y)), as.matrix(K_G), as.matrix(K_G > 0) * 1)

        if (type == "combined") {
          V_function = separate_svd_irlba
        } else {
          V_function = mvREHE:::svd_irlba
        }

        # Variance component model estimation (replace code)
        start.time = Sys.time()
        fit = mvREHE::mvREHE_cvDR(Y, D_list = D_list, r_seq = intersect(1:min(ncol(Y), nrow(Y)), 1:4 * 10), V_function = V_function)
        end.time = Sys.time()
        time.taken = end.time - start.time
        print(time.taken)

        # Heritability estimates
        print(sapply(fit$Sigma_r_hat, function(x) sum(diag(x)))/sum(sapply(fit$Sigma_r_hat, function(x) sum(diag(x)))))

        # Extract loadings of a PCA on the genetic covariance structure
        vv_gen = fit$V %*% eigen(fit$Sigma_r_hat[[2]], symmetric = TRUE)$vectors
        vv_gen = vv_gen %*% diag(sign(colMeans(vv_gen)))

        # Extract loadings of a PCA on the common environment covariance structure
        vv_common = fit$V %*% eigen(fit$Sigma_r_hat[[3]], symmetric = TRUE)$vectors
        vv_common = vv_common %*% diag(sign(colMeans(vv_common)))

        # Extract loadings of a PCA on the raw covariance structure
        vv_raw = irlba::irlba(Y)$v
        vv_raw = vv_raw %*% diag(sign(colMeans(vv_raw)))

        fun_indices = grep("fun", colnames(Y))
        str_indices = grep("str", colnames(Y))

        fun_groups = colnames(fun_connectomes[[1]])
        str_groups = colnames(str_connectomes[[1]])

        for (community in c(T, F)) {

        for (pc in 1:2) {

          figure = list()
          figure[[1]] = plot_connectome_vec(Y_mean[fun_indices], "Functional Mean", groups = fun_groups, community = community)
          figure[[2]] = plot_connectome_vec(Y_mean[str_indices], "Structural Mean", groups = str_groups, community = community)
          figure[[5]] = plot_connectome_vec(vv_gen[fun_indices, pc], paste0("Functional Genetic PC ", pc), groups = fun_groups, community = community)
          figure[[6]] = plot_connectome_vec(vv_gen[str_indices, pc], paste0("Structural Genetic PC ", pc), groups = str_groups, community = community)
          figure[[7]] = plot_connectome_vec(vv_common[fun_indices, pc], paste0("Functional Common Environ PC ", pc), groups = fun_groups, community = community)
          figure[[8]] = plot_connectome_vec(vv_common[str_indices, pc], paste0("Structural Common Environ PC ", pc), groups = str_groups, community = community)
          figure[[3]] = plot_connectome_vec(vv_raw[fun_indices, pc], paste0("Functional Raw PC ", pc), groups = fun_groups, community = community)
          figure[[4]] = plot_connectome_vec(vv_raw[str_indices, pc], paste0("Structural Raw PC ", pc), groups = str_groups, community = community)

          figure = figure[!sapply(figure, function(x) is.na(x)[1])]

          pdf(file.path(RESULT_PATH, paste0(type, "_average", average, "_fixed_effects", fixed_effects, "_scale_", scale_type, "_pc_", pc, "_community", community, "_loadings.pdf")), height = 4.3 / 2 * length(figure) / 4, width = 8)
          print(cowplot::plot_grid(plotlist = figure, ncol = 4, byrow = FALSE))
          dev.off()

        }

        }

      }

    }

  }

}
