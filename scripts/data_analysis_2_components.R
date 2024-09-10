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

# Plot rotated PC modes of variation
plot_connectome_vec = function(connectome_vec, title, groups, community = FALSE, breaks = NULL, colors = NULL, legend = FALSE, upper_triangle = FALSE) {
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
    cluster_rows = FALSE,
    cluster_columns = FALSE,
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

separate_svd_irlba = function(Y, r) {

  V1 = mvREHE:::svd_irlba(Y[, grep("fun", colnames(Y))], r)
  V2 = mvREHE:::svd_irlba(Y[, grep("str", colnames(Y))], r)

  V = as.matrix(Matrix::bdiag(V1, V2))
  attr(V, "groups") = c(rep(1, r), rep(2, r))

  V

}

cv_component_ridge_regression = function(Y, D_list, Sigma_hat, outcomes, covariates, estimator, n_lambda = 10, K = 2, folds = NULL) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y)/K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }

  components = 1:length(Sigma_hat)

  lambda_grid = array(NA, dim = c(length(components), n_lambda))
  cv_loss = array(0, dim = c(length(components), length(outcomes), n_lambda))

  for (k in 1:K) {

    fit_train = estimator(Y[-folds[[k]], ], D_list = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]]))
    fit_test = estimator(Y[folds[[k]], ], D_list = lapply(D_list, function(D) D[folds[[k]], folds[[k]]]))

    for (c in 1:length(components)) {

      component = components[c]

      cor_component_hat = cov2cor(Sigma_hat[[component]])
      cor_component_hat[is.na(cor_component_hat)] = 0

      max_eigenvalue = max(eigen(cor_component_hat)$val)
      lambda_grid[component, ] = seq(max_eigenvalue / 200, max_eigenvalue / 5, length.out = 10)

      cor_hat_train = cov2cor(fit_train$Sigma_hat[[component]])
      cor_hat_train[is.na(cor_hat_train)] = 0

      cor_hat_test = cov2cor(fit_test$Sigma_hat[[component]])
      cor_hat_test[is.na(cor_hat_test)] = 0

      for (o in 1:length(outcomes)) {

        outcome = outcomes[o]

        for (l in 1:n_lambda) {

          beta_hat = solve(cor_hat_train[covariates, covariates] + diag(lambda_grid[c, l], length(covariates), length(covariates))) %*% cor_hat_train[covariates, outcome]

          cv_loss[component, outcome, l] = cv_loss[component, outcome, l] - 2 * cor_hat_test[outcome, covariates] %*% beta_hat + t(beta_hat) %*% cor_hat_test[covariates, covariates] %*% beta_hat

        }

      }

    }

  }

  lambda = matrix(NA, length(components), length(outcomes))

  for (c in 1:length(components)) {

    for (o in 1:length(outcomes)) {

      lambda[c, o] = lambda_grid[c, ][which.min(cv_loss[c, o, ])]

    }

  }

  return(lambda)

}

cv_ridge_regression = function(Y, covariates, outcomes, lambda_seq, K = 5, folds = NULL) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y)/K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }

  cv_loss = numeric(length(lambda_seq))

  for (l in 1:length(lambda_seq)) {

    for (k in 1:K) {

      cor_hat_train = cor(Y[-folds[[k]], ])
      cor_hat_train[is.na(cor_hat_train)] = 0
      beta_hat = solve(cor_hat_train[covariates, covariates] + diag(lambda_seq[l], length(covariates), length(covariates))) %*% cor_hat_train[covariates, outcomes]

      cor_hat_test = cor(Y[folds[[k]], ])
      cor_hat_test[is.na(cor_hat_test)] = 0

      cv_loss[l] = cv_loss[l] - 2 * cor_hat_test[outcomes, covariates] %*% beta_hat + t(beta_hat) %*% cor_hat_test[covariates, covariates] %*% beta_hat

    }

  }

  lambda = lambda_seq[which.min(cv_loss)]
  attr(lambda, "cv_loss") = cv_loss

  return(lambda)

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
common_subjects = setdiff(common_subjects, id_subjects[366])
# keep only rows with common_subjects in each data frame

K_G = K_G[id_subjects %in% common_subjects,id_subjects %in% common_subjects] # Kinship matrix
rownames(K_G) = as.character(common_subjects);
colnames(K_G) = as.character(common_subjects);

fun_connectomes = fun_connectomes[id_subjects %in% common_subjects]
str_connectomes = str_connectomes[id_subjects %in% common_subjects]

#################################
##########  Analysis ############
#################################

X = read.csv("data/conf.csv")
rownames(X) = X$subject
X = X[common_subjects, c("Age", "Age.2", "Sex", "FS_IntraCranial_Vol..1.3.", "FS_BrainSeg_Vol..1.3.")]
X = cbind(rep(1, nrow(X)), X)

RESULT_PATH = "data_analysis_2_components"
dir.create(RESULT_PATH, recursive = TRUE)

pre_average = TRUE

groups = colnames(fun_connectomes[[1]])
P_community = t(sapply(unique(groups), FUN = function(x) as.numeric(groups == x)/sum(groups == x)))
fun_connectomes = lapply(fun_connectomes, function(x) {
  connectome = P_community %*% x %*% t(P_community)
  rownames(connectome) = colnames(connectome) = unique(groups)
  connectome[sort(rownames(connectome)), sort(rownames(connectome))]
})

groups = colnames(str_connectomes[[1]])
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

fun_groups = colnames(fun_connectomes[[1]])
str_groups = colnames(str_connectomes[[1]])

connection_names = outer(fun_groups, fun_groups, "paste")
connection_names = connection_names[lower.tri(connection_names, diag = TRUE)]

fixed_effects = TRUE
residuals = lsfit(X, Y, intercept = FALSE)$residuals
colnames(residuals) = colnames(Y)
Y = residuals

scale_type = "column"
Y = scale(Y)
Y[, attr(Y, "scaled:scale") == 0] = 0
# Y[, fun_indices] = Y[, fun_indices] / norm(Y[, fun_indices], "F")
# Y[, str_indices] = Y[, str_indices] / norm(Y[, str_indices], "F")

D_list = list(diag(nrow = nrow(Y), ncol = nrow(Y)), as.matrix(K_G))
genetic_component = 2
unique_env_component = 1

# Variance component model estimation (replace code)
start.time = Sys.time()
# fit = mvREHE::mvREHE_cvDR(Y, D_list = D_list, r_seq = 1:9 * 10, V_function = separate_svd_irlba)
# fit$Sigma_hat = lapply(fit$Sigma_r_hat, function(Sigma) fit$V %*% Sigma %*% t(fit$V))
fit = mvREHE::mvREHE(Y, D_list = D_list)
Sigma_hat = lapply(fit$Sigma_hat, function(Sigma) diag(attr(Y, "scaled:scale")) %*% Sigma %*% diag(attr(Y, "scaled:scale")))
end.time = Sys.time()
time.taken = end.time - start.time
print(time.taken)

saveRDS(fit, file.path(RESULT_PATH, "fit.rds"))
