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

cv_component_ridge_regression = function(Y, D_list, component, covariates, outcomes, lambda_seq, K = 5, folds = NULL, tolerance = 1e-3, max_iter = 100) {

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

      # fit_train = mvREHE::mvREHE_cvDR(Y[-folds[[k]], ], D_list = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]]), r_seq = 1:4 * 10, V_function = separate_svd_irlba, tolerance = tolerance, max_iter = max_iter)
      # fit_train$Sigma_hat = lapply(fit_train$Sigma_r_hat, function(Sigma) fit_train$V %*% Sigma %*% t(fit_train$V))
      fit_train = mvREHE::mvREHE(Y[-folds[[k]], ], D_list = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]]), tolerance = tolerance, max_iter = max_iter)
      cor_hat_train = cov2cor(fit_train$Sigma_hat[[component]])
      beta_hat = solve(cor_hat_train[covariates, covariates] + diag(lambda_seq[l], length(covariates), length(covariates))) %*% cor_hat_train[covariates, outcomes]

      # fit_test = mvREHE::mvREHE_cvDR(Y[folds[[k]], ], D_list = lapply(D_list, function(D) D[folds[[k]], folds[[k]]]),  r_seq = 1:4 * 10, V_function = separate_svd_irlba, tolerance = tolerance, max_iter = max_iter)
      # fit_test$Sigma_hat = lapply(fit_test$Sigma_r_hat, function(Sigma) fit_test$V %*% Sigma %*% t(fit_test$V))
      fit_test = mvREHE::mvREHE(Y[folds[[k]], ], D_list = lapply(D_list, function(D) D[folds[[k]], folds[[k]]]), tolerance = tolerance, max_iter = max_iter)
      cor_hat_test = cov2cor(fit_test$Sigma_hat[[component]])

      cv_loss[l] = cv_loss[l] - 2 * cor_hat_test[outcomes, covariates] %*% beta_hat + t(beta_hat) %*% cor_hat_test[covariates, covariates] %*% beta_hat

    }

  }

  lambda = lambda_seq[which.min(cv_loss)]
  attr(lambda, "cv_loss") = cv_loss

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

X = read.csv("data/conf.csv")
rownames(X) = X$subject
X = X[common_subjects, c("Age", "Age.2", "Sex", "FS_IntraCranial_Vol..1.3.", "FS_BrainSeg_Vol..1.3.")]
X = cbind(rep(1, nrow(X)), X)

RESULT_PATH = "data_analysis_results"
dir.create(RESULT_PATH, recursive = TRUE)

pre_average = TRUE

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
outcome = grep("SMM AUD", connection_names)[1]

fixed_effects = TRUE
residuals = lsfit(X, Y, intercept = FALSE)$residuals
colnames(residuals) = colnames(Y)
Y = residuals

scale_type = "column"
Y = scale(Y)
Y[, attr(Y, "scaled:scale") == 0] = 0
# Y[, fun_indices] = Y[, fun_indices] / norm(Y[, fun_indices], "F")
# Y[, str_indices] = Y[, str_indices] / norm(Y[, str_indices], "F")

D_list = list(diag(nrow = nrow(Y), ncol = nrow(Y)), as.matrix(K_G), as.matrix(K_G > 0) * 1)
genetic_component = 2
common_env_component = 3
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

# Heritability estimate
print(sapply(Sigma_hat, function(x) sum(diag(x)))/sum(sapply(Sigma_hat, function(x) sum(diag(x)))))

# Extract loadings of a PCA on the genetic correlation structure
cor_gen = cov2cor(Sigma_hat[[genetic_component]])
cor_gen[is.na(cor_gen)] = 0
vv_gen = eigen(cor_gen, symmetric = TRUE)$vectors
vv_gen = vv_gen %*% diag(sign(colMeans(vv_gen)))

# Extract loadings of a PCA on the raw covariance structure
cor_Y = cor(Y)
cor_Y[is.na(cor_Y)] = 0
vv_raw = eigen(cor_Y, symmetric = TRUE)$vectors
vv_raw = vv_raw %*% diag(sign(colMeans(vv_raw)))

community = TRUE

figure = list()
# figure[[1]] = plot_connectome_vec(Y_mean[fun_indices], "Functional Mean", groups = fun_groups, community = community)
# figure[[2]] = plot_connectome_vec(Y_mean[str_indices], "Structural Mean", groups = str_groups, community = community)
figure[[4]] = plot_connectome_vec(vv_gen[fun_indices, 1], "Genetic Component PC 1 - Functional", groups = fun_groups, community = community)
figure[[3]] = plot_connectome_vec(vv_gen[str_indices, 1], "Genetic Component PC 1 - Structural", groups = str_groups, community = community)
# figure[[7]] = plot_connectome_vec(vv_common[fun_indices, pc], paste0("Functional Common Environ PC ", pc), groups = fun_groups, community = community)
# figure[[8]] = plot_connectome_vec(vv_common[str_indices, pc], paste0("Structural Common Environ PC ", pc), groups = str_groups, community = community)
figure[[2]] = plot_connectome_vec(vv_raw[fun_indices, 1], "Observed Data PC 1 - Functional", groups = fun_groups, community = community)
figure[[1]] = plot_connectome_vec(vv_raw[str_indices, 1], "Observed Data PC 1 - Structural", groups = str_groups, community = community)
figure[[8]] = plot_connectome_vec(vv_gen[fun_indices, 2], "Genetic Component PC 2 - Functional", groups = fun_groups, community = community)
figure[[7]] = plot_connectome_vec(vv_gen[str_indices, 2], "Genetic Component PC 2 - Structural", groups = str_groups, community = community)
figure[[6]] = plot_connectome_vec(vv_raw[fun_indices, 2], "Observed Data PC 2 - Functional", groups = fun_groups, community = community)
figure[[5]] = plot_connectome_vec(vv_raw[str_indices, 2], "Observed Data PC 2 - Structural", groups = str_groups, community = community)
figure = figure[!sapply(figure, function(x) is.na(x)[1])]

pdf(file.path(RESULT_PATH, "pc_loadings_new.pdf"), height = 5.3, width = 10)
print(cowplot::plot_grid(plotlist = figure, ncol = 4, byrow = TRUE))
dev.off()

beta_components = list()
for (component in 1:3) {
  cor_component_hat = cov2cor(Sigma_hat[[component]])
  cor_component_hat[is.na(cor_component_hat)] = 0
  max_eigenvalue = max(eigen(cor_component_hat)$val)
  lambda = cv_component_ridge_regression(Y, D_list = D_list, component = component, covariates = str_indices, outcomes = outcome, lambda_seq = seq(max_eigenvalue / 200, max_eigenvalue / 5, length.out = 10), K = 2)
  beta_components[[component]] = solve(cor_component_hat[str_indices, str_indices] + diag(lambda, length(str_indices), length(str_indices))) %*% cor_component_hat[str_indices, outcome]
}

cor_raw_hat = cor(Y)
cor_raw_hat[is.na(cor_raw_hat)] = 0
lambda_raw = cv_ridge_regression(Y, covariates = str_indices, outcomes = outcome, lambda_seq = seq(1, 10, by = 0.5), K = 2)
beta_raw = solve(cor_raw_hat[str_indices, str_indices] + diag(lambda_raw, length(str_indices), length(str_indices))) %*% cor_raw_hat[str_indices, 9]

figure = list()
figure[[1]] = plot_connectome_vec(beta_raw, "Coefficients from Observed Data", groups = str_groups, community = community, breaks = c(-0.05, 0, 0.05), colors = c("blue", "white", "red"))
figure[[2]] = plot_connectome_vec(beta_components[[genetic_component]], "Coefficients from Genetic Component", groups = str_groups, community = community, breaks = c(-0.05, 0, 0.05), colors = c("blue", "white", "red"))
figure[[3]] = plot_connectome_vec(beta_components[[common_env_component]], "Coefficients from Common Env Component", groups = str_groups, community = community, breaks = c(-0.05, 0, 0.05), colors = c("blue", "white", "red"))
figure[[4]] = plot_connectome_vec(beta_components[[unique_env_component]], "Coefficients from Unique Env Component", groups = str_groups, community = community, breaks = c(-0.05, 0, 0.05), colors = c("blue", "white", "red"))

pdf(file.path(RESULT_PATH, "regression_coefficients_new.pdf"), height = 2.6, width = 10)
print(cowplot::plot_grid(plotlist = figure, ncol = 4, byrow = TRUE))
dev.off()
