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

library(Rcpp)

cppFunction('
NumericMatrix compute_A_for_Sigma_y_Xbeta2(NumericMatrix beta2, NumericMatrix V) {
  int p = beta2.nrow();
  int rank = beta2.ncol();
  int numCols = V.ncol();
  NumericMatrix result(1 + p * rank, numCols);

  // First row is V[0,]
  for(int j = 0; j < numCols; j++) {
    result(0, j) = V(0, j);
  }

  // Compute A = rbind(c(1, rep(0, p^2)),
  //                   cbind(rep(0, p * rank),
  //                         kronecker(t(beta2), diag(1, p, p)))) %*% V
  for(int r = 0; r < rank; r++) {
    for(int i = 0; i < p; i++) {
      for(int col = 0; col < numCols; col++) {
        double sum = 0.0;
        for(int k = 0; k < p; k++) {
          sum += beta2(k, r) * V(1 + i*p + k, col);
        }
        result(1 + r*p + i, col) = sum;
      }
    }
  }

  return result;
}')

cppFunction('
NumericMatrix compute_A_for_Sigma_y_Xtbeta1(NumericMatrix beta1, NumericMatrix V) {
  int p = beta1.nrow();
  int rank = beta1.ncol();
  int numCols = V.ncol();
  NumericMatrix result(1 + p * rank, numCols);

  // First row is V[0,]
  for(int j = 0; j < numCols; j++) {
    result(0, j) = V(0, j);
  }

  // Compute A = rbind(c(1, rep(0, p^2)),
  //                   cbind(rep(0, p * rank),
  //                         kronecker(diag(1, p, p), t(beta1)))) %*% V
  for(int r = 0; r < rank; r++) {
    for(int i = 0; i < p; i++) {
      for(int col = 0; col < numCols; col++) {
        double sum = 0.0;
        for(int k = 0; k < p; k++) {
          sum += beta1(k, r) * V(1 + i + k*p, col);
        }
        result(1 + r*p + i, col) = sum;
      }
    }
  }

  return result;
}')

cppFunction('
NumericMatrix compute_A_for_objective(NumericMatrix beta1,
                                      NumericMatrix beta2,
                                      NumericMatrix V) {
  int p = beta1.nrow();
  int rank = beta1.ncol();
  int numCols = V.ncol();
  NumericMatrix result(2, numCols);

  // First row of result is just V[0,]
  for(int j = 0; j < numCols; j++) {
    result(0, j) = V(0, j);
  }

  // Second row is the result of c(0, c(tcrossprod(beta1, beta2))) %*% V
  for(int col = 0; col < numCols; col++) {
    double sum = 0.0;

    // Skip the first element of V (corresponds to the 0 in the second row)
    for(int i = 0; i < p; i++) {
      for(int j = 0; j < p; j++) {
        // Calculate element-wise contribution from beta1 and beta2
        double beta_prod = 0.0;
        for(int r = 0; r < rank; r++) {
          beta_prod += beta1(i, r) * beta2(j, r);
        }

        // Add contribution to the sum
        sum += beta_prod * V(1 + i*p + j, col);
      }
    }

    result(1, col) = sum;
  }

  return result;
}')

matrix_regression_from_cov = function(Sigma_r_hat, V, rank, lambda, max_iter = 1000, tolerance = 1e-6, n_init = 5) {

  if (!is.null(V)) {
    p = sqrt(nrow(V) - 1)
  } else {
    p = sqrt(nrow(Sigma_r_hat) - 1)
  }
  best_objective = Inf
  best_solution = NULL

  for (init in 1:n_init) {

    beta1 = matrix(rnorm(p * rank, sd = 1), p, rank)
    beta2 = beta1

    objective = numeric(2 * max_iter)
    difference = numeric(max_iter)

    for (iter in 1:max_iter) {

      beta1_old = beta1
      beta2_old = beta2

      if (!is.null(V)) {
        # A = compute_A_for_Sigma_y_Xbeta2(beta2, V)
        A = rbind(c(1, rep(0, p^2)), cbind(rep(0, p * rank), kronecker(t(beta2), diag(1, p, p)))) %*% V
      } else {
        A = rbind(c(1, rep(0, p^2)), cbind(rep(0, p * rank), kronecker(t(beta2), diag(1, p, p))))
      }
      Sigma = A %*% Sigma_r_hat %*% t(A)
      beta1[] = solve(Sigma[-1, -1] + diag(lambda, ncol(Sigma) - 1, ncol(Sigma) - 1), Sigma[1, -1])

      if (!is.null(V)) {
        # A = compute_A_for_objective(beta1, beta2, V)
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2)))) %*% V
      } else {
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2))))
      }
      Sigma = A %*% Sigma_r_hat %*% t(A)
      objective[2 * iter - 1] = Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2] + lambda * (sum(beta1^2) + sum(beta2^2))

      if (!is.null(V)) {
        # A = compute_A_for_Sigma_y_Xtbeta1(beta1, V)
        A = rbind(c(1, rep(0, p^2)), cbind(rep(0, p * rank), kronecker(diag(1, p, p), t(beta1)))) %*% V
      } else {
        A = rbind(c(1, rep(0, p^2)), cbind(rep(0, p * rank), kronecker(diag(1, p, p), t(beta1))))
      }
      Sigma = A %*% Sigma_r_hat %*% t(A)
      beta2[] = matrix(solve(Sigma[-1, -1] + diag(lambda, ncol(Sigma) - 1, ncol(Sigma) - 1), Sigma[1, -1]), ncol = rank, byrow = TRUE)

      if (!is.null(V)) {
        # A = compute_A_for_objective(beta1, beta2, V)
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2)))) %*% V
      } else {
        A = rbind(c(1, rep(0, p^2)), c(0, c(tcrossprod(beta1, beta2))))
      }
      Sigma = A %*% Sigma_r_hat %*% t(A)
      objective[2 * iter] = Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2] + lambda * (sum(beta1^2) + sum(beta2^2))
      difference[iter] = norm(beta1 %*% t(beta2) - beta1_old %*% t(beta2_old), "F") / (norm(beta1 %*% t(beta2), "F") + 1e-12)

      if (iter > 5 && difference[iter] < tolerance) {
        break
      }
    }

    final_objective = objective[max(which(objective != 0))]
    if (final_objective < best_objective) {
      print(c(final_objective, best_objective))
      best_objective = final_objective
      best_solution = list(
        beta1 = beta1,
        beta2 = beta2,
        difference = difference[difference != 0],
        objective = objective[objective != 0],
        init_number = init,
        final_objective = final_objective
      )
    }
  }

  return(best_solution)
}

latent_matrix_regression = function(Y, D_list, component, outcome, covariates, rank, lambda, estimator, ...) {

  fit = estimator(Y, D_list)

  Sigma_r_hat = fit$Sigma_hat[[component]]

  covariates = vech_to_vec_indices(covariates)

  if (!is.null(fit$V)) {
    V <- fit$V[c(outcome, covariates), ]
  } else {
    V = NULL
    Sigma_r_hat = Sigma_r_hat[c(outcome, covariates), c(outcome, covariates)]
  }

  matrix_regression_from_cov(Sigma_r_hat, V, rank, lambda, ...)

}


matrix_regression = function(Y, outcome, covariates, rank, lambda, ...) {

  Sigma_r_hat = cov(Y)

  covariates = vech_to_vec_indices(covariates)

  V = NULL
  Sigma_r_hat = Sigma_r_hat[c(outcome, covariates), c(outcome, covariates)]

  matrix_regression_from_cov(Sigma_r_hat, V, rank, lambda, ...)

}

vech_to_vec_indices <- function(vech_indices) {

  q <- length(vech_indices)
  p <- (sqrt(8*q + 1) - 1)/2

  if(p != round(p)) {
    stop("Input vector length is not valid for a symmetric matrix vech representation")
  }

  mat = matrix(0, p, p)

  mat[lower.tri(mat, diag = TRUE)] = vech_indices

  c(mat + t(mat) - diag(diag(mat), p, p))

}

cv_latent_matrix_regression <- function(Y, D_list, outcomes, covariates, rank_seq, lambda_seq, estimator,
                                        K = 2, folds = NULL, cores = parallel::detectCores() - 1, ...) {

  require(parallel)

  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y) / K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }

  covariates = vech_to_vec_indices(covariates)

  cv_r2 = array(0, dim = c(length(D_list), length(outcomes), length(rank_seq), length(lambda_seq)))

  for (k in 1:K) {

    fit_train = estimator(Y[-folds[[k]], ], D_list = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]]), return_full = FALSE)
    fit_test = estimator(Y[folds[[k]], ], D_list = lapply(D_list, function(D) D[folds[[k]], folds[[k]]]), return_full = FALSE)

    for (c in 1:length(D_list)) {

      # Create full grid of outcome × rank × lambda
      grid <- expand.grid(
        o_idx = seq_along(outcomes),
        r_idx = seq_along(rank_seq),
        l_idx = seq_along(lambda_seq)
      )

      results_grid <- mclapply(seq_len(nrow(grid)), function(idx) {
        row <- grid[idx, ]
        o_idx <- row$o_idx
        r_idx <- row$r_idx
        l_idx <- row$l_idx

        o <- outcomes[o_idx]
        r <- rank_seq[r_idx]
        l <- lambda_seq[l_idx]

        # Training
        Sigma_r_hat_train <- fit_train$Sigma_r_hat[[c]]
        if (!is.null(fit_train$V)) {
          V_train <- fit_train$V[c(o, covariates), ]
        } else {
          V_train = NULL
          Sigma_r_hat_train = Sigma_r_hat_train[c(o, covariates), c(o, covariates)]
        }
        fit <- matrix_regression_from_cov(Sigma_r_hat_train, V_train, r, l, ...)

        # Testing
        Sigma_r_hat_test <- fit_test$Sigma_r_hat[[c]]
        if (!is.null(fit_test$V)) {
          V_test <- fit_test$V[c(o, covariates), ]
        } else {
          V_test = NULL
          Sigma_r_hat_test = Sigma_r_hat_test[c(o, covariates), c(o, covariates)]
        }

        A <- rbind(
          c(1, rep(0, length(covariates))),
          c(0, c(tcrossprod(fit$beta1, fit$beta2)))
        )

        if (!is.null(V_test)) {
          A = A %*% V_test
        }

        Sigma <- A %*% Sigma_r_hat_test %*% t(A)
        r2 <- 1 - (Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2]) / Sigma[1, 1]

        list(o_idx = o_idx, r_idx = r_idx, l_idx = l_idx, r2 = r2 / K)
      }, mc.cores = cores)

      # Accumulate results
      for (res in results_grid) {
        cv_r2[c, res$o_idx, res$r_idx, res$l_idx] <-
          cv_r2[c, res$o_idx, res$r_idx, res$l_idx] + res$r2
      }
    }
  }

  # Select best rank/lambda based on max CV R²
  lambda = matrix(NA, length(D_list), length(outcomes))
  rank = matrix(NA, length(D_list), length(outcomes))

  for (c in 1:length(D_list)) {
    for (o in 1:length(outcomes)) {
      idx = which(cv_r2[c, o, , , drop = FALSE] == max(cv_r2[c, o, , ]), arr.ind = TRUE)[1, ]
      rank[c, o] = rank_seq[idx[3]]
      lambda[c, o] = lambda_seq[idx[4]]
    }
  }

  return(list(
    lambda = lambda,
    rank = rank,
    cv_r2 = pmax(apply(cv_r2, c(1, 2), max), 0),
    cv_r2_full = cv_r2
  ))

}

cv_matrix_regression <- function(Y, outcomes, covariates, rank_seq, lambda_seq,
                                 K = 2, folds = NULL, cores = parallel::detectCores() - 1, ...) {

  require(parallel)

  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y) / K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }

  covariates = vech_to_vec_indices(covariates)

  cv_r2 = array(0, dim = c(length(outcomes), length(rank_seq), length(lambda_seq)))

  for (k in 1:K) {

    cov_train = cov(Y[-folds[[k]], ])
    cov_test = cov(Y[folds[[k]], ])

      # Create full grid of outcome × rank × lambda
      grid <- expand.grid(
        o_idx = seq_along(outcomes),
        r_idx = seq_along(rank_seq),
        l_idx = seq_along(lambda_seq)
      )

      results_grid <- mclapply(seq_len(nrow(grid)), function(idx) {
        row <- grid[idx, ]
        o_idx <- row$o_idx
        r_idx <- row$r_idx
        l_idx <- row$l_idx

        o <- outcomes[o_idx]
        r <- rank_seq[r_idx]
        l <- lambda_seq[l_idx]

        # Training

        fit <- matrix_regression_from_cov(cov_train[c(o, covariates), c(o, covariates)], NULL, r, l, ...)

        # Testing

        A <- rbind(
          c(1, rep(0, length(covariates))),
          c(0, c(tcrossprod(fit$beta1, fit$beta2)))
        )

        Sigma <- A %*% cov_test[c(o, covariates), c(o, covariates)] %*% t(A)
        r2 <- 1 - (Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2]) / Sigma[1, 1]

        list(o_idx = o_idx, r_idx = r_idx, l_idx = l_idx, r2 = r2 / K)
      }, mc.cores = cores)

      # Accumulate results
      for (res in results_grid) {
        cv_r2[res$o_idx, res$r_idx, res$l_idx] <-
          cv_r2[res$o_idx, res$r_idx, res$l_idx] + res$r2
      }

  }

  # Select best rank/lambda based on max CV R²
  lambda = numeric(length(outcomes))
  rank = numeric(length(outcomes))

    for (o in 1:length(outcomes)) {
      idx = which(cv_r2[o, , , drop = FALSE] == max(cv_r2[o, , ]), arr.ind = TRUE)[1, ]
      rank[o] = rank_seq[idx[2]]
      lambda[o] = lambda_seq[idx[3]]
    }

  return(list(
    lambda = lambda,
    rank = rank,
    cv_r2 = pmax(apply(cv_r2, 1, max), 0),
    cv_r2_full = cv_r2
  ))

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

# REORDER THESE LIKE IN SIMULATION!!!!
# 10 x 10 simulations, 68 x 68 application
# 10 x 10 use full covariance
# vary r between 1 and 10 for application
# try and visualize using python package

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

#################################
##########  Analysis ############
#################################

X = read.csv("data/conf.csv")
rownames(X) = X$subject
X = X[common_subjects[order], c("Age", "Age.2", "Sex", "FS_IntraCranial_Vol..1.3.", "FS_BrainSeg_Vol..1.3.")]
X = cbind(rep(1, nrow(X)), X)

RESULT_PATH = "data_analysis_3_components_68ROI_cov"
dir.create(RESULT_PATH, recursive = TRUE)

groups = read.csv("~/Documents/Downloads/communities_aparc.csv")
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
# Y[, fun_indices] = Y[, fun_indices] / norm(Y[, fun_indices], "F")
# Y[, str_indices] = Y[, str_indices] / norm(Y[, str_indices], "F")

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

latent_matrix_regression_fit = cv_latent_matrix_regression(Y, D_list, fun_indices, str_indices, 1:5, 10^seq(0, -2, length.out = 10), mvREHE, 2, n_init = 1)
raw_matrix_regression_fit = cv_matrix_regression(Y, fun_indices, str_indices, 1:5, 10^seq(0, -2, length.out = 10), 2, n_init = 1)

r2_components = lapply(1:length(D_list), function(i) latent_matrix_regression_fit$cv_r2[i, ])
r2_raw = raw_matrix_regression_fit$cv_r2

community = FALSE
figure = list()
figure[[1]] = plot_connectome_vec(r2_raw, bquote("Observed"~R^2~"\n(Min" == .(round(min(r2_raw), 2))~", Max" == .(round(max(r2_raw), 2))~")"), groups = fun_groups, community = community, breaks = c(0, 0.5), colors = c("white", "red"), upper_triangle = TRUE)
figure[[2]] = plot_connectome_vec(r2_components[[genetic_component]], bquote("Genetic"~R^2~"\n(Min" == .(round(min(r2_components[[genetic_component]]), 2))~", Max" == .(round(max(r2_components[[genetic_component]]), 2))~")"), groups = fun_groups, community = community, breaks = c(0, 0.5), colors = c("white", "red"), upper_triangle = TRUE)
figure[[3]] = plot_connectome_vec(r2_components[[common_env_component]], bquote("Common Env"~R^2~"\n(Min" == .(round(min(r2_components[[common_env_component]]), 2))~", Max" == .(round(max(r2_components[[common_env_component]]), 2))~")"), groups = fun_groups, community = community, breaks = c(0, 0.5), colors = c("white", "red"), upper_triangle = TRUE)
figure[[4]] = plot_connectome_vec(r2_components[[unique_env_component]], bquote("Unique Env"~R^2~"\n(Min" == .(round(min(r2_components[[unique_env_component]]), 2))~", Max" == .(round(max(r2_components[[unique_env_component]]), 2))~")"), groups = fun_groups, community = community, breaks = c(0, 0.5), colors = c("white", "red"), upper_triangle = TRUE)

pdf(file.path(RESULT_PATH, "r2_matrix.pdf"), height = 2.9, width = 10)
patchwork::wrap_plots(figure, ncol = 4, byrow = TRUE) +
  patchwork::plot_annotation(tag_levels = list(c("a", "b", "c", "d")))
dev.off()

outcome = fun_indices[which.max(latent_matrix_regression_fit$cv_r2[2, ])]
covariates = str_indices
latent_genetic_matrix_regression_fit_max = latent_matrix_regression(Y, D_list, 2, outcome, covariates, latent_matrix_regression_fit$rank[2, outcome], latent_matrix_regression_fit$lambda[2, outcome], function(Y, D_list) mvREHE(Y, D_list, return_full = FALSE))
latent_unique_environment_matrix_regression_fit_max = latent_matrix_regression(Y, D_list, 1, outcome, covariates, latent_matrix_regression_fit$rank[1, outcome], latent_matrix_regression_fit$lambda[1, outcome], function(Y, D_list) mvREHE(Y, D_list, return_full = FALSE))
latent_common_environment_matrix_regression_fit_max = latent_matrix_regression(Y, D_list, 3, outcome, covariates, latent_matrix_regression_fit$rank[3, outcome], latent_matrix_regression_fit$lambda[3, outcome], function(Y, D_list) mvREHE(Y, D_list, return_full = FALSE))
raw_matrix_regression_fit_max = matrix_regression(Y, outcome, covariates, raw_matrix_regression_fit$rank[outcome], raw_matrix_regression_fit$lambda[outcome])


get_vec_beta = function(beta1, beta2) {
  beta = tcrossprod(beta1, beta2) + tcrossprod(beta2, beta1) - diag(diag(tcrossprod(beta1, beta2)))
  beta[lower.tri(beta, diag = TRUE)]
}

sqrt_matrix = function(A) {
  eig = eigen(A)
  eig$vec %*% diag(sqrt(pmax(eig$val, 0))) %*% t(eig$vec)
}

cond_num = 100
covariates = str_indices

Sigma_hat = fit$Sigma_hat
for (k in 1:length(Sigma_hat)) {
  eig = eigen(Sigma_hat[[k]][covariates, covariates])
  diag(Sigma_hat[[k]])[covariates] = diag(Sigma_hat[[k]])[covariates] + eig$val[1] / (cond_num - 1)
  attr(Sigma_hat[[k]], "sqrt") = sqrt_matrix(Sigma_hat[[k]])
}

beta_0 = get_vec_beta(latent_unique_environment_matrix_regression_fit_max$beta1, latent_unique_environment_matrix_regression_fit_max$beta2)
beta_1 = get_vec_beta(latent_genetic_matrix_regression_fit_max$beta1, latent_genetic_matrix_regression_fit_max$beta2)
beta_2 = get_vec_beta(latent_common_environment_matrix_regression_fit_max$beta1, latent_common_environment_matrix_regression_fit_max$beta2)

beta = list(beta_0, beta_1, beta_2)

r2 = c(0.01, 0.5, 0.1)

Sigma_hat = lapply(1:3, function(i) {
  A = rbind(beta[[i]], diag(1, length(beta[[i]])))
  cov = A %*% Sigma_hat[[i]][covariates, covariates] %*% t(A)
  sigma2 = (1 - r2[i]) * cov[1, 1] / r2[i]
  cov[1, 1] = cov[1, 1] + sigma2
  cov
})

# outcome = 1
# covariates = -1
# beta = beta[[2]]
# k = 2
# 1 - (Sigma_hat[[k]][outcome, outcome] - 2 * Sigma_hat[[k]][outcome, covariates] %*% beta + t(beta) %*% Sigma_hat[[k]][covariates, covariates] %*% beta) / Sigma_hat[[k]][outcome, outcome]

fit$Sigma_hat_sim = Sigma_hat
fit$beta = beta

saveRDS(fit, file.path(RESULT_PATH, "fit.rds"))

pdf(file.path(RESULT_PATH, "r2_matrix_hist.pdf"), height = 3, width = 8)
par(mfrow = c(1, 3))
for (i in c(2, 3, 1)) {
  hist(latent_matrix_regression_fit$cv_r2[i, ])
}
dev.off()


