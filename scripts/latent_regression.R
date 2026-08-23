library(igraph)
library(parallel)
Rcpp::sourceCpp("scripts/lasso_from_cov.cpp")
Rcpp::sourceCpp("scripts/matrix_regression_from_cov.cpp")

family_ids_from_kinship = function(kinship) {
  g = graph_from_adjacency_matrix(kinship != 0, mode = "undirected", diag = FALSE)
  components(g)$membership
}
make_family_folds = function(fam, K) {
  families = sample(unique(fam))
  stopifnot(length(families) >= K)
  fold_of_family = rep(seq_len(K), length.out = length(families))[match(fam, families)]
  split(seq_along(fam), fold_of_family)
}
weights_cov_to_cor = function(Sigma) {
  variances = diag(Sigma)
  variances[variances < 1e-10] = 1
  1 / sqrt(variances)
}
weight_Sigma = function(Sigma, s) {
  sweep(sweep(Sigma, 1, s, "*"), 2, s, "*")
}
sigma_diag = function(Sigma, V) {
  if (!is.null(V)) rowSums((V %*% Sigma) * V) else diag(Sigma)
}
ridge_regression_from_cov = function(Sigma, outcomes, covariates, tuning_grids, cores = 1, V = NULL) {
  lambda_seq = tuning_grids[[1]]$lambda
  vars = c(covariates, outcomes)
  weights = weights_cov_to_cor(Sigma[vars, vars])
  corr = weight_Sigma(Sigma[vars, vars], weights)
  p = length(covariates); m = length(outcomes)
  corr_X = corr[1:p, 1:p]
  corr_XY = corr[1:p, (p + 1):(p + m), drop = FALSE]
  s_cov = weights[1:p]
  s_out = weights[(p + 1):(p + m)]

  eig = eigen(corr_X, symmetric = TRUE)
  eig_V = eig$vectors; d = eig$values
  VtXY = crossprod(eig_V, corr_XY)

  beta = array(0, dim = c(p, m, length(lambda_seq)))
  for (i in seq_along(lambda_seq)) {
    B = eig_V %*% (VtXY / (d + lambda_seq[i]))
    beta[, , i] = sweep(B, 2, s_out, "/") * s_cov
  }
  beta
}
ridge_tuning_grid_fn = function(fit, outcomes, covariates, n_lambda = 10) {
  lapply(seq_along(fit$Sigma_hat), function(comp) {
    Sigma = fit$Sigma_hat[[comp]][covariates, covariates]
    weights = weights_cov_to_cor(Sigma)
    corr = weight_Sigma(Sigma, weights)
    max_ev = max(RSpectra::eigs_sym(corr, k = 1)$val)
    grid = data.frame(lambda = 10 ^ seq(log10(max_ev), log10(max_ev / 10000),
                                        length.out = n_lambda))
    rep(list(grid), length(outcomes))
  })
}
cv_latent_regression = function(Y, D_list, family_ids, fit, outcomes, covariates,
                                tuning_grid_fn, estimator_fn, regression_fn,
                                K = 5, cores = 8) {
  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)
  n = nrow(Y)
  folds = make_family_folds(family_ids, K)
  tuning_grids = tuning_grid_fn(fit, outcomes, covariates)
  n_tune = nrow(tuning_grids[[1]][[1]])

  predictions = lapply(1:length(D_list), function(comp)
    array(NA_real_, dim = c(n, length(outcomes), n_tune)))
  for (k in 1:K) {
    print(k)
    train = setdiff(1:n, folds[[k]])
    test = folds[[k]]
    fit_train = estimator_fn(Y[train, ], D_list = lapply(D_list, function(D)
      D[train, train]))
    Xtest = Y[test, covariates, drop = FALSE]
    for (component in 1:length(D_list)) {
      beta = regression_fn(fit_train$Sigma_hat[[component]], outcomes, covariates,
                           tuning_grids[[component]], cores = cores, V = fit_train$V)
      for (l in 1:n_tune)
        predictions[[component]][test, , l] = Xtest %*% beta[, , l]
    }
  }

  eval_grid = expand.grid(component = 1:length(D_list), l = 1:n_tune)
  degenerate = apply(Y[, outcomes, drop = FALSE], 2, sd) < 1e-8
  m = length(outcomes)
  eval_out = mclapply(1:nrow(eval_grid), function(i) {
    component = eval_grid$component[i]; l = eval_grid$l[i]
    resid = Y[, outcomes] - predictions[[component]][, , l]
    r2_fit = estimator_fn(cbind(Y, resid), D_list = D_list)
    d = sigma_diag(r2_fit$Sigma_hat[[component]], r2_fit$V)
    var_y = d[outcomes]
    var_r = d[(ncol(Y) + 1):(ncol(Y) + m)]
    r2 = 1 - var_r / var_y
    r2[degenerate | !is.finite(var_y) | var_y < 1e-10 | !is.finite(var_r)] = NA_real_
    r2
  }, mc.cores = cores)

  cv_r2 = array(NA_real_, dim = c(length(D_list), length(outcomes), n_tune))
  for (i in 1:nrow(eval_grid)) cv_r2[eval_grid$component[i], , eval_grid$l[i]] = eval_out[[i]]

  ## Selection with null-model fallback.
  ##   selected == 0L  -> null model (beta = 0): either no finite CV values, or
  ##                      the peak falls at the largest lambda (index 1), meaning
  ##                      more regularization is always better -> null model optimal.
  ##   selected >= 1L  -> index of the chosen tuning value.
  selected = matrix(NA_integer_, length(D_list), length(outcomes))
    for (component in 1:length(D_list)) {
    for (o in 1:length(outcomes)) {
      v = cv_r2[component, o, ]
      print(v)
      if (!any(is.finite(v))) {
        selected[component, o] = 0L
      } else {
        best = which.max(replace(v, !is.finite(v), -Inf))
        selected[component, o] = if (max(v, na.rm = TRUE) <= 0) 0L else best
      }
    }
  }

  list(cv_r2 = cv_r2, selected = selected, tuning_grids = tuning_grids)
}
nested_cv_latent_regression = function(Y, D_list, family_ids, fit, outcomes, covariates,
                                       tuning_grid_fn, estimator_fn, regression_fn,
                                       K_outer = 10, K_inner = 5, cores = 8,
                                       checkpoint_path = NULL) {
  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)
  n = nrow(Y)
  folds = make_family_folds(family_ids, K_outer)

  predictions = lapply(1:length(D_list), function(comp)
    matrix(NA_real_, n, length(outcomes)))
  selected = array(NA_integer_, dim = c(K_outer, length(D_list), length(outcomes)))

  if (!is.null(checkpoint_path))
    dir.create(checkpoint_path, recursive = TRUE, showWarnings = FALSE)

  for (k in 1:K_outer) {
    print(k)

    ckpt_file = if (!is.null(checkpoint_path))
      file.path(checkpoint_path, sprintf("fold_%d.rds", k)) else NULL

    if (!is.null(ckpt_file) && file.exists(ckpt_file)) {
      cat(sprintf("Resuming from checkpoint: outer fold %d\n", k))
      ckpt = readRDS(ckpt_file)
      selected[k, , ] = ckpt$selected_k
      for (component in seq_along(D_list))
        predictions[[component]][ckpt$test_indices, ] = ckpt$predictions_k[[component]]
      next
    }

    train = setdiff(1:n, folds[[k]])

    inner = cv_latent_regression(
      Y[train, ], lapply(D_list, function(D) D[train, train]), family_ids[train],
      fit, outcomes, covariates, tuning_grid_fn, estimator_fn, regression_fn,
      K = K_inner, cores = cores)
    selected[k, , ] = inner$selected

    fit_train = estimator_fn(Y[train, ], D_list = lapply(D_list, function(D)
      D[train, train]))
    Xtest = Y[folds[[k]], covariates, drop = FALSE]
    for (component in 1:length(D_list)) {
      beta = regression_fn(fit_train$Sigma_hat[[component]], outcomes, covariates,
                           inner$tuning_grids[[component]], cores = cores, V = fit_train$V)
      for (o in seq_along(outcomes)) {
        sel = inner$selected[component, o]
        ## sel == 0L is the null model: beta = 0, so the prediction is 0 and the
        ## residual equals the outcome (R^2 = 0).
        predictions[[component]][folds[[k]], o] =
          if (sel == 0L) 0 else Xtest %*% beta[, o, sel]
      }
    }

    if (!is.null(ckpt_file))
      saveRDS(list(
        selected_k    = inner$selected,
        predictions_k = lapply(predictions, function(p) p[folds[[k]], , drop = FALSE]),
        test_indices  = folds[[k]]
      ), ckpt_file)
  }

  cv_r2 = matrix(NA_real_, length(D_list), length(outcomes))
  var_outcome = matrix(NA_real_, length(D_list), length(outcomes))
  var_resid = matrix(NA_real_, length(D_list), length(outcomes))
  degenerate = apply(Y[, outcomes, drop = FALSE], 2, sd) < 1e-8
  m = length(outcomes)
  for (component in 1:length(D_list)) {
    resid = Y[, outcomes] - predictions[[component]]
    r2_fit = estimator_fn(cbind(Y, resid), D_list = D_list)
    d = sigma_diag(r2_fit$Sigma_hat[[component]], r2_fit$V)
    var_y = d[outcomes]
    var_r = d[(ncol(Y) + 1):(ncol(Y) + m)]
    var_outcome[component, ] = var_y
    var_resid[component, ] = var_r
    r2 = 1 - var_r / var_y
    r2[degenerate | !is.finite(var_y) | var_y < 1e-10 | !is.finite(var_r)] = NA_real_
    cv_r2[component, ] = r2
  }

  ## Fraction of outer folds in which each (component, outcome) fell back to the
  ## null model (selected == 0L). High for the common-environment component is
  ## expected: many connections have no positive inner C signal.
  prop_null_model = apply(selected == 0L, c(2, 3), mean)

  list(cv_r2 = cv_r2, var_outcome = var_outcome, var_resid = var_resid,
       prop_null_model = prop_null_model,
       selected = selected, predictions = predictions)
}
lasso_tuning_grid_fn = function(fit, outcomes, covariates, n_lambda = 10) {
  lapply(seq_along(fit$Sigma_hat), function(comp) {
    Sig = fit$Sigma_hat[[comp]]
    vars = c(covariates, outcomes)
    weights = weights_cov_to_cor(Sig[vars, vars])
    s_cov = weights[seq_along(covariates)]
    s_out = weights[seq_along(outcomes) + length(covariates)]
    lapply(seq_along(outcomes), function(o) {
      outcome = outcomes[o]
      lambda_max = if (Sig[outcome, outcome] > 1e-10)
        max(abs(Sig[covariates, outcome] * s_cov * s_out[o]))
      else 0
      if (lambda_max > 0)
        data.frame(lambda = exp(seq(log(lambda_max), log(lambda_max * 1e-4),
                                   length.out = n_lambda)))
      else
        data.frame(lambda = rep(NA_real_, n_lambda))
    })
  })
}
lasso_regression_from_cov = function(Sigma, outcomes, covariates, tuning_grids, cores = 1, V = NULL) {
  vars = c(covariates, outcomes)
  weights = weights_cov_to_cor(Sigma[vars, vars])
  corr = weight_Sigma(Sigma[vars, vars], weights)
  p = length(covariates); m = length(outcomes)
  s_cov = weights[1:p]
  s_out = weights[(p + 1):(p + m)]

  n_lambda = nrow(tuning_grids[[1]])
  covariates_nonzero = which(diag(Sigma[covariates, covariates]) > 1e-8)

  beta_list = mclapply(seq_len(m), function(o) {
    lambda_seq = tuning_grids[[o]]$lambda
    if (any(is.na(lambda_seq)))
      return(matrix(0, p, n_lambda))
    idx = c(p + o, covariates_nonzero)
    b_std = matrix(0, p, length(lambda_seq))
    b_std[covariates_nonzero, ] = lasso_from_cov(corr[idx, idx], lambda_seq)
    b_std * (s_cov / s_out[o])
  }, mc.cores = cores)

  beta = array(0, dim = c(p, m, n_lambda))
  for (o in seq_len(m)) beta[, o, ] = beta_list[[o]]
  beta
}
vech_to_vec_indices = function(vech_indices) {
  q = length(vech_indices)
  p = (sqrt(8 * q + 1) - 1) / 2
  mat = matrix(0, p, p)
  mat[lower.tri(mat, diag = TRUE)] = vech_indices
  c(mat + t(mat) - diag(diag(mat), p, p))
}
matrix_tuning_grid_fn = function(fit, outcomes, covariates, rank_seq = 1:5, lambda_seq = NULL) {
  grid = if (is.null(lambda_seq)) data.frame(rank = rank_seq)
         else expand.grid(rank = rank_seq, lambda = lambda_seq)
  lapply(seq_along(fit$Sigma_hat), function(comp) {
    lapply(seq_along(outcomes), function(o) grid)
  })
}
low_rank_regression_from_cov = function(Sigma, outcomes, covariates, tuning_grids,
                                        V = NULL, lambda = 1e-6, cores = 1,
                                        return_factors = FALSE) {
  full_cov = covariates[vech_to_vec_indices(seq_along(covariates))]
  p_ltr = length(covariates)
  p_str = round((sqrt(8 * p_ltr + 1) - 1) / 2)
  m = length(outcomes)
  grid = tuning_grids[[1]]
  rank_seq = grid$rank
  lambda_seq = if (!is.null(grid$lambda)) grid$lambda else rep(lambda, length(rank_seq))
  n_tune = length(rank_seq)

  # mclapply over outcomes; V_sub and B_sym eigdecomp computed once per unique lambda,
  # top-r eigenvectors used as init for each (rank, lambda) (avoiding redundant V*Sigma_hat*V' calls)
  beta_list = mclapply(seq_len(m), function(o) {
    outcome = outcomes[o]
    idx = c(outcome, full_cov)
    if (!is.null(V)) {
      V_sub = V[idx, ]
      Sigma_sub = Sigma
      Sigma_eff = V_sub %*% Sigma_sub %*% t(V_sub)
    } else {
      V_sub = NULL
      Sigma_sub = Sigma[idx, idx]
      Sigma_eff = Sigma_sub
    }
    # Unconstrained OLS B and eigdecomp — depends on lambda only, so computed
    # once per unique lambda and shared across ranks.
    p2 = p_str^2
    Sigma_XX = Sigma_eff[2:(p2 + 1), 2:(p2 + 1)]
    Sigma_XY = Sigma_eff[2:(p2 + 1), 1]
    unique_lambdas = unique(lambda_seq)
    init_by_lambda = lapply(unique_lambdas, function(lam) {
      Sxx = Sigma_XX
      diag(Sxx) = diag(Sxx) + lam
      B_mat = matrix(solve(Sxx, Sigma_XY), p_str, p_str)
      B_sym_mat = 0.5 * (B_mat + t(B_mat))
      eig = eigen(B_sym_mat, symmetric = TRUE)
      ord = order(abs(eig$values), decreasing = TRUE)
      list(ev = eig$values[ord], evec = eig$vectors[, ord, drop = FALSE])
    })
    names(init_by_lambda) = as.character(unique_lambdas)

    beta_o = matrix(0, p_ltr, n_tune)
    factors_o = vector("list", n_tune)
    for (t_idx in seq_len(n_tune)) {
      r = min(rank_seq[t_idx], p_str)
      lam = lambda_seq[t_idx]
      init = init_by_lambda[[as.character(lam)]]
      ev = init$ev; evec = init$evec
      sqrt_ev = sqrt(abs(ev[seq_len(r)]))
      b1_init = t(t(evec[, seq_len(r), drop = FALSE]) * sqrt_ev)
      b2_init = t(t(evec[, seq_len(r), drop = FALSE]) * (sqrt_ev * sign(ev[seq_len(r)])))
      fit = matrix_regression_from_cov(Sigma_sub, V_sub, rank_seq[t_idx], lam,
                                       tolerance = 1e-2,
                                       beta1_init = b1_init,
                                       beta2_init = b2_init)
      B = fit$beta1 %*% t(fit$beta2)
      Bsym = B + t(B) - diag(diag(B))
      beta_o[, t_idx] = Bsym[lower.tri(Bsym, diag = TRUE)]
      factors_o[[t_idx]] = list(beta1 = fit$beta1, beta2 = fit$beta2)
    }
    list(beta = beta_o, factors = factors_o)
  }, mc.cores = cores)

  beta = array(0, dim = c(p_ltr, m, n_tune))
  for (o in seq_len(m)) beta[, o, ] = beta_list[[o]]$beta
  if (!return_factors) return(beta)
  factors = lapply(seq_len(m), function(o) beta_list[[o]]$factors)
  list(beta = beta, factors = factors)
}
raw_tuning_grid_fn = function(Sigma, outcomes, covariates, n_lambda = 10) {
  block = Sigma[covariates, covariates]
  weights = weights_cov_to_cor(block)
  corr = weight_Sigma(block, weights)
  max_ev = max(RSpectra::eigs_sym(corr, k = 1)$val)
  grid = data.frame(lambda = 10 ^ seq(log10(max_ev), log10(max_ev / 10000),
                                      length.out = n_lambda))
  rep(list(grid), length(outcomes))
}
cv_raw_regression = function(Y, family_ids, outcomes, covariates,
                             tuning_grid_fn, regression_fn, K = 5, cores = 1,
                             use_svd = FALSE) {
  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)
  n = nrow(Y)
  folds = make_family_folds(family_ids, K)
  if (use_svd) {
    s = svd(Y)
    Sigma_init = diag(s$d^2 / (n - 1))
    tuning_grids = tuning_grid_fn(Sigma_init, outcomes, covariates)
  } else {
    tuning_grids = tuning_grid_fn(cov(Y), outcomes, covariates)
  }
  n_tune = nrow(tuning_grids[[1]])

  predictions = array(NA_real_, dim = c(n, length(outcomes), n_tune))
  for (k in 1:K) {
    train = setdiff(1:n, folds[[k]])
    test = folds[[k]]
    Xtest = Y[test, covariates, drop = FALSE]
    if (use_svd) {
      s_train = svd(Y[train, ])
      Sigma_train = diag(s_train$d^2 / (length(train) - 1))
      V_train = s_train$v
      beta = regression_fn(Sigma_train, outcomes, covariates, tuning_grids, V = V_train, cores = cores)
    } else {
      Sigma_train = cov(Y[train, ])
      beta = regression_fn(Sigma_train, outcomes, covariates, tuning_grids, cores = cores)
    }
    for (l in 1:n_tune) predictions[test, , l] = Xtest %*% beta[, , l]
  }

  var_y = apply(Y[, outcomes, drop = FALSE], 2, var)
  cv_r2 = matrix(NA_real_, length(outcomes), n_tune)
  for (o in seq_along(outcomes)) {
    for (l in 1:n_tune) {
      cv_r2[o, l] = 1 - var(Y[, outcomes[o]] - predictions[, o, l]) / var_y[o]
    }
  }

  ## selected == 0L is the null model (beta = 0); see cv_latent_regression.
  selected = integer(length(outcomes))
  for (o in seq_along(outcomes)) {
    v = cv_r2[o, ]
    if (!any(is.finite(v))) {
      selected[o] = 0L
    } else {
      best = which.max(replace(v, !is.finite(v), -Inf))
      selected[o] = if (best == 1L && max(v, na.rm = TRUE) <= 0) 0L else best
    }
  }

  list(cv_r2 = cv_r2, selected = selected, tuning_grids = tuning_grids)
}
nested_cv_raw_regression = function(Y, family_ids, outcomes, covariates,
                                    tuning_grid_fn, regression_fn,
                                    K_outer = 10, K_inner = 5, cores = 1,
                                    use_svd = FALSE) {
  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)
  n = nrow(Y)
  folds = make_family_folds(family_ids, K_outer)

  predictions = matrix(NA_real_, n, length(outcomes))
  selected = matrix(NA_integer_, K_outer, length(outcomes))

  for (k in 1:K_outer) {
    print(k)
    train = setdiff(1:n, folds[[k]])

    inner = cv_raw_regression(Y[train, ], family_ids[train], outcomes, covariates,
                              tuning_grid_fn, regression_fn, K = K_inner, cores = cores,
                              use_svd = use_svd)
    selected[k, ] = inner$selected

    Xtest = Y[folds[[k]], covariates, drop = FALSE]
    if (use_svd) {
      s_train = svd(Y[train, ])
      Sigma_train = diag(s_train$d^2 / (length(train) - 1))
      V_train = s_train$v
      beta = regression_fn(Sigma_train, outcomes, covariates, inner$tuning_grids, V = V_train, cores = cores)
    } else {
      Sigma_train = cov(Y[train, ])
      beta = regression_fn(Sigma_train, outcomes, covariates, inner$tuning_grids, cores = cores)
    }
    for (o in seq_along(outcomes)) {
      sel = inner$selected[o]
      predictions[folds[[k]], o] = if (sel == 0L) 0 else Xtest %*% beta[, o, sel]
    }
  }

  var_y = apply(Y[, outcomes, drop = FALSE], 2, var)
  cv_r2 = numeric(length(outcomes))
  for (o in seq_along(outcomes)) {
    cv_r2[o] = 1 - var(Y[, outcomes[o]] - predictions[, o]) / var_y[o]
  }

  list(cv_r2 = cv_r2, selected = selected, predictions = predictions)
}
