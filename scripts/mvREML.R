library(glmmTMB)
library(reshape2)
library(Matrix)

mvREML = function(Y, D_list) {

  if (length(D_list) == 2) {
    estimate = mvREML_inner(Y, D_list[1:2])
  } else if (length(D_list) == 3) {
    estimate = mvREML_inner(Y, D_list)
  }

  if (any(sapply(estimate$Sigma_hat, function(Sigma) any(is.na(Sigma))))) {
    estimate$Sigma_hat = lapply(estimate$Sigma_hat, function(Sigma) {
      Sigma[] = NA
      Sigma
    })
  }

  estimate

}

resolve_dups = function(D, digits = 10) {
  row_keys = apply(round(D, digits), 1, paste, collapse = "\t")
  unique_mask = !duplicated(row_keys)
  ids = match(row_keys, row_keys[unique_mask])
  list(ids = ids, D = D[unique_mask, unique_mask, drop = FALSE])
}

to_sparse_sym = function(M) {
  as(as(as(M, "dMatrix"), "symmetricMatrix"), "CsparseMatrix")
}

ensure_psd = function(M, eps = 1e-8) {
  e = eigen(M, symmetric = TRUE)
  e$vectors %*% diag(pmax(e$values, eps)) %*% t(e$vectors)
}

mvREML_inner = function(Y, D_list) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)
  q = ncol(Y)
  n = nrow(Y)
  colnames(Y) = paste0("y.", 1:q)

  # D_list[[1]] is the noise covariance (not used — handled by id_noise term)
  # D_list[[2]] is the genetic kinship
  # D_list[[3]] (if present) is the common env kinship

  D_1_raw = as.matrix(D_list[[2]])
  res_1 = resolve_dups(D_1_raw)
  D_1 = res_1$D
  rownames(D_1) = colnames(D_1) = seq_len(nrow(D_1))
  id_genetic = factor(res_1$ids, levels = rownames(D_1))

  has_common_env = length(D_list) >= 3
  if (has_common_env) {
    D_2_raw = as.matrix(D_list[[3]])
    res_2 = resolve_dups(D_2_raw)
    D_2 = res_2$D
    rownames(D_2) = colnames(D_2) = seq_len(nrow(D_2))
    id_common_env = factor(res_2$ids, levels = rownames(D_2))
  }

  id_noise = factor(1:n)

  if (has_common_env) {
    Ydata = data.frame(Y, id_genetic = id_genetic, id_common_env = id_common_env, id_noise = id_noise)
    mYdata = melt(Ydata, id.var = c("id_genetic", "id_common_env", "id_noise"), variable.name = "variable")
  } else {
    Ydata = data.frame(Y, id_genetic = id_genetic, id_noise = id_noise)
    mYdata = melt(Ydata, id.var = c("id_genetic", "id_noise"), variable.name = "variable")
  }

  n_gen = nlevels(id_genetic)
  L_1 = t(chol(to_sparse_sym(D_1)))
  L_kron_1 = kronecker(L_1, Diagonal(q))

  if (has_common_env) {
    n_com = nlevels(id_common_env)
    L_2 = t(chol(to_sparse_sym(D_2)))
    L_kron_2 = kronecker(L_2, Diagonal(q))
  }

  if (has_common_env) {
    form = value ~ 0 +
      (0 + variable | id_genetic) +
      (0 + variable | id_common_env) +
      (0 + variable | id_noise)
    if (ncol(Y) == 1) {
      form = value ~ 0 +
        (1 | id_genetic) +
        (1 | id_common_env) +
        (1 | id_noise)
    }
  } else {
    form = value ~ 0 +
      (0 + variable | id_genetic) +
      (0 + variable | id_noise)
    if (ncol(Y) == 1) {
      form = value ~ 0 +
        (1 | id_genetic) +
        (1 | id_noise)
    }
  }

  m0 = glmmTMB(form, data = mYdata, REML = TRUE, doFit = FALSE, dispformula = ~0)

  cols_gen = 1:(n_gen * q)
  m0$data.tmb$Z[, cols_gen] = m0$data.tmb$Z[, cols_gen] %*% L_kron_1

  if (has_common_env) {
    cols_com = (n_gen * q + 1):(n_gen * q + n_com * q)
    m0$data.tmb$Z[, cols_com] = m0$data.tmb$Z[, cols_com] %*% L_kron_2
  }

  best_fit = tryCatch(suppressWarnings(glmmTMB:::fitTMB(m0, doOptim = TRUE)), error = function(e) NULL)

  na_mat = matrix(NA_real_, q, q)
  if (is.null(best_fit) || is.null(best_fit$fit)) {
    if (has_common_env)
      return(list(Sigma_hat = list(na_mat, na_mat, na_mat)))
    else
      return(list(Sigma_hat = list(na_mat, na_mat)))
  }

  vc = VarCorr(best_fit)
  Sigma_noise   = as.matrix(vc$cond$id_noise)
  Sigma_genetic = as.matrix(vc$cond$id_genetic)

  if (has_common_env) {
    Sigma_common_env = as.matrix(vc$cond$id_common_env)
    return(list(Sigma_hat = list(Sigma_noise, Sigma_genetic, Sigma_common_env)))
  } else {
    return(list(Sigma_hat = list(Sigma_noise, Sigma_genetic)))
  }

}
