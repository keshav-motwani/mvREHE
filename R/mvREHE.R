#' mvREHE
#'
#' @param Y
#' @param D_list
#' @param lambda
#' @param tolerance
#' @param max_iter
#' @param Sigma_init_list
#' @param W_list
#' @param Q
#' @param row_indices
#' @param col_indices
#'
#' @return
#' @export
#'
#' @examples
mvREHE = function(Y, D_list, lambda = NULL, tolerance = 1e-3, max_iter = 100, Sigma_init_list = NULL, W_list = NULL, Q = NULL, row_indices = NULL, col_indices = NULL) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  n = nrow(Y)
  q = ncol(Y)
  K = length(D_list)
  objective = numeric(max_iter)

  if (is.null(lambda)) lambda = rep(0, K)
  if (is.null(Sigma_init_list)) {
    Sigma_list = lapply(1:length(D_list), function(i) diag(1, q, q))
  } else if (is.character(Sigma_init_list) && Sigma_init_list == "mvHE") {
    Sigma_list = mvHE(Y, D_list)$Sigma_hat
  } else {
    Sigma_list = Sigma_init_list
  }
  if (is.null(row_indices) && (!is.null(tolerance) | is.null(W_list))) {
    indices = do.call(`+`, D_list) > 0
    indices = indices & lower.tri(indices, diag = TRUE)
    row_indices = which(indices, arr.ind = TRUE)[, 1] - 1
    col_indices = which(indices, arr.ind = TRUE)[, 2] - 1
  }
  if (is.null(W_list)) {
    W_list = lapply(1:K, function(x) matrix(0, q, q))
    compute_W_list(Y, D_list, W_list, row_indices, col_indices)
  }
  if (is.null(Q)) {
    Q = compute_Q(D_list)
  }

  for (iter in 1:max_iter) {

    Sigma_list_old = Sigma_list

    for (z in 1:K) {
      mat = W_list[[z]]
      for (k in setdiff(1:K, z)) {
        mat = mat - Sigma_list[[k]] * Q[k, z]
      }
      eig = eigen(mat, symmetric = TRUE)
      Sigma_list[[z]] = eig$vectors %*% (t(eig$vectors) * pmax(eig$values / (Q[z, z] + (lambda[z] * n^2)), 0))
    }

    if (!is.null(tolerance)) {
      # objective[iter] = loss(Y, D_list, Sigma_list, lambda, row_indices, col_indices)
      # if (iter > 1 && abs(objective[iter - 1] - objective[iter]) / objective[iter - 1] < tolerance) {
      #   break
      # }
      objective[iter] = mean(mapply(Sigma_list_old, Sigma_list, FUN = function(x, y) norm(x - y, "F") / norm(x, "F")), na.rm = TRUE)
      if (iter > 1 && objective[iter] < tolerance) {
        break
      }
    }

  }

  return(list(Sigma_hat = Sigma_list,
              objective = objective[objective != 0]))

}

#' mvREHE_DR
#'
#' @param Y
#' @param D_list
#' @param r
#' @param lambda
#' @param tolerance
#' @param max_iter
#' @param Sigma_init_list
#' @param Q
#' @param row_indices
#' @param col_indices
#'
#' @return
#' @export
#'
#' @examples
mvREHE_DR = function(Y, D_list, V = NULL, tolerance = NULL, max_iter = 100, Sigma_init_list = NULL, Q = NULL, row_indices = NULL, col_indices = NULL) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  if (is.null(V)) {

    V = floor(ncol(Y) / 2)

  }

  if (!is.matrix(V)) {

    V = min(V, min(nrow(Y), ncol(Y)))

    s = svd(Y)
    V = s$v[, 1:V]

  }

  Yr = Y %*% V

  estimate = mvREHE(Yr, D_list, lambda = NULL, tolerance, max_iter, Sigma_init_list, W_list = NULL, Q, row_indices, col_indices)

  estimate$Sigma_hat = lapply(estimate$Sigma_hat, function(x) V %*% x %*% t(V))

  estimate

}

#' mvREHE_cvDR
#'
#' @param Y
#' @param D_list
#' @param V_seq
#' @param K
#' @param n_rank
#' @param rank_max
#' @param rank_min
#' @param tolerance
#' @param max_iter
#' @param Sigma_init_list
#'
#' @return
#' @export
#'
#' @examples
mvREHE_cvDR = function(Y, D_list, V_seq, K = 5, folds = NULL, tolerance = 1e-3, max_iter = 100) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y)/K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }

  q = ncol(Y)

  precomputed_values = precompute_cv(Y, D_list, folds, compute_W = FALSE)

  cv = with(precomputed_values, {

    loss = numeric(length(V_seq))

    cv_loss = function(V) {
      loss_l = 0
      for (k in 1:K) {
        Sigma_hat = mvREHE_DR(Y[-folds[[k]], , drop = FALSE], D_list_mk_list[[k]], V = V, tolerance, max_iter, Sigma_init_list = NULL, Q = Q_mk_list[[k]], row_indices = row_indices_mk_list[[k]], col_indices = col_indices_mk_list[[k]])$Sigma_hat
        loss_l = loss_l + length(folds[[k]])^2 * loss(Y[folds[[k]], , drop = FALSE], D_list_k_list[[k]], Sigma_hat, row_indices_k_list[[k]], col_indices_k_list[[k]], rep(0, length(D_list)))
      }
      return(loss_l)
    }

    for (l in 1:length(V_seq)) {
      cat(l)
      loss[l] = cv_loss(V_seq[[l]])
    }

    l = which.min(loss)
    V = V_seq[[l]]

    cat("\n")
    print(V)

    list(V = V, loss = loss)

  })

  V = cv$V
  loss = cv$loss

  result = mvREHE_DR(Y, D_list, V = V, tolerance, max_iter, Sigma_init_list = NULL)
  result$cv_loss = loss
  result$V = V

  result

}

#' mvREHE_cvL2
#'
#' @param Y
#' @param D_list
#' @param K
#' @param n_rank
#' @param rank_max
#' @param rank_min
#' @param tolerance
#' @param max_iter
#' @param Sigma_init_list
#'
#' @return
#' @export
#'
#' @examples
mvREHE_cvL2 = function(Y, D_list, K = 5, folds = NULL, grid = TRUE, n_lambda = 5, lambda_max = 1e-2, lambda_min = 1e-6, tolerance = 1e-3, max_iter = 100, Sigma_init_list = NULL) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y)/K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }

  q = ncol(Y)

  Sigma_init_list = mvREHE(Y, D_list, lambda = NULL, tolerance, max_iter, Sigma_init_list)$Sigma_hat

  precomputed_values = precompute_cv(Y, D_list, folds)

  cv = with(precomputed_values, {

    lambda_seq = c(10^seq(log10(lambda_max), log10(lambda_min), length.out = n_lambda), 0)
    loss = numeric(length(lambda_seq))

    cv_loss = function(lambda) {
      loss_l = 0
      for (k in 1:K) {
        cat(k)
        Sigma_hat = mvREHE(Y[-folds[[k]], , drop = FALSE], D_list_mk_list[[k]], lambda = lambda, tolerance, max_iter, Sigma_init_list, W_list = W_list_mk_list[[k]], Q = Q_mk_list[[k]], row_indices = row_indices_mk_list[[k]], col_indices = col_indices_mk_list[[k]])$Sigma_hat
        loss_l = loss_l + length(folds[[k]])^2 * loss(Y[folds[[k]], , drop = FALSE], D_list_k_list[[k]], Sigma_hat, row_indices_k_list[[k]], col_indices_k_list[[k]], rep(0, length(D_list)))
      }
      return(loss_l)
    }

    if (grid) {

      for (l in 1:length(lambda_seq)) {
        cat(l)
        loss[l] = cv_loss(rep(lambda_seq[l], length(D_list)))
      }

      l = which.min(loss)
      lambda = rep(lambda_seq[l], length(D_list))

    } else {

      lambda = 10^coef(optimx::optimx(rep(-5, length(D_list)), function(lambda) cv_loss(10^lambda), method = "L-BFGS-B", control = list(trace = 0, kkt =FALSE, starttests = FALSE)))

    }

    cat("\n")
    print(lambda)

    list(lambda_seq = lambda_seq, lambda = lambda, loss = loss)

  })

  lambda_seq = cv$lambda_seq
  lambda = cv$lambda
  loss = cv$loss

  result = mvREHE(Y, D_list, lambda = lambda, tolerance, max_iter, Sigma_init_list)
  result$cv_loss = loss
  result$lambda = lambda
  if (grid) result$lambda_seq = lambda_seq

  result

}

precompute_cv = function(Y, D_list, folds, compute_W = TRUE) {

  K = length(folds)
  q = ncol(Y)

  D_list_mk_list = vector(mode = "list", length = K)
  W_list_mk_list = vector(mode = "list", length = K)
  Q_mk_list = vector(mode = "list", length = K)
  D_list_k_list = vector(mode = "list", length = K)
  row_indices_mk_list = vector(mode = "list", length = K)
  col_indices_mk_list = vector(mode = "list", length = K)
  row_indices_k_list = vector(mode = "list", length = K)
  col_indices_k_list = vector(mode = "list", length = K)
  for (k in 1:K) {
    D_list_k_list[[k]] = lapply(D_list, function(D) D[folds[[k]], folds[[k]]])
    D_list_mk_list[[k]] = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]])
    Q_mk_list[[k]] = compute_Q(D_list_mk_list[[k]])
    indices = do.call(`+`, D_list_k_list[[k]]) > 0
    indices = indices & lower.tri(indices, diag = TRUE)
    row_indices_k_list[[k]] = which(indices, arr.ind = TRUE)[, 1] - 1
    col_indices_k_list[[k]] = which(indices, arr.ind = TRUE)[, 2] - 1
    indices = do.call(`+`, D_list_mk_list[[k]]) > 0
    indices = indices & lower.tri(indices, diag = TRUE)
    row_indices_mk_list[[k]] = which(indices, arr.ind = TRUE)[, 1] - 1
    col_indices_mk_list[[k]] = which(indices, arr.ind = TRUE)[, 2] - 1
    if (compute_W) {
      W_list_mk_list[[k]] = lapply(1:length(D_list), function(x) matrix(0, q, q))
      compute_W_list(Y[-folds[[k]], , drop = FALSE], D_list_mk_list[[k]], W_list_mk_list[[k]], row_indices_mk_list[[k]], col_indices_mk_list[[k]])
    }
  }

  return(list(D_list_mk_list = D_list_mk_list,
              W_list_mk_list = W_list_mk_list,
              Q_mk_list = Q_mk_list,
              D_list_k_list = D_list_k_list,
              row_indices_mk_list = row_indices_mk_list,
              col_indices_mk_list = col_indices_mk_list,
              row_indices_k_list = row_indices_k_list,
              col_indices_k_list = col_indices_k_list))

}

compute_Q = function(D_list) {

  K = length(D_list)

  Q = matrix(NA, K, K)

  for (i in 1:K) {
    for (j in 1:i) {
      Q[i, j] = Q[j, i] = sum(D_list[[i]] * D_list[[j]])
    }
  }

  return(Q)

}
