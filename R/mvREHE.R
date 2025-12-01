#' mvREHE
#'
#' @param Y
#' @param D_list
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
mvREHE = function(Y, D_list, tolerance = 1e-6, max_iter = 1000, return_full = TRUE, Sigma_init_list = NULL, W_list = NULL, Q = NULL, row_indices = NULL, col_indices = NULL) {

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  n = nrow(Y)
  q = ncol(Y)
  K = length(D_list)
  difference = numeric(max_iter)

  highdim = q > n

  if (highdim) {

    s = svd(Y)
    Y = Y %*% s$v
    q = ncol(Y)

  }

  if (is.null(Sigma_init_list)) {
    Sigma_list = lapply(1:length(D_list), function(i) matrix(0, q, q))
  } else if (is.character(Sigma_init_list) && Sigma_init_list == "mvHE") {
    Sigma_list = mvHE(Y, D_list)$Sigma_hat
  } else {
    Sigma_list = Sigma_init_list
  }
  if (is.null(row_indices) && (!is.null(tolerance) | is.null(W_list))) {
    indices = abs(Reduce(`+`, D_list)) > .Machine$double.eps
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

    print(iter)

    Sigma_list_old = Sigma_list

    for (z in 1:K) {
      mat = W_list[[z]]
      for (k in setdiff(1:K, z)) {
        mat = mat - Sigma_list[[k]] * Q[k, z]
      }
      eig = eigen(mat)
      Sigma_list[[z]] = eig$vectors %*% (t(eig$vectors) * pmax(c(eig$values) / Q[z, z], 0))
    }

    if (!is.null(tolerance)) {
      # objective[iter] = loss(Y, D_list, Sigma_list, row_indices, col_indices)
      # if (iter > 1 && abs(objective[iter - 1] - objective[iter]) / objective[iter - 1] < tolerance) {
      #   break
      # }
      difference[iter] = mean(mapply(Sigma_list_old, Sigma_list, FUN = function(x, y) norm(x - y, "F") / norm(x, "F")), na.rm = TRUE)
      if (iter > 1 && is.na(difference[iter])) {
        break
      }
      if (iter > 1 && difference[iter] < tolerance) {
        break
      }
    }

  }

  result = list(difference = difference[difference != 0])

  if (highdim) {

    result$Sigma_hat = Sigma_list
    result$V = s$v

    if (return_full) {
      result$Sigma_hat = lapply(Sigma_list, function(Sigma_r) s$v %*% Sigma_r %*% t(s$v))
      result$V = NULL
    }

  } else {

    result$Sigma_hat = Sigma_list

  }

  return(result)

}

precompute_cv = function(Y, D_list, folds, compute_W = TRUE, V_function = NULL, r = NULL) {

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
  V_mk_list = vector(mode = "list", length = K)
  for (k in 1:K) {
    D_list_k_list[[k]] = lapply(D_list, function(D) D[folds[[k]], folds[[k]]])
    D_list_mk_list[[k]] = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]])
    Q_mk_list[[k]] = compute_Q(D_list_mk_list[[k]])
    indices = Reduce(`+`, D_list_k_list[[k]]) > 0
    indices = indices & lower.tri(indices, diag = TRUE)
    row_indices_k_list[[k]] = which(indices, arr.ind = TRUE)[, 1] - 1
    col_indices_k_list[[k]] = which(indices, arr.ind = TRUE)[, 2] - 1
    indices = Reduce(`+`, D_list_mk_list[[k]]) > 0
    indices = indices & lower.tri(indices, diag = TRUE)
    row_indices_mk_list[[k]] = which(indices, arr.ind = TRUE)[, 1] - 1
    col_indices_mk_list[[k]] = which(indices, arr.ind = TRUE)[, 2] - 1
    if (compute_W) {
      W_list_mk_list[[k]] = lapply(1:length(D_list), function(x) matrix(0, q, q))
      compute_W_list(Y[-folds[[k]], , drop = FALSE], D_list_mk_list[[k]], W_list_mk_list[[k]], row_indices_mk_list[[k]], col_indices_mk_list[[k]])
    }
    if (!is.null(V_function) & !is.null(r)) {
      V_mk_list[[k]] = V_function(Y[-folds[[k]], , drop = FALSE], r)
    }
  }

  return(list(D_list_mk_list = D_list_mk_list,
              W_list_mk_list = W_list_mk_list,
              Q_mk_list = Q_mk_list,
              D_list_k_list = D_list_k_list,
              row_indices_mk_list = row_indices_mk_list,
              col_indices_mk_list = col_indices_mk_list,
              row_indices_k_list = row_indices_k_list,
              col_indices_k_list = col_indices_k_list,
              V_mk_list = V_mk_list))

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
