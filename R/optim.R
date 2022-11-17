fit_optim = function(Y_tilde_list, X_tilde, lambda, max_iter = 10000, L_init_list = NULL, algorithm = "L-BFGS-B") {

  library(optimx)

  q = length(Y_tilde_list)
  K = ncol(X_tilde)

  loss_wrapper = function(par) {

    dim(par) = c(q*(q+1)/2,K)
    L_list = apply(par,2,function(cov_vec){
      cov_mat <- matrix(0,q,q)
      cov_mat[lower.tri(cov_mat, diag=TRUE)] <- cov_vec
      cov_mat
    }, simplify = F)

    value = loss(Y_tilde_list, X_tilde, L_list, lambda)

    return(value)
  }

  gradient_wrapper = function(par) {

    dim(par) = c(q*(q+1)/2,K)
    L_list = apply(par, 2, function(cov_vec){
      cov_mat <- matrix(0,q,q)
      cov_mat[lower.tri(cov_mat, diag=TRUE)] <- cov_vec
      cov_mat
    }, simplify = F)

    gradient_list = replicate(K, matrix(0, q, q), simplify = FALSE)
    gradient_full(Y_tilde_list, X_tilde, L_list, gradient_list, lambda)

    sapply(gradient_list, function(cov_mat){
      cov_mat[lower.tri(cov_mat, diag=TRUE)]
    })

  }

  if (is.null(L_init_list)) {
    L_init_list = lapply(1:K, function(i) t(chol(clusterGeneration::rcorrmatrix(q))))
  } else if (is.character(L_init_list) && L_init_list == "mvHE") {
    mvHE_estimate = mvHE(Y, D_list)$Sigma_hat
    if (all(sapply(mvHE_estimate, function(x) attr(x, "truncated")) == 0)) return(list(Sigma_hat = mvHE_estimate))
    L_init_list = lapply(1:length(mvHE_estimate), function(i) t(chol(mvHE_estimate[[i]], pivot = TRUE)))
  }

  par = c(sapply(L_init_list,function(cov_mat){
    cov_mat[lower.tri(cov_mat, diag=TRUE)]
  }))

  opt <- optimx(par, loss_wrapper,
                gr = gradient_wrapper,
                method = algorithm,
                control = list(trace = 0, maxit = max_iter, kkt =FALSE, starttests = FALSE) )
  par_opt <- coef(opt)

  dim(par_opt) = c(q*(q+1)/2,K)
  L_list = apply(par_opt,2,function(cov_vec){
    cov_mat <- matrix(0,q,q)
    cov_mat[lower.tri(cov_mat, diag=TRUE)] <- cov_vec
    cov_mat
  }, simplify = F)

  result = list(L_list = L_list, objective = loss_wrapper(par_opt))

  return(result)

}
