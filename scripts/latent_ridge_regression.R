cov2cor_NA0 = function(cov) {
  cov_hat = cov2cor(cov)
  cov_hat[diag(cov) < .Machine$double.eps, ] = 0
  cov_hat[, diag(cov) < .Machine$double.eps] = 0
  cov_hat
}

latent_ridge_regression = function(fit,
                                   component,
                                   outcome,
                                   covariates,
                                   lambda) {
  if (is.null(fit$V)) {
    Sigma_hat = cov2cor_NA0(fit$Sigma_hat[[component]])
  } else {
    Sigma_hat = cov2cor_NA0(fit$V %*% fit$Sigma_hat[[component]] %*% t(fit$V))
  }
  Sigma_hat = Sigma_hat[c(outcome, covariates), c(outcome, covariates)]
  
  return(solve(Sigma_hat[-1, -1] + diag(lambda, length(covariates)), Sigma_hat[-1, 1]))
  
}

raw_ridge_regression = function(Y, outcome, covariates, lambda) {
  Sigma_hat = cor(Y)
  Sigma_hat[is.na(Sigma_hat)] = 0
  
  Sigma_hat = Sigma_hat[c(outcome, covariates), c(outcome, covariates)]
  
  return(solve(Sigma_hat[-1, -1] + diag(lambda, length(covariates)), Sigma_hat[-1, 1]))
  
}

cv_latent_ridge_regression = function(Y,
                                      D_list,
                                      fit,
                                      outcomes,
                                      covariates,
                                      n_lambda,
                                      estimator,
                                      K = 2,
                                      folds = NULL,
                                      cores = 8) {
  require(parallel)
  
  if (!is.matrix(Y))
    Y = matrix(Y, ncol = 1)
  
  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y) / K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }
  
  components = 1:length(D_list)
  
  lambda_grid = array(NA, dim = c(length(components), n_lambda))
  cv_r2 = array(0, dim = c(length(components), length(outcomes), n_lambda))
  
  for (c in 1:length(components)) {
    component = components[c]
    
    cov_hat = cov2cor_NA0(fit$Sigma_hat[[component]])
    
    max_eigenvalue = max(RSpectra::eigs_sym(cov_hat[covariates, covariates], k = 1)$val)
    lambda_grid[component, ] = 10 ^ seq(log10(max_eigenvalue / 10000),
                                        log10(max_eigenvalue),
                                        length.out = n_lambda)
  }
  
  for (k in 1:K) {
    fit_train = estimator(Y[-folds[[k]], ], D_list = lapply(D_list, function(D)
      D[-folds[[k]], -folds[[k]]]))
    fit_test = estimator(Y[folds[[k]], ], D_list = lapply(D_list, function(D)
      D[folds[[k]], folds[[k]]]))
    
    for (c in 1:length(components)) {
      component = components[c]
      
      cov_hat_train = cov2cor_NA0(fit_train$Sigma_hat[[component]])
      cov_hat_test = cov2cor_NA0(fit_test$Sigma_hat[[component]])
      
      grid = expand.grid(o_idx = seq_along(outcomes),
                         l_idx = seq_len(n_lambda))

      results_grid = mclapply(seq_len(nrow(grid)), function(idx) {
        row = grid[idx,]
        o_idx = row$o_idx
        l_idx = row$l_idx
        
        outcome = outcomes[o_idx]
        lambda_val = lambda_grid[c, l_idx]
        
        beta_hat = solve(cov_hat_train[covariates, covariates] +
                           diag(lambda_val, length(covariates), length(covariates)),
                         cov_hat_train[covariates, outcome])
        
        r2 = 1 - (
          cov_hat_test[outcome, outcome] -
            2 * cov_hat_test[outcome, covariates] %*% beta_hat +
            t(beta_hat) %*% cov_hat_test[covariates, covariates] %*% beta_hat
        ) / cov_hat_test[outcome, outcome]
        
        list(o_idx = o_idx,
             l_idx = l_idx,
             r2 = as.numeric(r2))
        
      }, mc.cores = cores)
      
      for (res in results_grid) {
        cv_r2[c, res$o_idx, res$l_idx] = cv_r2[c, res$o_idx, res$l_idx] + res$r2 / K
      }
      
    }
    
  }
  
  lambda = matrix(NA, length(components), length(outcomes))
  for (c in 1:length(components)) {
    for (o in 1:length(outcomes)) {
      idx = which.max(cv_r2[c, o, ])
      if (length(idx) > 0) {
        lambda[c, o] = lambda_grid[c, idx]
      } else {
        lambda[c, o] = NA
      }
    }
  }
  
  return(list(
    lambda = lambda,
    cv_r2 = pmax(apply(cv_r2, c(1, 2), max), 0),
    cv_r2_full = cv_r2
  ))
  
}

cv_raw_ridge_regression = function(Y,
                                   outcomes,
                                   covariates,
                                   n_lambda,
                                   K = 2,
                                   folds = NULL,
                                   cores = 8) {
  require(parallel)
  
  if (!is.matrix(Y))
    Y = matrix(Y, ncol = 1)
  
  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y) / K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }
  
  
  lambda_grid = numeric(n_lambda)
  cv_r2 = array(0, dim = c(length(outcomes), n_lambda))
  
  
  cov_hat = cor(Y)
  cov_hat[is.na(cov_hat)] = 0
  
  max_eigenvalue = max(eigen(cov_hat[covariates, covariates])$val)
  lambda_grid = 10 ^ seq(log10(max_eigenvalue / 10000),
                         log10(max_eigenvalue),
                         length.out = n_lambda)
  
  
  
  for (k in 1:K) {
    cov_hat_train = cor(Y[-folds[[k]], ])
    cov_hat_test = cor(Y[folds[[k]], ])
    
    cov_hat_train[is.na(cov_hat_train)] = 0
    cov_hat_test[is.na(cov_hat_test)] = 0
    
    
    grid = expand.grid(o_idx = seq_along(outcomes),
                       l_idx = seq_len(n_lambda))
    
    results_grid = mclapply(seq_len(nrow(grid)), function(idx) {
      row = grid[idx,]
      o_idx = row$o_idx
      l_idx = row$l_idx
      
      outcome = outcomes[o_idx]
      lambda_val = lambda_grid[l_idx]
      
      beta_hat = solve(cov_hat_train[covariates, covariates] +
                         diag(lambda_val, length(covariates), length(covariates)),
                       cov_hat_train[covariates, outcome])
      
      r2 = 1 - (
        cov_hat_test[outcome, outcome] -
          2 * cov_hat_test[outcome, covariates] %*% beta_hat +
          t(beta_hat) %*% cov_hat_test[covariates, covariates] %*% beta_hat
      ) / cov_hat_test[outcome, outcome]
      
      list(o_idx = o_idx,
           l_idx = l_idx,
           r2 = as.numeric(r2))
      
    }, mc.cores = cores)
    
    for (res in results_grid) {
      cv_r2[res$o_idx, res$l_idx] = cv_r2[res$o_idx, res$l_idx] + res$r2 / K
    }
    
  }
  
  lambda = matrix(NA, length(outcomes))
  for (o in 1:length(outcomes)) {
    idx = which.max(cv_r2[o, ])
    if (length(idx) > 0) {
      lambda[o] = lambda_grid[idx]
    } else {
      lambda[o] = NA
    }
  }
  
  return(list(
    lambda = lambda,
    cv_r2 = pmax(apply(cv_r2, 1, max), 0),
    cv_r2_full = cv_r2
  ))
  
}
