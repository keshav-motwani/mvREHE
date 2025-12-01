Rcpp::sourceCpp("scripts/lasso_from_cov.cpp")

latent_lasso_regression = function(fit, component, outcome, covariates, lambda_seq, ...) {
  
  full_beta = matrix(0, length(covariates), length(lambda_seq))
  
  indices = diag(fit$Sigma_hat[[component]])[covariates] > 1e-12
  
  covariates = covariates[indices]
  
  if (is.null(fit$V)) {
    Sigma_hat = cov2cor_NA0(fit$Sigma_hat[[component]])
  } else {
    Sigma_hat = cov2cor_NA0(fit$V %*% fit$Sigma_hat[[component]] %*% t(fit$V))
  }
  Sigma_hat = Sigma_hat[c(outcome, covariates), c(outcome, covariates)]
  
  beta = lasso_from_cov(Sigma_hat, lambda_seq, ...)
  
  full_beta[indices, ] = beta
  
  return(full_beta)
  
}

raw_lasso_regression = function(Y, outcome, covariates, lambda_seq, ...) {
  
  covariates = covariates[matrixStats::colVars(Y[, covariates]) > 1e-12]
  
  Sigma_hat = cov2cor_NA0(cov(Y))

  Sigma_hat = Sigma_hat[c(outcome, covariates), c(outcome, covariates)]
  
  return(lasso_from_cov(Sigma_hat, lambda_seq, ...))
  
}

cv_latent_lasso_regression = function(Y, D_list, fit, outcomes, covariates, n_lambda, estimator, K = 2, folds = NULL, cores = 8, ...) {
  
  require(parallel)
  
  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)
  
  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y)/K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }
  
  covariates = covariates[matrixStats::colVars(Y[, covariates]) > 1e-12]
  
  components = 1:length(D_list)
  
  lambda_grid = array(NA, dim = c(length(components), length(outcomes), n_lambda))
  cv_r2 = array(0, dim = c(length(components), length(outcomes), n_lambda))
  
  for (c in 1:length(components)) {
    
    component = components[c]
    
    cov_hat = cov2cor_NA0(fit$Sigma_hat[[component]])

    for (o in 1:length(outcomes)) {
      
      outcome = outcomes[o]
      
      if (cov_hat[outcome, outcome] > 0) {
        lambda_max = max(abs(cov_hat[covariates, outcome]))
        lambda_grid[c, o, ] = exp(seq(log(lambda_max), log(lambda_max * 1e-4), length.out = n_lambda))
      } else {
        lambda_grid[c, o, ] = NA
      }
      
    }
    
  }
  
  for (k in 1:K) {
    
    fit_train = estimator(Y[-folds[[k]], ], D_list = lapply(D_list, function(D) D[-folds[[k]], -folds[[k]]]))
    fit_test = estimator(Y[folds[[k]], ], D_list = lapply(D_list, function(D) D[folds[[k]], folds[[k]]]))
    
    for (c in 1:length(components)) {
      
      component = components[c]
      
      cov_hat_train = cov2cor_NA0(fit_train$Sigma_hat[[component]])
      cov_hat_test = cov2cor_NA0(fit_test$Sigma_hat[[component]])

      results_outcomes = mclapply(1:length(outcomes), function(o_idx) {
        
        outcome = outcomes[o_idx]
        
        if (cov_hat_train[outcome, outcome] > 0) { 
        
          beta_path = lasso_from_cov(
            cov_hat_train[c(outcome, covariates), c(outcome, covariates)],
            lambda_grid[c, o_idx, ],
            ...
          )
          
          r2_values = numeric(n_lambda)
          for (l in 1:n_lambda) {
            beta_hat = beta_path[, l]
            r2 = 1 - (cov_hat_test[outcome, outcome] -
                        2 * cov_hat_test[outcome, covariates] %*% beta_hat +
                        t(beta_hat) %*% cov_hat_test[covariates, covariates] %*% beta_hat) /
              cov_hat_test[outcome, outcome]
            r2_values[l] = r2
          }
          
          list(o_idx = o_idx, r2_values = r2_values)
          
        } else {
          
          list(o_idx = o_idx, r2_values = NA)
          
        }
        
      }, mc.cores = cores)
      
      for (res in results_outcomes) {
        cv_r2[c, res$o_idx, ] = cv_r2[c, res$o_idx, ] + res$r2_values / K
      }
      
    }
    
  }
  
  lambda = matrix(NA, length(components), length(outcomes))
  for (c in 1:length(components)) {
    for (o in 1:length(outcomes)) {
      idx = which.max(cv_r2[c, o, ])
      if (length(idx) > 0) {
        lambda[c, o] = lambda_grid[c, o, ][idx]
      } else {
        lambda[c, o] = NA
      }
    }
  }
  
  return(list(lambda = lambda, cv_r2 = pmax(apply(cv_r2, c(1, 2), max), 0), cv_r2_full = cv_r2))
  
}


cv_raw_lasso_regression = function(Y, outcomes, covariates, n_lambda, K = 2, folds = NULL, cores = 8, ...) {
  
  require(parallel)
  
  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)
  
  if (is.null(folds)) {
    folds = split(1:nrow(Y), rep(1:K, each = ceiling(nrow(Y)/K)))
  } else {
    stopifnot(length(setdiff(1:nrow(Y), unlist(folds))) == 0)
    stopifnot(length(setdiff(unlist(folds), 1:nrow(Y))) == 0)
  }
  
  covariates = covariates[matrixStats::colVars(Y[, covariates]) > 1e-12]
  
  lambda_grid = array(NA, dim = c(length(outcomes), n_lambda))
  cv_r2 = array(0, dim = c(length(outcomes), n_lambda))
  
  cov_hat = cov2cor_NA0(cov(Y))

  for (o in 1:length(outcomes)) {
    
    outcome = outcomes[o]
    
    if (cov_hat[outcome, outcome] > 0) {
      lambda_max = max(abs(cov_hat[covariates, outcome]))
      lambda_grid[o, ] = exp(seq(log(lambda_max), log(lambda_max * 1e-4), length.out = n_lambda))
    } else {
      lambda_grid[o, ] = NA
    }
    
  }
  
  
  for (k in 1:K) {
    
    
    cov_hat_train = cov2cor_NA0(cov(Y[-folds[[k]], ]))
    cov_hat_test = cov2cor_NA0(cov(Y[folds[[k]], ]))

    
    results_outcomes = mclapply(seq_along(outcomes), function(o_idx) {
      
      outcome = outcomes[o_idx]
      
      if (cov_hat_train[outcome, outcome] > 0) {
      
        beta_path = lasso_from_cov(
          cov_hat_train[c(outcome, covariates), c(outcome, covariates)],
          lambda_grid[o_idx, ],
          ...
        )
        
        r2_values = numeric(n_lambda)
        for (l in 1:n_lambda) {
          beta_hat = beta_path[, l]
          r2 = 1 - (cov_hat_test[outcome, outcome] -
                      2 * cov_hat_test[outcome, covariates] %*% beta_hat +
                      t(beta_hat) %*% cov_hat_test[covariates, covariates] %*% beta_hat) /
            cov_hat_test[outcome, outcome]
          r2_values[l] = r2
        }
        
        list(o_idx = o_idx, r2_values = r2_values)
      
      } else {
        
        list(o_idx = o_idx, r2_values = NA)
        
      }
      
    }, mc.cores = cores)
    
    for (res in results_outcomes) {
      cv_r2[res$o_idx, ] = cv_r2[res$o_idx, ] + res$r2_values / K
    }
    
  }
  
  
  lambda = numeric(length(outcomes))
  
  for (o in 1:length(outcomes)) {
    idx = which.max(cv_r2[o, ])
    if (length(idx) > 0) {
      lambda[o] = lambda_grid[o, ][idx]
    } else {
      lambda[o] = NA
    }
  }
  
  
  return(list(lambda = lambda, cv_r2 = pmax(apply(cv_r2, 1, max), 0), cv_r2_full = cv_r2))
  
}
