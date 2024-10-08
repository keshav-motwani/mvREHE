library(devtools)
library(optimx)
library(reshape2)
library(Matrix)

library(lme4)
attach(loadNamespace("lme4"), name = "lme4_all")

source("scripts/lme4 - multi/functions/mkBlistMV.R")
source("scripts/lme4 - multi/functions/mkReTrmsMV.R")
source("scripts/lme4 - multi/functions/lformulaMV.R")

mvREML = function(Y, D_list) {

  if (length(D_list) == 2) {
    estimate = mvREML_2(Y, D_list)
  } else if (length(D_list) == 3) {
    estimate = mvREML_3(Y, D_list)
  }

  if (any(sapply(estimate$Sigma_hat, function(Sigma) any(is.na(Sigma))))) {
    estimate$Sigma_hat = lapply(estimate$Sigma_hat, function(Sigma) {
      Sigma[] = NA
      Sigma
    })
  }

  estimate

}

mvREML_3 = function(Y, D_list) {

  D_0 = D_list[[1]]
  D_1 = D_list[[2]]
  D_2 = D_list[[3]]

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  colnames(Y) = paste("y.", 1:ncol(Y), sep = "")

  expand = diag(1, nrow(D_1), nrow(D_1))
  expand = expand[, !duplicated(D_1)]
  expand[cbind(which(duplicated(D_1)), which(duplicated(D_1)) - 1:sum(duplicated(D_1)))] = 1
  id_genetic = apply(expand, 1, function(x) which(x == 1))
  D_1 = D_1[!duplicated(D_1), !duplicated(D_1)]
  colnames(D_1) = rownames(D_1) = unique(id_genetic)

  expand = diag(1, nrow(D_2), nrow(D_2))
  expand = expand[, !duplicated(D_2)]
  expand[cbind(which(duplicated(D_2)), which(duplicated(D_2)) - 1:sum(duplicated(D_2)))] = 1
  id_common_env = apply(expand, 1, function(x) which(x == 1))
  D_2 = D_2[!duplicated(D_2), !duplicated(D_2)]
  colnames(D_2) = rownames(D_2) = unique(id_common_env)

  id_noise = 1:nrow(Y)

  Ydata = data.frame(Y, id_genetic = id_genetic, id_common_env = id_common_env, id_noise = id_noise)
  mYdata = melt(data.frame(Ydata, obs = seq(nrow(Ydata))), id.var = c("obs", "id_genetic", "id_common_env", "id_noise"), variable.name = "variable")
  mYdata = cbind(id_unique = 1:nrow(mYdata),mYdata)

  corr = NULL
  corr$id_genetic = D_1
  corr$id_common_env = D_2

  lf = lFormulaMV(value ~ 0 + (0 + variable | id_genetic) + (0 + variable | id_common_env) + (0 + variable | id_noise), data = mYdata, cov_re = corr)
  dfun <- do.call(mkLmerDevfun, lf)
  opt <- optimizeLmer(dfun)

  fit <- mkMerMod(environment(dfun), opt, lf$reTrms,fr = lf$fr) # Fit model
  vc_lme = VarCorr(fit)                                         # Get Variance-Covariance

  corr_noise = as.matrix(attr(vc_lme$id_noise,"correlation"));
  std_dev_noise = as.vector(attr(vc_lme$id_noise,"stddev"));
  res_noise = as.numeric(attr(vc_lme,"sc"));
  cov_noise = std_dev_noise%*%t(std_dev_noise)*corr_noise + diag(res_noise,nrow(corr_noise),ncol(corr_noise))^2;

  corr_genetic = as.matrix(attr(vc_lme$id_genetic,"correlation"));
  std_dev_genetic = as.vector(attr(vc_lme$id_genetic,"stddev"));
  cov_genetic = std_dev_genetic%*%t(std_dev_genetic)*corr_genetic;

  corr_common_env = as.matrix(attr(vc_lme$id_common_env,"correlation"));
  std_dev_common_env = as.vector(attr(vc_lme$id_common_env,"stddev"));
  cov_common_env = std_dev_common_env%*%t(std_dev_common_env)*corr_common_env;

  return(list(Sigma_hat = list(cov_noise, cov_genetic, cov_common_env)))

}

mvREML_2 = function(Y, D_list) {

  D_0 = D_list[[1]]
  D_1 = D_list[[2]]

  if (!is.matrix(Y)) Y = matrix(Y, ncol = 1)

  colnames(Y) = paste("y.", 1:ncol(Y), sep = "")

  expand = diag(1, nrow(D_1), nrow(D_1))
  expand = expand[, !duplicated(D_1)]
  expand[cbind(which(duplicated(D_1)), which(duplicated(D_1)) - 1:sum(duplicated(D_1)))] = 1
  id_genetic = apply(expand, 1, function(x) which(x == 1))
  D_1 = D_1[!duplicated(D_1), !duplicated(D_1)]
  colnames(D_1) = rownames(D_1) = unique(id_genetic)

  id_envir = 1:nrow(Y)

  Ydata = data.frame(Y, id_genetic = id_genetic, id_envir = id_envir)
  mYdata = melt(data.frame(Ydata, obs = seq(nrow(Ydata))), id.var = c("obs", "id_genetic", "id_envir"), variable.name = "variable")
  mYdata = cbind(id_unique = 1:nrow(mYdata),mYdata)

  corr = NULL;
  corr$id_genetic = D_1;

  lf = lFormulaMV(value ~ 0 + (0 + variable | id_genetic) + (0 + variable | id_envir), data = mYdata, cov_re = corr)
  dfun <- do.call(mkLmerDevfun, lf)
  opt <- optimizeLmer(dfun)

  fit <- mkMerMod(environment(dfun), opt, lf$reTrms,fr = lf$fr) # Fit model
  vc_lme = VarCorr(fit)                                         # Get Variance-Covariance

  corr_envir = as.matrix(attr(vc_lme$id_envir,"correlation"));
  std_dev_envir = as.vector(attr(vc_lme$id_envir,"stddev"));
  res_envir = as.numeric(attr(vc_lme,"sc"));
  cov_envir = std_dev_envir%*%t(std_dev_envir)*corr_envir + diag(res_envir,nrow(corr_envir),ncol(corr_envir))^2;

  corr_genetic = as.matrix(attr(vc_lme$id_genetic,"correlation"));
  std_dev_genetic = as.vector(attr(vc_lme$id_genetic,"stddev"));
  cov_genetic = std_dev_genetic%*%t(std_dev_genetic)*corr_genetic;

  return(list(Sigma_0_hat = cov_envir, Sigma_1_hat = cov_genetic, Sigma_hat = list(cov_envir, cov_genetic)))

}
