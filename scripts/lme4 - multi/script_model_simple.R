###########################################################################################################
# Estimation of phenotypic and additive genetic covariance functions for aligned function-valued traits   #
###########################################################################################################
rm(list=ls(all=TRUE))

library(devtools)
library(optimx)
library(reshape2)
library(Matrix)

load_all("lme4")

source("functions/mkBlistMV.R")
source("functions/mkReTrmsMV.R")
source("functions/lformulaMV.R")


n_genes = 500       # Number of individuals
n_per_genes = 1   # Repeated obervations (sharing same individual's random effects across repetition) (think of quantity measured in time with mean as random effect)
p = 3   # "Repeated obervations" *not* sharing same individual's random effects across repetition, i.e. multivariate response
                  # (think of different quantities measured in time)

n_families = n_genes*n_per_genes;
id_genetic = gl(n_genes, n_per_genes)  ## generate genetic identifiers (same for genetically identical subjects)
id_envir =   gl(n_families,n_genes*n_per_genes/n_families) ## generate enviromental identifier
#id_envir = sample(id_envir)

#K = read.csv("kinship.csv", header = F) # This is 2*K with Solar-Eclipse notation
#K = Matrix(K, sparse = TRUE)
#colnames(K) = NULL
#rownames(K) = NULL

# Genetic Covariance structure across individuals (between data matrix ROWS)
#K_G = as.matrix(K[c(1,2,535,3,824,4,5,24,6,7),c(1,2,535,3,824,4,5,24,6,7)])*0.8 + 0.2*diag(n_genes)   # Covariance between **individuals** (n_genes*n_genes)
#K_G = K

#set.seed(0)
n_DZ = 500
K_G = Diagonal(n_genes)
rownames(K_G) = as.character(levels(id_genetic))
colnames(K_G) = as.character(levels(id_genetic))
ind = sample(1:n_genes, n_DZ, replace=F)
ind_pairs_DZ = cbind(ind[1:(n_DZ/2)],ind[(n_DZ/2+1):n_DZ]);
K_G[ind_pairs_DZ] = 1-10e-10;
K_G[cbind(ind_pairs_DZ[,2],ind_pairs_DZ[,1])] = 1-10e-10;

Sq_K_G = t(chol(K_G))
# Sq_K_G%*%t(Sq_K_G)

# Genetic Covariance structure between p multivariate random effects (between data matrix COLUMNS)
#Sigma_G = diag(p)   #Sigma_G = matrix(c() , p, p)
Sigma_G = diag(0.7, p) + 0.3
Sq_Sigma_G = t(chol(Sigma_G))
# Sq_Sigma_G%*%t(Sq_Sigma_G)

# Environmental Covariance structure across individuals (between data matrix ROWS)
K_E = Diagonal(n_families)
rownames(K_E) = as.character(levels(id_envir))
colnames(K_E) = as.character(levels(id_envir))
# n_E = 1000;
# ind_E = sample(1:n_families, n_E, replace=F)
# ind_pairs_E = cbind(ind_E[1:(n_E/2)],ind_E[(n_E/2+1):n_E]);
# K_E[ind_pairs_E] = 1-10e-10;
# K_E[cbind(ind_pairs_E[,2],ind_pairs_E[,1])] =1-10e-10;
# # Sq_K_E%*%t(Sq_K_E)
Sq_K_E = t(chol(K_E))

# Environmental Covariance structure between p multivariate random effects (between data matrix COLUMNS)
#Sigma_E = diag(p)   #matrix( c(1,1,1,1), p, p)
Sigma_E = diag(0.3, p) + 0.7
Sq_Sigma_E = t(chol(Sigma_E))  #C=R*R^T
# Sq_Sigma_E%*%t(Sq_Sigma_E)

# ?Covariance model: kronecker(Sigma_G,K_G) + kronecker(Sigma_E,K_E)

#Generate data
U_G_st = matrix(rnorm(n_genes*p,mean=0,sd=1),n_genes)
U_E_st = matrix(rnorm(n_families*p,mean=0,sd=1),n_families)
E_st = matrix(rnorm(n_genes*n_per_genes*p,mean=0,sd=2),n_genes*n_per_genes)

U_G = as.matrix(Sq_K_G) %*% U_G_st %*% t(Sq_Sigma_G)        # Genetic component      (with covariance strtucture)
U_E = as.matrix(Sq_K_E) %*% U_E_st %*% t(Sq_Sigma_E)                   # Enviromental component (no covariance strtucture)


Y_no_noise = 1*U_G[id_genetic,] + 1*U_E[id_envir,]
Y          = Y_no_noise                           #+ 1*E_st
colnames(Y) = paste("y.", 1:p, sep = "")

Ydata = data.frame(Y, id_genetic=id_genetic, id_envir=id_envir)  ##Generate multivariate dataset
mYdata = melt(data.frame(Ydata, obs = seq(nrow(Ydata))), id.var = c("obs", "id_genetic", "id_envir"), variable.name = "variable")
mYdata = cbind(id_unique = 1:nrow(mYdata),mYdata)

####### Define the mixed effects model #########
corr = NULL;
corr$id_genetic = K_G;

#corr$id_envir = K_E;

#lf = lFormulaMV(value ~ 0 + (0 + variable | id_genetic), data = mYdata)
#lf = lFormulaMV(value ~ 0 + (0 + variable | id_genetic), data = mYdata, cov_re = corr)
lf = lFormulaMV(value ~ 0 + (0 + variable | id_genetic) + (0 + variable | id_envir), data = mYdata, cov_re = corr)
lf_simple = lFormulaMV(value ~ 0 + (0 + variable | id_genetic) + (0 + variable | id_envir), data = mYdata)

#lf_check <- lmer(lf$formula, data = mYdata)
#lf_check
#lf_model_param = names(getME(lf_check, "theta"))   #Get number of parameters to optimize on
# lf$reTrms$Zt = t(Z)       ## Modify accordingly lmer object
dfun_simple <- do.call(mkLmerDevfun,lf_simple)     # Obtain deviance function
dfun <- do.call(mkLmerDevfun,lf)                   # Obtain deviance function
#
# #### lmer optimazation
#opt <- optimizeLmer(dfun, optimizer = "bobyqa")    # Optimize function
opt_simple <- optimizeLmer(dfun_simple)
opt <- optimizeLmer(dfun)

fit <- mkMerMod(environment(dfun), opt, lf$reTrms,fr = lf$fr) # Fit model
vc_lme = VarCorr(fit)                                         # Get Variance-Covariance
#print(vc_lme)

# Construct covariance matrix - environment
corr_envir = as.matrix(attr(vc_lme$id_envir,"correlation"));
std_dev_envir = as.vector(attr(vc_lme$id_envir,"stddev"));
res_envir = as.numeric(attr(vc_lme,"sc"));
cov_envir = std_dev_envir%*%t(std_dev_envir)*corr_envir + diag(res_envir,nrow(corr_envir),ncol(corr_envir))^2;

# Construct covariance matrix - genetics
corr_genetic = as.matrix(attr(vc_lme$id_genetic,"correlation"));
std_dev_genetic = as.vector(attr(vc_lme$id_genetic,"stddev"));
cov_genetic = std_dev_genetic%*%t(std_dev_genetic)*corr_genetic;

print(cov_envir)
print(cov_genetic)

fit_simple <- mkMerMod(environment(dfun_simple), opt, lf_simple$reTrms,fr = lf_simple$fr) # Fit model
vc_lme_simple= VarCorr(fit_simple)                                         # Get Variance-Covariance
#print(vc_lme_simple)

# #### Alternative optimazation
# # opt1 <- optimx(par = rep(1,length(lf_model_param)), fn = dfun, method = "bobyqa")
# #
# # theta <- as.numeric(opt1[seq(length(lf_model_param))])     ## estimated Cholesky factors
# # pwrss <- opt1$value  ## penalized weighted residual sum of sq
# # n_mData <- nrow(mYdata)
# # cnms <- lf_check@cnms ## named list of names of variables per RE
# # vc <- lme4:::mkVarCorr(sc=pwrss/n_mData,
# #                        cnms=cnms,
# #                        nc=sapply(cnms,length), ## vars per random effect term
# #                        theta=theta,
# #                        nms=names(cnms))
# # attr(vc,"useSc") <- TRUE
# # class(vc) <- "VarCorr.merMod"
# # Sigma_g_est<-as.matrix(vc$id)[1:p,1:p]
#
# #vc_lme
#
# # Compare "true" (column) correlation with estimated
# #attr(vc$id,"correlation")                   # Estimated with manual minimazation
# cov_g_rmle = attr(vc_lme$id_genetic,"correlation")               # Estimated with lme minimazation
# cov_g_data = cov2cor(t(U_G)%*%solve(K_G)%*%U_G) #MLE estimation of covariance
#
# cov_e_rmle = attr(vc_lme$id_envir,"correlation")               # Estimated with lme minimazation
# cov_e_rmle_mod = cov2cor(cov_e_rmle + diag(attr(vc_lme, "sc")/attr(vc_lme$id_envir, "stddev")))
# cov_e_data = cor(U_E)
#
# print("Genetic Component - RMLE estimates\n")
# print(cov_g_rmle)
# print(cov_g_data)
# print(cor(U_G))
# cat("\n||cor(hat{U}_G)-cor(U_G)|| = ", norm(cov_g_rmle-cov_g_data, "F"),"\n")
#
# print("Enviroment Component - RMLE estimates\n")
# print(cov_e_rmle_mod)
# print(cov_e_data)
# cat("\n||cor(hat{U}_E)-cor(U_E)|| = ", norm(cov_e_rmle_mod-cov_e_data, "F"),"\n")
#
# print("Random Effects\n")
# cat("||hat{Y}-Y|| =               ", norm((Y_no_noise-fitted(fit))/dim(Y_no_noise)[1]))
#
#
# # Compare "true" (column) covariance with estimated
# # Attention difference of a constant factor (in covariance) between the two
# # vc$id[1:p, 1:p]       # Estimated with manual minimazation
# # vc_lme$id[1:p, 1:p]   # Estimated with lme minimazation
# # cov(U_G)
