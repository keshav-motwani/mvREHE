MixedEM <- function(data, K_G) {
  
  # Labels subjects
  K_G_labels = as.factor(rownames(K_G))
  
  weak_strong_cor = getWeakStrongCorDecomp(K_G);
  
  ## genetic identifiers (in theory same for genetically identical subjects)
  id_genetic = weak_strong_cor$id_genetic;  
  
  ## generate enviromental identifier (left over variability from the genetic one)                          
  id_envir =   gl(length(K_G_labels), 1)  
  #id_envir = sample(id_envir)
  
  #Prepare dataset for lme4
  Ydata = data.frame(data, id_genetic=id_genetic, id_envir=id_envir)  ##Generate multivariate dataset
  mYdata = melt(data.frame(Ydata, obs = seq(nrow(Ydata))), id.var = c("obs", "id_genetic", "id_envir"), variable.name = "variable")
  mYdata = cbind(id_unique = 1:nrow(mYdata),mYdata)
  
  # Structure that contains weak row correlations (\rho < 1) due to genetic and environmental factors
  corr = NULL;
  corr$id_genetic = weak_strong_cor$K_G_weak;
  
  lf = lFormulaMV(value ~ 0 + (0 + variable | id_genetic) + (0 + variable | id_envir), data = mYdata, cov_re = corr)
  
  #lf_check <- lmer(lf$formula, data = mYdata)
  #lf_check
  #lf_model_param = names(getME(lf_check, "theta"))   #Get number of parameters to optimize on
  # lf$reTrms$Zt = t(Z)       ## Modify accordingly lmer object
  control = lmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                  optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE));
  dfun <- do.call(mkLmerDevfun,c(lf, list(control = control)))          # Obtain deviance function
  
  
  # #### lmer optimization
  opt <- optimizeLmer(dfun, optimizer = control$optimizer, control = control$optCtrl, calc.derivs=control$calc.derivs)    # Optimize function
  # For faster options (see https://stats.stackexchange.com/questions/132841/default-lme4-optimizer-requires-lots-of-iterations-for-high-dimensional-data)
  # and https://cran.r-project.org/web/packages/lme4/vignettes/lmerperf.html
  # opt <- optimizeLmer(dfun, optimizer = "nloptwrap2")
  
  fit <- mkMerMod(environment(dfun), opt, lf$reTrms,fr = lf$fr) # Fit model
  vc_lme = VarCorr(fit)                                         # Get Variance-Covariance
  #print(vc_lme)
  
  # Construct covariance matrix - environment
  corr_envir_no_const_diag = as.matrix(attr(vc_lme$id_envir,"correlation"));
  std_dev_envir = as.vector(attr(vc_lme$id_envir,"stddev"));
  res_envir = as.numeric(attr(vc_lme,"sc"));
  cov_envir = std_dev_envir%*%t(std_dev_envir)*corr_envir_no_const_diag + diag(res_envir,nrow(corr_envir_no_const_diag),ncol(corr_envir_no_const_diag))^2;
  corr_envir = cov2cor(cov_envir);
  
  # Construct covariance matrix - genetics
  corr_genetic = as.matrix(attr(vc_lme$id_genetic,"correlation"));
  std_dev_genetic = as.vector(attr(vc_lme$id_genetic,"stddev"));
  cov_genetic = std_dev_genetic%*%t(std_dev_genetic)*corr_genetic;
  
  colnames(cov_genetic) = colnames(data);
  rownames(cov_genetic) = colnames(data);
  colnames(corr_genetic) = colnames(data);
  rownames(corr_genetic) = colnames(data);
  colnames(cov_envir) = colnames(data);
  rownames(cov_envir) = colnames(data);
  colnames(corr_envir) = colnames(data);
  rownames(corr_envir) = colnames(data);
  
  out_list = list(corr_genetic=corr_genetic, cov_genetic=cov_genetic, corr_envir=corr_envir, cov_envir=cov_envir,VarCorr=vc_lme, Z = t(lf$reTrms$Zt))
  return(out_list)
}


MixedEM_approx <- function(data, K_G, epsilon = 10e-4) {
  # Here we approximate strong correlations (rho=1) with
  # (rho=1-epsilon) to simplify implementation
  
  # Labels subjects
  K_G_labels = as.factor(rownames(K_G))
  
  # Kinship correlations approximating matrix
  K_G[K_G==1] = 1-epsilon;
  K_G = K_G + diag(epsilon,dim(K_G));
  
  ## generate genetic identifiers (in theory same for genetically identical subjects), in this implementation
  # we we will instead model genetically identical subjects with correlations (rho=1-epsilon)
  id_genetic = K_G_labels;  
  
  ## generate enviromental identifier (left over variability from the genetic one)                          
  id_envir =   gl(length(K_G_labels), 1)  
  #id_envir = sample(id_envir)
  
  #Prepare dataset for lme4
  Ydata = data.frame(data, id_genetic=id_genetic, id_envir=id_envir)  ##Generate multivariate dataset
  mYdata = melt(data.frame(Ydata, obs = seq(nrow(Ydata))), id.var = c("obs", "id_genetic", "id_envir"), variable.name = "variable")
  mYdata = cbind(id_unique = 1:nrow(mYdata),mYdata)
  
  # Structure that contains weak row correlations (\rho < 1) due to genetic and environmental factors
  corr = NULL;
  corr$id_genetic = K_G;
  
  lf = lFormulaMV(value ~ 0 + (0 + variable | id_genetic) + (0 + variable | id_envir), data = mYdata, cov_re = corr)
  
  #lf_check <- lmer(lf$formula, data = mYdata)
  #lf_check
  #lf_model_param = names(getME(lf_check, "theta"))   #Get number of parameters to optimize on
  # lf$reTrms$Zt = t(Z)       ## Modify accordingly lmer object
  dfun <- do.call(mkLmerDevfun,lf)                   # Obtain deviance function
  # 
  # #### lmer optimization
  opt <- optimizeLmer(dfun, optimizer = "bobyqa")    # Optimize function
  # Faster option (see https://stats.stackexchange.com/questions/132841/default-lme4-optimizer-requires-lots-of-iterations-for-high-dimensional-data)
  #opt <- optimizeLmer(dfun, optimizer = "nloptwrap2")
  
  fit <- mkMerMod(environment(dfun), opt, lf$reTrms,fr = lf$fr) # Fit model
  vc_lme = VarCorr(fit)                                         # Get Variance-Covariance
  #print(vc_lme)
  
  # Construct covariance matrix - environment
  corr_envir_no_const_diag = as.matrix(attr(vc_lme$id_envir,"correlation"));
  std_dev_envir = as.vector(attr(vc_lme$id_envir,"stddev"));
  res_envir = as.numeric(attr(vc_lme,"sc"));
  cov_envir = std_dev_envir%*%t(std_dev_envir)*corr_envir_no_const_diag + diag(res_envir,nrow(corr_envir_no_const_diag),ncol(corr_envir_no_const_diag))^2;
  corr_envir = cov2cor(cov_envir);
  
  # Construct covariance matrix - genetics
  corr_genetic = as.matrix(attr(vc_lme$id_genetic,"correlation"));
  std_dev_genetic = as.vector(attr(vc_lme$id_genetic,"stddev"));
  cov_genetic = std_dev_genetic%*%t(std_dev_genetic)*corr_genetic;
  
  # Naming
  colnames(cov_genetic) = colnames(data);
  rownames(cov_genetic) = colnames(data);
  colnames(corr_genetic) = colnames(data);
  rownames(corr_genetic) = colnames(data);
  colnames(cov_envir) = colnames(data);
  rownames(cov_envir) = colnames(data);
  colnames(corr_envir) = colnames(data);
  rownames(corr_envir) = colnames(data);
  
  out_list = list(corr_genetic=corr_genetic, cov_genetic=cov_genetic, corr_envir=corr_envir, cov_envir=cov_envir,VarCorr=vc_lme)
  return(out_list)
}

MixedEM_sbat <- function(data, K_G) {
  # This uses Franchini & Co implementation

  write.table(as.matrix(K_G),file="sbat/kinship_sim.txt",row.names=FALSE, col.names = FALSE, sep = " ")
  write.table(data,file="sbat/phenotypes_sim.txt",row.names=FALSE, col.names = TRUE, sep = " ")

  system("sbat/executables/sbat_centos6.8_static --covariances-only -p sbat/phenotypes_sim.txt -k sbat/kinship_sim.txt -o sbat/output_sim.txt", wait=TRUE)
  
  prec_genetic = read.table("sbat/output_sim.txt.C.mpmm.txt", header = FALSE)
  prec_envir = read.table("sbat/output_sim.txt.D.mpmm.txt", header = FALSE)
  
  cov_genetic = solve(prec_genetic);
  corr_genetic = cov2cor(cov_genetic);
  
  cov_envir = solve(prec_envir);
  corr_envir = cov2cor(cov_envir);
  
  # Naming
  colnames(cov_genetic) = colnames(data);
  rownames(cov_genetic) = colnames(data);
  colnames(corr_genetic) = colnames(data);
  rownames(corr_genetic) = colnames(data);
  colnames(cov_envir) = colnames(data);
  rownames(cov_envir) = colnames(data);
  colnames(corr_envir) = colnames(data);
  rownames(corr_envir) = colnames(data);
  
  out_list = list(corr_genetic=corr_genetic, cov_genetic=cov_genetic, corr_envir=corr_envir, cov_envir=cov_envir)
  return(out_list)
}

defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-16,maxeval=1e6)
nloptwrap2 <- function(fn,par,lower,upper,control=list(),...) {
  print("nloptwrap2")
  for (n in names(defaultControl)) 
    if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
    res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
    with(res,list(par=solution,
                  fval=objective,
                  feval=iterations,
                  conv=if (status>0) 0 else status,
                  message=message))
}


plotCorr <- function(corr) {
  ## ggplot2 plot of correlation
  melted_corr <- melt(corr)
  ggplot(data = melted_corr, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()
}

getWeakStrongCorDecomp <- function(K_G){
  # Kinship correlations approximating matrix
  K_G_strong = K_G;
  K_G_strong[K_G_strong<(0.98)] = 0;
  # K_G_strong
  # unique(c(as.matrix(K_G_strong)))
  
  g = graph_from_adjacency_matrix(K_G_strong, mode = "undirected", diag = FALSE)
  #plot(g)
  cl <- components(g)
  cl
  
  # loop through to extract common vertices
  #lapply(seq_along(cl$csize)[cl$csize > 1], function(x) V(g)$name[cl$membership %in% x])
  # loop through to one vertex per group
  list_conn_comp = lapply(seq_along(cl$csize), function(x) cl$membership[cl$membership %in% x])
  list_rep_comp = unlist(lapply(list_conn_comp, function(x) names(x[1])),use.names=FALSE)
  K_G_weak = K_G[list_rep_comp,list_rep_comp];
  colnames(K_G_weak) = as.character(seq_along(cl$csize));
  rownames(K_G_weak) = as.character(seq_along(cl$csize));
  # unique(c(as.matrix(triu(K_G_weak, k = 1))))
  
  ## id_genetic:  genetic identifiers (in theory same for genetically identical subjects)
  return(list(K_G_weak = K_G_weak, id_genetic = as.factor(cl$membership)))
}


getReplicatedKinship<-function(K_G,n_repl){
  
  mlist <- replicate(n_repl, K_G, simplify = FALSE)
  K_G_large = bdiag(mlist);
  replicate(n_repl, K_G, simplify = FALSE)
  K_G_large_labels = make.unique(rep(rownames(K_G), n_repl), sep = ".");
  rownames(K_G_large) = as.character(K_G_large_labels);
  colnames(K_G_large) = as.character(K_G_large_labels);
  
  return(K_G_large)
}

# #### Alternative optimization
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