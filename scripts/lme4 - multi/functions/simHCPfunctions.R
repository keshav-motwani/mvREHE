oneSim <-function(p, mult_sample_size, K_G)
{
  # Generate large row genetic covariance
  K_G_large = getReplicatedKinship(K_G,mult_sample_size);
  
  # Decompose genetic covariance in weak and strong correlations
  weak_strong_cor = getWeakStrongCorDecomp(K_G_large);
  K_G_weak = weak_strong_cor$K_G_weak;
  
  id_genetic = weak_strong_cor$id_genetic;  ## generate genetic identifiers (same for genetically identical subjects)
  id_envir =   gl(nrow(K_G_large), 1)  ## generate enviromental identifier
  

  # Genetic Covariance structure between p multivariate random effects (between data matrix COLUMNS)
  Sigma_G = rcorrmatrix(p,alphad=1);
  # Environmental Covariance structure between p multivariate random effects (between data matrix COLUMNS)
  Sigma_E = rcorrmatrix(p,alphad=1);
  
  data = generateData(K_G_weak, id_genetic, id_envir, Sigma_G, Sigma_E)
  
  # Approx
  start.time <- Sys.time()
  mixed_out_approx = MixedEM_approx(data, K_G_large, epsilon = 10e-6);
  end.time <- Sys.time();
  time.taken <- end.time - start.time;
  # Exact
  start.time <- Sys.time()
  mixed_out = MixedEM(data, K_G_large);
  end.time <- Sys.time();
  time.taken <- end.time - start.time;
  print(time.taken)
  
  # SBAT
  #start.time <- Sys.time()
  #mixed_out_sbat = MixedEM_sbat(data, K_G_large);
  #end.time <- Sys.time();
  #time.taken <- end.time - start.time;
  #print(time.taken)
  
  
  err_gen = norm(mixed_out$corr_genetic-Sigma_G, "F");
  err_env = norm(mixed_out$corr_envir-Sigma_E, "F");
  err_gen_approx = norm(mixed_out_approx$corr_genetic-Sigma_G, "F");
  err_env_approx = norm(mixed_out_approx$corr_envir-Sigma_E, "F");
  #err_gen_sbat = norm(mixed_out_sbat$corr_genetic-Sigma_G, "F");
  #err_env_sbat = norm(mixed_out_sbat$corr_envir-Sigma_E, "F");
  
  #err_gen_sbat,err_env_sbat
  # N_rs = 10000
  # distr_random_gen = replicate(N_rs, randomEstimation(Sigma_G))
  # distr_random_env = replicate(N_rs, randomEstimation(Sigma_E))
  
  # qtl_gen = sum(distr_random_gen<=err_gen)/N_rs
  # qtl_env = sum(distr_random_env<=err_env)/N_rs
  # qtl_gen_approx = sum(distr_random_gen<=err_gen_approx)/N_rs
  # qtl_env_approx = sum(distr_random_env<=err_env_approx)/N_rs
  
  cat("Err genetic: ", err_gen,"\n")
  cat("Err envir",     err_env,"\n")
  cat("Err genetic approx: ", err_gen_approx,"\n")
  cat("Err envir approx",     err_env_approx,"\n")
  # cat("Qtl genetic: ", qtl_gen,"\n")
  # cat("Qtl envir",     qtl_env,"\n")
  # cat("Qtl genetic approx: ", qtl_gen_approx,"\n")
  # cat("Qtl envir approx",     qtl_env_approx,"\n")
  #cat("Err genetic sbat: ", err_gen_sbat,"\n")
  #cat("Err envir sbat",     err_env_sbat,"\n")
  
  return(c(err_gen,err_env,err_gen_approx,err_env_approx));
}

oneSim_exact_approx <-function(p, mult_sample_size, K_G)
{
  # Generate large row genetic covariance
  K_G_large = getReplicatedKinship(K_G,mult_sample_size);
  
  # Decompose genetic covariance in weak and strong correlations
  weak_strong_cor = getWeakStrongCorDecomp(K_G_large);
  K_G_weak = weak_strong_cor$K_G_weak;
  
  id_genetic = weak_strong_cor$id_genetic;  ## generate genetic identifiers (same for genetically identical subjects)
  id_envir =   gl(nrow(K_G_large), 1)  ## generate enviromental identifier
  
  
  # Genetic Covariance structure between p multivariate random effects (between data matrix COLUMNS)
  Sigma_G = rcorrmatrix(p,alphad=1);
  # Environmental Covariance structure between p multivariate random effects (between data matrix COLUMNS)
  Sigma_E = rcorrmatrix(p,alphad=1);
  
  data = generateData(K_G_weak, id_genetic, id_envir, Sigma_G, Sigma_E)
  
  # Approx
  start.time <- Sys.time()
  mixed_out_approx = MixedEM_approx(data, K_G_large, epsilon = 10e-6);
  end.time <- Sys.time();
  time.taken <- end.time - start.time;
  # Exact
  start.time <- Sys.time()
  mixed_out = MixedEM(data, K_G_large);
  end.time <- Sys.time();
  time.taken <- end.time - start.time;
  print(time.taken)
  
  err_gen = norm(mixed_out$corr_genetic-Sigma_G, "F");
  err_env = norm(mixed_out$corr_envir-Sigma_E, "F");
  err_gen_approx = norm(mixed_out_approx$corr_genetic-Sigma_G, "F");
  err_env_approx = norm(mixed_out_approx$corr_envir-Sigma_E, "F");
  
  cat("Err genetic: ", err_gen,"\n")
  cat("Err envir",     err_env,"\n")
  cat("Err genetic approx: ", err_gen_approx,"\n")
  cat("Err envir approx",     err_env_approx,"\n")
  
  return(c(err_gen,err_env, err_gen_approx, err_env_approx));
}

randomEstimation <-function(Sigma_G){
  return(norm(Sigma_G - rcorrmatrix(p,alphad=1),"F"))
}

generateData <-function(K_G_weak, id_genetic, id_envir, Sigma_G, Sigma_E){
  
  n_g = length(levels(id_genetic))
  n_e = length(levels(id_envir))
  
  if(any(dim(Sigma_G) != dim(Sigma_E))) stop('Sigma_G and Sigma_E have different dimensions')
  p = nrow(Sigma_G);
  
  Sq_Sigma_G = t(chol(Sigma_G))
  Sq_Sigma_E = t(chol(Sigma_E))  #C=R*R^T
  
  Sq_K_G = t(chol(K_G_weak)) #Try another one
  
  #Generate data
  U_G_st = matrix(rnorm(n_g*p,mean=0,sd=1),n_g)
  U_E_st = matrix(rnorm(n_e*p,mean=0,sd=1),n_e)
  U_G = as.matrix(Sq_K_G) %*% U_G_st %*% t(Sq_Sigma_G)        # Genetic component      (with covariance strtucture)
  U_E = U_E_st %*% t(Sq_Sigma_E)                   # Enviromental component (no covariance strtucture)
  
  Y = 1*U_G[id_genetic,] + U_E[id_envir,];
  
  colnames(Y) = paste("y.", 1:p, sep = "");
  return(Y)
}



