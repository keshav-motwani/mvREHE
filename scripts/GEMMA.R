read_lower_tri = function(file, q) {

  mat = matrix(0, q, q)
  mat[upper.tri(mat, diag = TRUE)] = scan(file)
  mat = mat + t(mat)
  diag(mat) = diag(mat) / 2

  mat

}

GEMMA = function(Y, D_list) {

  K = D_list[[2]]
  K_file = tempfile(fileext = ".txt")
  write.table(K, K_file, col.names = FALSE, row.names = FALSE)

  Y_file = tempfile(fileext = ".txt")
  write.table(Y, Y_file, col.names = FALSE, row.names = FALSE)

  geno = cbind(data.frame(A = "SNP", B = "A", C = "T"), as.list(rep_len(c(0, 1), nrow(Y))))
  geno_file = tempfile(fileext = ".txt")
  write.table(geno, geno_file, row.names = FALSE, col.names = FALSE, sep = ", ")

  out_name = basename(tempfile())

  gemma = paste0("/apps/gemma/0.98.5/bin/gemma -g ", geno_file, " -p ", Y_file, " -n ", paste(1:ncol(Y), collapse = " "), " -k ", K_file, " -lmm -o ", out_name)
  system(gemma)
  extract_G = paste0("sed -n '", 25, ",", 25 + ncol(Y) - 1, " p' ", "output/", out_name, ".log.txt > output/Sigma_G_", out_name, ".txt")
  system(extract_G)
  extract_E = paste0("sed -n '", 25 + 2 * ncol(Y) + 2, ",", 25 + 3 * ncol(Y) + 1, " p' ", "output/", out_name, ".log.txt > output/Sigma_E_", out_name, ".txt")
  system(extract_E)

  Sigma_G = read_lower_tri(paste0("output/Sigma_G_", out_name, ".txt"), ncol(Y))
  Sigma_E = read_lower_tri(paste0("output/Sigma_E_", out_name, ".txt"), ncol(Y))

  return(list(Sigma_hat = list(Sigma_E, Sigma_G)))

}

