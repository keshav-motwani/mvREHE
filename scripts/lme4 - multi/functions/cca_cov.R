cca_cov <- function(Sxy,Sx,Sy)
{
  p <- ncol(Sx) 
  q <- ncol(Sy) 
  
  
  Sxeig <- eigen(Sx, symmetric=TRUE)
  Sxisqrt <- Sxeig$vectors %*% diag(1/sqrt(Sxeig$values)) %*% t(Sxeig$vectors)
  Syeig <- eigen(Sy, symmetric=TRUE)
  Syisqrt <- Syeig$vectors %*% diag(1/sqrt(Syeig$values)) %*% t(Syeig$vectors)
  Xmat <- Sxisqrt %*% Sxy %*% solve(Sy) %*% t(Sxy) %*% Sxisqrt
  Ymat <- Syisqrt %*% t(Sxy) %*% solve(Sx) %*% Sxy %*% Syisqrt
  Xeig <- eigen(Xmat, symmetric=TRUE)
  Yeig <- eigen(Ymat, symmetric=TRUE)
  
  # correlations (alternatively sqrt(Yeig$values[1:min(p,q)]))
  rho <- sqrt(Xeig$values[1:min(p,q)])
  #rho
  
  # compare linear combinations (different!)
  Ahat <- Sxisqrt %*% Xeig$vectors
  Bhat <- Syisqrt %*% Yeig$vectors
  
  return(list(Ahat=Ahat,Bhat=Bhat,rho=rho))
}