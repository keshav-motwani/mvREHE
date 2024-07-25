#' Perform massively univariate analysis
#'
#' @param Y
#' @param D_list
#' @param estimator
#'
#' @return
#' @export
#'
#' @examples
univariate = function(Y, D_list, estimator) {

  q = ncol(Y)

  estimates = apply(Y, 2, FUN = estimator, D_list = D_list, simplify = FALSE)

  Sigma_hat = lapply(1:length(D_list), function(i) matrix(0, q, q))

  for (j in 1:q) {
    for (k in 1:length(D_list)) {
      Sigma_hat[[k]][j, j] = estimates[[j]]$Sigma_hat[[k]][1, 1]
    }
  }

  return(list(Sigma_hat = Sigma_hat))

}
