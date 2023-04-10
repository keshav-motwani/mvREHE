// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double loss3(const arma::mat Y, const arma::mat & X_tilde, const List & Sigma_list, arma::vec lambda) {

  double value = 0;

  R_xlen_t K = X_tilde.n_cols;
  R_xlen_t q = Y.n_cols;

  for (R_xlen_t j = 0; j < q; j++) {

    for (R_xlen_t m = 0; m <= j; m++) {

      arma::mat YtY = Y.col(j) * Y.col(m).t();
      arma::vec YtY_vec = YtY.as_col();

      arma::vec sigma(K);
      sigma.zeros();
      for (R_xlen_t k = 0; k < K; k++) {
        NumericMatrix Sigma_ = Sigma_list[k];
        arma::mat Sigma(Sigma_.begin(), Sigma_.nrow(), Sigma_.ncol(), false) ;
        sigma[k] = Sigma(j, m);
      }

      if (j == m) {
        value += arma::accu(arma::pow(YtY_vec - (X_tilde * sigma), 2));
      } else {
        value += 2 * arma::accu(arma::pow(YtY_vec - (X_tilde * sigma), 2));
      }

    }

  }

  value = value / (X_tilde.n_rows * q * q);

  for (R_xlen_t k = 0; k < K; k++) {
    NumericMatrix Sigma_ = Sigma_list[k];
    arma::mat Sigma(Sigma_.begin(), Sigma_.nrow(), Sigma_.ncol(), false) ;
    value += lambda(k) * arma::accu(arma::square(Sigma));
  }

  return value;

}

// [[Rcpp::export]]
void compute_W_list(const arma::mat & Y, const List & D_list, List & W_list) {

  R_xlen_t n = Y.n_rows;
  R_xlen_t K = D_list.size();

  for (R_xlen_t k = 0; k < K; k++) {
    NumericMatrix D_ = D_list[k];
    arma::mat D(D_.begin(), D_.nrow(), D_.ncol(), false);
    NumericMatrix W_ = W_list[k];
    arma::mat W(W_.begin(), W_.nrow(), W_.ncol(), false);
    for (R_xlen_t i = 0; i < n; i++) {
      for (R_xlen_t l = 0; l < n; l++) {
        if (D(i, l) > 1e-10) {
          W += D(i, l) * Y.row(i).t() * Y.row(l);
        }
      }
    }
  }

}
