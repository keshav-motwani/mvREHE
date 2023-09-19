// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double loss(const arma::mat Y, const List & D_list, const List & Sigma_list, const arma::vec row_indices, const arma::vec col_indices, arma::vec lambda) {

  double value = 0;

  R_xlen_t K = D_list.size();
  R_xlen_t q = Y.n_cols;
  R_xlen_t s = row_indices.size();

  arma::mat YtY(q, q, arma::fill::zeros);

  int j = 0;
  int m = 0;

  for (R_xlen_t i = 0; i < s; i++) {

    j = row_indices(i);
    m = col_indices(i);

    YtY = Y.row(j).t() * Y.row(m);

    for (R_xlen_t k = 0; k < K; k++) {
      NumericMatrix Sigma_ = Sigma_list[k];
      arma::mat Sigma(Sigma_.begin(), Sigma_.nrow(), Sigma_.ncol(), false);
      NumericMatrix D_ = D_list[k];
      arma::mat D(D_.begin(), D_.nrow(), D_.ncol(), false);
      YtY = YtY - D(j, m) * Sigma;
    }

    if (j == m) {
      value += arma::accu(arma::square(YtY));
    } else {
      value += 2 * arma::accu(arma::square(YtY));
    }

  }

  value = value / (Y.n_rows * Y.n_rows);

  for (R_xlen_t k = 0; k < K; k++) {
    NumericMatrix Sigma_ = Sigma_list[k];
    arma::mat Sigma(Sigma_.begin(), Sigma_.nrow(), Sigma_.ncol(), false) ;
    value += lambda(k) * arma::accu(arma::square(Sigma));
  }

  return value;

}

// [[Rcpp::export]]
void compute_W_list(const arma::mat & Y, const List & D_list, List & W_list, const arma::vec row_indices, const arma::vec col_indices) {

  R_xlen_t s = row_indices.size();
  R_xlen_t q = Y.n_cols;
  R_xlen_t K = D_list.size();

  arma::mat YtY(q, q, arma::fill::zeros);

  int i = 0;
  int l = 0;

  for (R_xlen_t k = 0; k < K; k++) {

    NumericMatrix D_ = D_list[k];
    arma::mat D(D_.begin(), D_.nrow(), D_.ncol(), false);
    NumericMatrix W_ = W_list[k];
    arma::mat W(W_.begin(), W_.nrow(), W_.ncol(), false);

    for (R_xlen_t j = 0; j < s; j++) {

      i = row_indices(j);
      l = col_indices(j);

      YtY = Y.row(i).t() * Y.row(l);

      if (i == l) {
        W += D(i, l) * YtY;
      } else {
        W += D(i, l) * (YtY + YtY.t());
      }
    }

  }

}

// [[Rcpp::export]]
arma::vec compute_Y_tilde(const arma::mat & Y, const arma::vec row_indices, const arma::vec col_indices, int j, int m) {

  R_xlen_t s = row_indices.size();

  arma::vec Y_tilde(s, arma::fill::zeros);

  for (R_xlen_t i = 0; i < s; i++) {
    Y_tilde(i) = Y(row_indices(i), j) * Y(col_indices(i), m);
  }

  return Y_tilde;

}
