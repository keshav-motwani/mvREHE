// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double loss(const arma::mat Y, const List & D_list, const List & Sigma_list, arma::vec lambda) {

  double value = 0;

  R_xlen_t K = D_list.size();
  R_xlen_t q = Y.n_cols;
  R_xlen_t n = Y.n_rows;

  arma::mat YtY(q, q, arma::fill::zeros);

  for (R_xlen_t j = 0; j < n; j++) {

    for (R_xlen_t m = 0; m <= j; m++) {

      int all0 = 1;
      for (R_xlen_t k = 0; k < K; k++) {
        NumericMatrix D_ = D_list[k];
        arma::mat D(D_.begin(), D_.nrow(), D_.ncol(), false);
        if (D(j, m) > 1e-10) {
          all0 = 0;
          break;
        }
      }

      if (all0 == 0) {

        YtY = Y.row(j).t() * Y.row(m);

        for (R_xlen_t k = 0; k < K; k++) {
          NumericMatrix Sigma_ = Sigma_list[k];
          arma::mat Sigma(Sigma_.begin(), Sigma_.nrow(), Sigma_.ncol(), false);
          NumericMatrix D_ = D_list[k];
          arma::mat D(D_.begin(), D_.nrow(), D_.ncol(), false);
          if (D(j, m) > 1e-10) {
            YtY = YtY - D(j, m) * Sigma;
          }
        }

        if (j == m) {
          value += arma::accu(arma::square(YtY));
        } else {
          value += 2 * arma::accu(arma::square(YtY));
        }

      }

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
void compute_W_list(const arma::mat & Y, const List & D_list, List & W_list) {

  R_xlen_t n = Y.n_rows;
  R_xlen_t q = Y.n_cols;
  R_xlen_t K = D_list.size();

  arma::mat YtY(q, q, arma::fill::zeros);

  for (R_xlen_t i = 0; i < n; i++) {
    for (R_xlen_t l = 0; l <= i; l++) {

      int all0 = 1;
      for (R_xlen_t k = 0; k < K; k++) {
        NumericMatrix D_ = D_list[k];
        arma::mat D(D_.begin(), D_.nrow(), D_.ncol(), false);
        if (D(i, l) > 1e-10) {
          all0 = 0;
          break;
        }
      }

      if (all0 == 0) {
        YtY = Y.row(i).t() * Y.row(l);
        for (R_xlen_t k = 0; k < K; k++) {
          NumericMatrix D_ = D_list[k];
          arma::mat D(D_.begin(), D_.nrow(), D_.ncol(), false);
          NumericMatrix W_ = W_list[k];
          arma::mat W(W_.begin(), W_.nrow(), W_.ncol(), false);
          if (D(i, l) > 1e-10) {
            if (i == l) {
              W += D(i, l) * YtY;
            } else {
              W += D(i, l) * (YtY + YtY.t());
            }
          }
        }
      }

    }
  }

}

// [[Rcpp::export]]
arma::vec compute_Y_tilde(const arma::mat & Y, const::arma::vec indices, const arma::vec row_indices, const arma::vec col_indices, int j, int m) {

  R_xlen_t s = indices.size();

  arma::vec Y_tilde(s, arma::fill::zeros);

  for (R_xlen_t i = 0; i < s; i++) {
      Y_tilde(i) = Y(row_indices(i), j) * Y(col_indices(i), m);
  }

  return Y_tilde;

}
