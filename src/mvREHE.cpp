// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double frobenius_inner_product(S4 A, S4 B) {

  IntegerVector A_p = A.slot("p");
  IntegerVector A_i = A.slot("i");
  NumericVector A_x = A.slot("x");

  IntegerVector B_p = B.slot("p");
  IntegerVector B_i = B.slot("i");
  NumericVector B_x = B.slot("x");

  int ncol = A_p.size() - 1;
  double result = 0.0;

  for (int j = 0; j < ncol; j++) {
    int A_start = A_p[j];
    int A_end = A_p[j + 1];

    int B_start = B_p[j];
    int B_end = B_p[j + 1];

    int i_A = A_start;
    int i_B = B_start;

    while (i_A < A_end && i_B < B_end) {
      if (A_i[i_A] == B_i[i_B]) {
        double prod = A_x[i_A] * B_x[i_B];
        if (A_i[i_A] == j) {
          result += prod;
        } else {
          result += 2.0 * prod;
        }
        i_A++;
        i_B++;
      } else if (A_i[i_A] < B_i[i_B]) {
        i_A++;
      } else {
        i_B++;
      }
    }
  }

  return result;
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
