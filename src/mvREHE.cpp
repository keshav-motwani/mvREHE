// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double loss(const List & Y_tilde_list, const arma::mat & X_tilde, const List & L_list, arma::vec lambda) {

  double value = 0;

  R_xlen_t K = X_tilde.n_cols;
  R_xlen_t q = Y_tilde_list.size();

  for (R_xlen_t j = 0; j < q; j++) {

    List Y_tilde_list_j = Y_tilde_list[j];

    for (R_xlen_t m = 0; m <= j; m++) {

      NumericVector Y_ = Y_tilde_list_j[m];
      arma::vec Y(Y_.begin(), Y_.size(), false, true);

      arma::vec sigma(K);
      sigma.zeros();
      for (R_xlen_t k = 0; k < K; k++) {
        NumericMatrix L_ = L_list[k];
        arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false) ;
        sigma[k] = arma::dot(L.row(j).subvec(0, std::min(j, m)), L.row(m).subvec(0, std::min(j, m)));
      }

      if (j == m) {
        value += arma::accu(arma::pow(Y - (X_tilde * sigma), 2));
      } else {
        value += 2 * arma::accu(arma::pow(Y - (X_tilde * sigma), 2));
      }

    }

  }

  value = value / (X_tilde.n_rows * q * q);

  for (R_xlen_t k = 0; k < K; k++) {
    NumericMatrix L_ = L_list[k];
    arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false) ;
    value += lambda(k) / 2 * arma::accu(arma::square(L));
  }

  return value;

}

// [[Rcpp::export]]
double loss2(const List & Y_tilde_list, const arma::mat & X_tilde, const List & Sigma_list) {

  double value = 0;

  R_xlen_t K = X_tilde.n_cols;
  R_xlen_t q = Y_tilde_list.size();

  for (R_xlen_t j = 0; j < q; j++) {

    List Y_tilde_list_j = Y_tilde_list[j];

    for (R_xlen_t m = 0; m <= j; m++) {

      NumericVector Y_ = Y_tilde_list_j[m];
      arma::vec Y(Y_.begin(), Y_.size(), false, true);

      arma::vec sigma(K);
      sigma.zeros();
      for (R_xlen_t k = 0; k < K; k++) {
        NumericMatrix Sigma_ = Sigma_list[k];
        arma::mat Sigma(Sigma_.begin(), Sigma_.nrow(), Sigma_.ncol(), false) ;
        sigma[k] = Sigma(j, m);
      }

      if (j == m) {
        value += arma::accu(arma::pow(Y - (X_tilde * sigma), 2));
      } else {
        value += 2 * arma::accu(arma::pow(Y - (X_tilde * sigma), 2));
      }

    }

  }

  return value / (X_tilde.n_rows * q * q);

}

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
void gradient_full(const List & Y_tilde_list, const arma::mat & X_tilde, const List & L_list, List & gradient_list, arma::vec lambda) {

  R_xlen_t K = X_tilde.n_cols;
  R_xlen_t q = Y_tilde_list.size();

  arma::vec sigma(K);
  sigma.zeros();
  arma::vec R(X_tilde.n_rows);
  R.zeros();
  double ZtR = 0;

  for (R_xlen_t a = 0; a < q; a++) {
    for (R_xlen_t k = 0; k < K; k++) {
      NumericMatrix L_ = L_list[k];
      arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false) ;
      sigma[k] = arma::dot(L.row(a).subvec(0, a), L.row(a).subvec(0, a));
    }
    List Y_tilde_list_a = Y_tilde_list[a];
    NumericVector Y_ = Y_tilde_list_a[a];
    arma::vec Y(Y_.begin(), Y_.size(), false, true);
    R = Y - X_tilde * sigma;
    for (R_xlen_t z = 0; z < K; z++) {
      ZtR = arma::dot(X_tilde.col(z), R);
      for (R_xlen_t b = 0; b <= a; b++) {
        NumericMatrix L_ = L_list[z];
        arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false);
        NumericMatrix grad_ = gradient_list[z];
        arma::mat grad(grad_.begin(), grad_.nrow(), grad_.ncol(), false);
        grad(a, b) = -4 * L(a, b) * ZtR / (X_tilde.n_rows * q * q);
      }
    }
  }

  for (R_xlen_t a = 0; a < q; a++) {
    List Y_tilde_list_a = Y_tilde_list[a];
    for (R_xlen_t m = 0; m <= a - 1; m++) {
      for (R_xlen_t k = 0; k < K; k++) {
        NumericMatrix L_ = L_list[k];
        arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false) ;
        sigma[k] = arma::dot(L.row(a).subvec(0, m), L.row(m).subvec(0, m));
      }
      NumericVector Y_ = Y_tilde_list_a[m];
      arma::vec Y(Y_.begin(), Y_.size(), false, true);
      R = Y - X_tilde * sigma;
      for (R_xlen_t z = 0; z < K; z++) {
        ZtR = arma::dot(X_tilde.col(z), R);
        for (R_xlen_t b = 0; b <= a; b++) {
          if (m >= b) {
            NumericMatrix L_ = L_list[z];
            arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false);
            NumericMatrix grad_ = gradient_list[z];
            arma::mat grad(grad_.begin(), grad_.nrow(), grad_.ncol(), false);
            grad(a, b) += -4 * L(m, b) * ZtR / (X_tilde.n_rows * q * q);
          }
        }
      }
    }
  }

  for (R_xlen_t a = 0; a < q; a++) {
    for (R_xlen_t j = a + 1; j < q; j++) {
      for (R_xlen_t k = 0; k < K; k++) {
        NumericMatrix L_ = L_list[k];
        arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false) ;
        sigma[k] = arma::dot(L.row(j).subvec(0, a), L.row(a).subvec(0, a));
      }
      List Y_tilde_list_j = Y_tilde_list[j];
      NumericVector Y_ = Y_tilde_list_j[a];
      arma::vec Y(Y_.begin(), Y_.size(), false, true);
      R = Y - X_tilde * sigma;
      for (R_xlen_t z = 0; z < K; z++) {
        ZtR = arma::dot(X_tilde.col(z), R);
        for (R_xlen_t b = 0; b <= a; b++) {
          NumericMatrix L_ = L_list[z];
          arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false);
          NumericMatrix grad_ = gradient_list[z];
          arma::mat grad(grad_.begin(), grad_.nrow(), grad_.ncol(), false);
          grad(a, b) += -4 * L(j, b) * ZtR / (X_tilde.n_rows * q * q);
        }
      }
    }
  }

  for (R_xlen_t k = 0; k < K; k++) {
    NumericMatrix grad_ = gradient_list[k];
    arma::mat grad(grad_.begin(), grad_.nrow(), grad_.ncol(), false);
    NumericMatrix L_ = L_list[k];
    arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false);
    grad += lambda(k) * L;
  }

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
