// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double loss(const List & Y_tilde_list, const arma::mat & X_tilde, const List & L_list, double lambda) {

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
    value += lambda / 2 * arma::accu(arma::square(L));
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
void gradient_full(const List & Y_tilde_list, const arma::mat & X_tilde, const List & L_list, List & gradient_list, double lambda) {

  double value = 0;

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
    grad += lambda * L;
  }

}

// [[Rcpp::export]]
arma::vec fit_GD(const List & Y_tilde_list, const arma::mat & X_tilde, double lambda, int max_iter, double tolerance, List & L_list, List & gradient_list) {

  arma::colvec objective(max_iter);
  objective.zeros();

  R_xlen_t q = Y_tilde_list.size();
  R_xlen_t K = X_tilde.n_cols;

  List L_list_old(Rcpp::clone(L_list));

  double loss_old = 0;
  double rhs = 0;
  double lhs = 0;

  double step_size = 1e-3;

  for (int i = 0; i < max_iter; i++) {

    for (R_xlen_t k = 0; k < K; k++) {
      NumericMatrix L_ = L_list[k];
      arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false);
      NumericMatrix L_old_ = L_list_old[k];
      arma::mat L_old(L_old_.begin(), L_old_.nrow(), L_old_.ncol(), false);
      L_old = L;
    }

    loss_old = loss(Y_tilde_list, X_tilde, L_list_old, lambda);
    gradient_full(Y_tilde_list, X_tilde, L_list_old, gradient_list, lambda);

    step_size = step_size * 2;
    bool line_search = true;

    while(line_search) {

      rhs = loss_old;
      for (R_xlen_t z = 0; z < K; z++) {
        NumericMatrix L_ = L_list[z];
        arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false);
        NumericMatrix L_old_ = L_list_old[z];
        arma::mat L_old(L_old_.begin(), L_old_.nrow(), L_old_.ncol(), false);
        NumericMatrix grad_ = gradient_list[z];
        arma::mat grad(grad_.begin(), grad_.nrow(), grad_.ncol(), false);
        L = L_old - step_size * grad;
        rhs = rhs - 0.5 * step_size * arma::accu(arma::pow(grad, 2));
      }
      lhs = loss(Y_tilde_list, X_tilde, L_list, lambda);
      if (std::isnan(lhs) | std::isnan(rhs) | (lhs > rhs)) {
        step_size = 0.5 * step_size;
      } else {
        line_search = false;
      }

    }

    objective(i) = loss(Y_tilde_list, X_tilde, L_list, lambda);

    if (i > 1 && ((objective(i - 1) - objective(i)) / objective(i - 1)) < tolerance) {
      break;
    }

  }

  return objective;

}

// // [[Rcpp::export]]
// double loss_update(const List & Y_tilde_list, const arma::mat & X_tilde, const List & L_list, int a, int b) {
//
//   double value = 0;
//
//   R_xlen_t K = X_tilde.n_cols;
//   R_xlen_t q = Y_tilde_list.size();
//
//   List Y_tilde_list_a = Y_tilde_list[a];
//
//   for (R_xlen_t m = b; m <= a; m++) {
//
//     NumericVector Y_ = Y_tilde_list_a[m];
//     arma::vec Y(Y_.begin(), Y_.size(), false, true);
//
//     arma::vec sigma(K);
//     sigma.zeros();
//     for (R_xlen_t k = 0; k < K; k++) {
//       NumericMatrix L_ = L_list[k];
//       arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false) ;
//       sigma[k] = arma::dot(L.row(a).subvec(0, m), L.row(m).subvec(0, m));
//     }
//
//     value += arma::accu(arma::pow(Y - (X_tilde * sigma), 2));
//
//   }
//
//   for (R_xlen_t j = a + 1; j < q; j++) {
//
//     List Y_tilde_list_j = Y_tilde_list[j];
//
//     NumericVector Y_ = Y_tilde_list_j[a];
//     arma::vec Y(Y_.begin(), Y_.size(), false, true);
//
//     arma::vec sigma(K);
//     sigma.zeros();
//     for (R_xlen_t k = 0; k < K; k++) {
//       NumericMatrix L_ = L_list[k];
//       arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false) ;
//       sigma[k] = arma::dot(L.row(j).subvec(0, a), L.row(a).subvec(0, a));
//     }
//
//     value += arma::accu(arma::pow(Y - (X_tilde * sigma), 2));
//
//   }
//
//   return value;
//
// }

// // [[Rcpp::export]]
// double gradient(const List & Y_tilde_list, const arma::mat & X_tilde, const List & L_list, int a, int b, int z) {
//
//   double value = 0;
//
//   R_xlen_t K = X_tilde.n_cols;
//   R_xlen_t q = Y_tilde_list.size();
//
//   List Y_tilde_list_a = Y_tilde_list[a];
//   NumericVector Y_ = Y_tilde_list_a[a];
//   arma::vec Y(Y_.begin(), Y_.size(), false, true);
//
//   arma::vec sigma(K);
//   sigma.zeros();
//   for (R_xlen_t k = 0; k < K; k++) {
//     NumericMatrix L_ = L_list[k];
//     arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false) ;
//     sigma[k] = arma::dot(L.row(a).subvec(0, a), L.row(a).subvec(0, a));
//   }
//
//   arma::mat L = L_list[z];
//
//   value += -4 * L(a, b) * arma::dot(X_tilde.col(z), Y - (X_tilde * sigma));
//
//   for (R_xlen_t m = b; m <= a - 1; m++) {
//
//     NumericVector Y_ = Y_tilde_list_a[m];
//     arma::vec Y(Y_.begin(), Y_.size(), false, true);
//
//     arma::vec sigma(K);
//     sigma.zeros();
//     for (R_xlen_t k = 0; k < K; k++) {
//       NumericMatrix L_ = L_list[k];
//       arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false) ;
//       sigma[k] = arma::dot(L.row(a).subvec(0, m), L.row(m).subvec(0, m));
//     }
//
//     arma::mat L = L_list[z];
//
//     value += -2 * L(m, b) * arma::dot(X_tilde.col(z), Y - (X_tilde * sigma));
//
//   }
//
//   for (R_xlen_t j = a + 1; j < q; j++) {
//
//     List Y_tilde_list_j = Y_tilde_list[j];
//
//     NumericVector Y_ = Y_tilde_list_j[a];
//     arma::vec Y(Y_.begin(), Y_.size(), false, true);
//
//     arma::vec sigma(K);
//     sigma.zeros();
//     for (R_xlen_t k = 0; k < K; k++) {
//       NumericMatrix L_ = L_list[k];
//       arma::mat L(L_.begin(), L_.nrow(), L_.ncol(), false) ;
//       sigma[k] = arma::dot(L.row(j).subvec(0, a), L.row(a).subvec(0, a));
//     }
//
//     arma::mat L = L_list[z];
//
//     value += -2 * L(j, b) * arma::dot(X_tilde.col(z), Y - (X_tilde * sigma));
//
//   }
//
//   return value;
//
// }



