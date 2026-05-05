// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
double frobenius_inner_product(S4 A, S4 B, Rcpp::Nullable<Rcpp::S4> C = R_NilValue) {

  IntegerVector A_p = A.slot("p");
  IntegerVector A_i = A.slot("i");
  NumericVector A_x = A.slot("x");

  IntegerVector B_p = B.slot("p");
  IntegerVector B_i = B.slot("i");
  NumericVector B_x = B.slot("x");

  int ncol = A_p.size() - 1;
  double result = 0.0;

  bool has_C = C.isNotNull();

  if (!has_C) {
    // ---------------------------------------------------
    // Original 2-Matrix Logic: <A, B>
    // ---------------------------------------------------
    for (int j = 0; j < ncol; j++) {
      int A_start = A_p[j], A_end = A_p[j + 1];
      int B_start = B_p[j], B_end = B_p[j + 1];

      int i_A = A_start, i_B = B_start;

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
  } else {
    // ---------------------------------------------------
    // 3-Matrix Logic: <A, B \circ C>
    // ---------------------------------------------------
    S4 C_obj(C.get());
    IntegerVector C_p = C_obj.slot("p");
    IntegerVector C_i = C_obj.slot("i");
    NumericVector C_x = C_obj.slot("x");

    for (int j = 0; j < ncol; j++) {
      int A_start = A_p[j], A_end = A_p[j + 1];
      int B_start = B_p[j], B_end = B_p[j + 1];
      int C_start = C_p[j], C_end = C_p[j + 1];

      int i_A = A_start, i_B = B_start, i_C = C_start;

      while (i_A < A_end && i_B < B_end && i_C < C_end) {
        int rA = A_i[i_A];
        int rB = B_i[i_B];
        int rC = C_i[i_C];

        // If all three matrices have a non-zero entry at this row/col
        if (rA == rB && rB == rC) {
          double prod = A_x[i_A] * B_x[i_B] * C_x[i_C];
          if (rA == j) {
            result += prod;
          } else {
            result += 2.0 * prod;
          }
          i_A++; i_B++; i_C++;
        } else {
          // Find the maximum row index among the three pointers
          int max_r = rA;
          if (rB > max_r) max_r = rB;
          if (rC > max_r) max_r = rC;

          // Increment any pointer that is lagging behind the maximum
          if (rA < max_r) i_A++;
          if (rB < max_r) i_B++;
          if (rC < max_r) i_C++;
        }
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
arma::mat compute_WDY(S4 D, S4 W, const arma::mat& Y) {

  IntegerVector D_p = D.slot("p");
  IntegerVector D_i = D.slot("i");
  NumericVector D_x = D.slot("x");

  IntegerVector W_p = W.slot("p");
  IntegerVector W_i = W.slot("i");
  NumericVector W_x = W.slot("x");

  int N = Y.n_rows;
  int q = Y.n_cols;
  int ncol = D_p.size() - 1;

  // Initialize the output matrix Z = (W * D) %*% Y
  arma::mat Z(N, q, arma::fill::zeros);

  for (int l = 0; l < ncol; l++) {
    int D_start = D_p[l], D_end = D_p[l + 1];
    int W_start = W_p[l], W_end = W_p[l + 1];

    int idx_D = D_start, idx_W = W_start;

    // Find matching (row, col) entries in the sparse structures
    while (idx_D < D_end && idx_W < W_end) {
      int rD = D_i[idx_D];
      int rW = W_i[idx_W];

      if (rD == rW) {
        double val = D_x[idx_D] * W_x[idx_W];
        int r = rD; // Row index
        int c = l;  // Col index

        if (r == c) {
          // Diagonal elements: Add once
          for (int m = 0; m < q; m++) {
            Z(r, m) += val * Y(c, m);
          }
        } else {
          // Off-diagonal elements: Mirror for symmetry
          for (int m = 0; m < q; m++) {
            Z(r, m) += val * Y(c, m);
            Z(c, m) += val * Y(r, m);
          }
        }

        idx_D++;
        idx_W++;
      } else if (rD < rW) {
        idx_D++;
      } else {
        idx_W++;
      }
    }
  }

  return Z;
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

// [[Rcpp::export]]
S4 compute_W_sparse(S4 R, S4 R2, NumericVector R_diag, NumericVector R2_diag,
                    double h2, double rho_g, double rho_e, double N, double M_snps) {

  IntegerVector R_p = R.slot("p");
  IntegerVector R_i = R.slot("i");
  NumericVector R_x = R.slot("x");

  IntegerVector R2_p = R2.slot("p");
  IntegerVector R2_i = R2.slot("i");
  NumericVector R2_x = R2.slot("x");

  int M_cols = R_p.size() - 1;

  std::vector<int> W_p(M_cols + 1, 0);
  std::vector<int> W_i;
  std::vector<double> W_x;

  int estimated_nnz = std::max(R_i.size(), R2_i.size());
  W_i.reserve(estimated_nnz);
  W_x.reserve(estimated_nnz);

  for (int c = 0; c < M_cols; c++) {
    int ptr_R = R_p[c], end_R = R_p[c + 1];
    int ptr_R2 = R2_p[c], end_R2 = R2_p[c + 1];

    while (ptr_R < end_R || ptr_R2 < end_R2) {

      int r;
      double val_R = 0.0, val_R2 = 0.0;

      bool use_R = (ptr_R < end_R);
      bool use_R2 = (ptr_R2 < end_R2);

      if (use_R && use_R2) {
        if (R_i[ptr_R] == R2_i[ptr_R2]) {
          r = R_i[ptr_R];
          val_R = R_x[ptr_R++];
          val_R2 = R2_x[ptr_R2++];
        } else if (R_i[ptr_R] < R2_i[ptr_R2]) {
          r = R_i[ptr_R];
          val_R = R_x[ptr_R++];
        } else {
          r = R2_i[ptr_R2];
          val_R2 = R2_x[ptr_R2++];
        }
      } else if (use_R) {
        r = R_i[ptr_R];
        val_R = R_x[ptr_R++];
      } else {
        r = R2_i[ptr_R2];
        val_R2 = R2_x[ptr_R2++];
      }

      double term_s = (N / M_snps) * R2_diag[r] * h2 + R_diag[r] * (1.0 - h2);
      double term_t = (N / M_snps) * R2_diag[c] * h2 + R_diag[c] * (1.0 - h2);
      double term_st = (N / M_snps) * val_R2 * rho_g * h2 + val_R * rho_e * (1.0 - h2);

      double denom = (R2_diag[r] * R2_diag[c]) * (term_s * term_t + term_st * term_st);

      if (denom > 0) {
        W_i.push_back(r);
        W_x.push_back(1.0 / denom);
      }
    }

    W_p[c + 1] = W_i.size();
  }

  S4 W = clone(R);
  W.slot("p") = IntegerVector(W_p.begin(), W_p.end());
  W.slot("i") = IntegerVector(W_i.begin(), W_i.end());
  W.slot("x") = NumericVector(W_x.begin(), W_x.end());

  return W;
}
