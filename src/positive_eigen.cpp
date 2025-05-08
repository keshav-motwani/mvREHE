// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include "Spectra/SymEigsSolver.h"
#include "Spectra/MatOp/DenseSymMatProd.h"

using namespace Rcpp;

// [[Rcpp::export]]
List positive_eigen(const arma::mat& A, double tol = 1e-10, double pos_tol = 1e-12, int buffer = 10) {

  using namespace Spectra;

  int n = A.n_rows;
  int max_k = n - 1;  // Spectra can compute up to n - 1
  int ncv = std::min(n, std::max(2 * max_k + 1, 20));

  Eigen::Map<const Eigen::MatrixXd> A_eigen(A.memptr(), n, n);
  DenseSymMatProd<double> op(A_eigen);
  SymEigsSolver<DenseSymMatProd<double>> eigs(op, max_k, ncv);
  eigs.init();

  eigs.compute(SortRule::LargestAlge, 1000, tol);

  if (eigs.info() != CompInfo::Successful) {
    stop("Spectra failed to converge.");
  }

  int nconv = eigs.eigenvalues().size();
  if (nconv == 0) {
    return List::create(
      Named("values") = arma::vec(),
      Named("vectors") = arma::mat(),
      Named("n_iter") = 0
    );
  }

  arma::vec ritz_vals(eigs.eigenvalues().data(), nconv);
  arma::mat ritz_vecs(eigs.eigenvectors().data(), n, nconv);

  std::vector<int> pos_indices;
  bool seen_positive = false;
  int buffer_counter = 0;

  for (int i = 0; i < nconv; ++i) {
    double val = ritz_vals[i];
    if (val > pos_tol) {
      seen_positive = true;
      buffer_counter = 0;
      pos_indices.push_back(i);
    } else if (seen_positive) {
      buffer_counter++;
      if (buffer_counter >= buffer) break;
    }
  }

  if (pos_indices.empty()) {
    return List::create(
      Named("values") = arma::vec(),
      Named("vectors") = arma::mat(),
      Named("n_iter") = 0
    );
  }

  arma::uvec pos_idx = arma::conv_to<arma::uvec>::from(pos_indices);
  arma::vec pos_vals = arma::vec(ritz_vals.elem(pos_idx));
  arma::mat pos_vecs = arma::mat(ritz_vecs.cols(pos_idx));
  arma::uvec sorted_idx = arma::sort_index(pos_vals, "descend");

  arma::vec sorted_vals = pos_vals(sorted_idx);
  arma::mat sorted_vecs = pos_vecs.cols(sorted_idx);

  return List::create(
    Named("values") = sorted_vals,
    Named("vectors") = sorted_vecs,
    Named("n_iter") = static_cast<unsigned int>(sorted_vals.n_elem)
  );

}
