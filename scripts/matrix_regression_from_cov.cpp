// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat compute_A1B_arma(const arma::mat& beta2, const arma::mat& B) {
  int p = beta2.n_rows;
  int rank = beta2.n_cols;
  int numCols = B.n_cols;
  arma::mat result(1 + p * rank, numCols, fill::zeros);
  
  // First row is B[0,]
  result.row(0) = B.row(0);
  
  // Compute the rest using the kronecker structure
  for(int r = 0; r < rank; r++) {
    for(int i = 0; i < p; i++) {
      for(int col = 0; col < numCols; col++) {
        double sum = 0.0;
        for(int k = 0; k < p; k++) {
          sum += beta2(k, r) * B(1 + i*p + k, col);
        }
        result(1 + r*p + i, col) = sum;
      }
    }
  }
  
  return result;
}

// [[Rcpp::export]]
arma::mat compute_A2B_arma(const arma::mat& beta1, const arma::mat& B) {
  int p = beta1.n_rows;
  int rank = beta1.n_cols;
  int numCols = B.n_cols;
  arma::mat result(1 + p * rank, numCols, fill::zeros);
  
  // First row is B[0,]
  result.row(0) = B.row(0);
  
  // Compute the rest using the kronecker structure
  for(int r = 0; r < rank; r++) {
    for(int i = 0; i < p; i++) {
      for(int col = 0; col < numCols; col++) {
        double sum = 0.0;
        for(int k = 0; k < p; k++) {
          sum += beta1(k, r) * B(1 + i*p + k, col);
        }
        result(1 + i*rank + r, col) = sum;
      }
    }
  }
  
  return result;
}

// // [[Rcpp::export]]
// List matrix_regression_from_cov(const arma::mat& Sigma_hat,
//                                     Nullable<arma::mat> V_nullable = R_NilValue,
//                                     int rank = 1,
//                                     double lambda = 0.1,
//                                     int max_iter = 1000,
//                                     double tolerance = 1e-6,
//                                     int n_init = 5) {
  
//   // Handle nullable V
//   bool has_V = V_nullable.isNotNull();
//   arma::mat V;
//   if (has_V) {
//     V = as<arma::mat>(V_nullable);
//   }
  
//   // Determine p
//   int p;
//   if (has_V) {
//     p = std::sqrt(V.n_rows - 1);
//   } else {
//     p = std::sqrt(Sigma_hat.n_rows - 1);
//   }
  
//   double best_objective = datum::inf;
//   arma::mat best_beta1, best_beta2;
//   arma::vec best_difference, best_objective_vec;
//   int best_init = -1;
  
//   // Multiple initializations
//   for (int init = 0; init < n_init; init++) {
//     // Initialize beta1 and beta2 with random normal values
//     arma::mat beta1 = randn(p, rank);
//     arma::mat beta2 = beta1;
    
//     arma::vec objective(max_iter, fill::zeros);
//     arma::vec difference(max_iter, fill::zeros);
    
//     int iter;
//     for (iter = 0; iter < max_iter; iter++) {
//       arma::mat beta1_old = beta1;
//       arma::mat beta2_old = beta2;
      
//       // Update beta1
//       arma::mat A, Sigma, temp;
//       if (has_V) {
//         A = compute_A1B_arma(beta2, V);
//         Sigma = A * Sigma_hat * A.t();
//       } else {
//         temp = compute_A1B_arma(beta2, Sigma_hat);
//         Sigma = compute_A1B_arma(beta2, temp.t()).t();
//       }
      
//       // Extract submatrices: Sigma[-1,-1] and Sigma[1,-1]
//       arma::mat Sigma_sub = Sigma.submat(1, 1, Sigma.n_rows-1, Sigma.n_cols-1);
//       arma::vec Sigma_row = Sigma.submat(0, 1, 0, Sigma.n_cols-1).t();
      
//       // Add regularization
//       Sigma_sub.diag() += lambda;
      
//       // Solve for beta1
//       arma::vec beta1_vec = arma::solve(Sigma_sub, Sigma_row, arma::solve_opts::fast);
//       beta1 = reshape(beta1_vec, p, rank);
      
//       // Update beta2
//       if (has_V) {
//         A = compute_A2B_arma(beta1, V);
//         Sigma = A * Sigma_hat * A.t();
//       } else {
//         temp = compute_A2B_arma(beta1, Sigma_hat);
//         Sigma = compute_A2B_arma(beta1, temp.t()).t();
//       }
      
//       // Extract submatrices again
//       Sigma_sub = Sigma.submat(1, 1, Sigma.n_rows-1, Sigma.n_cols-1);
//       Sigma_row = Sigma.submat(0, 1, 0, Sigma.n_cols-1).t();
      
//       // Add regularization
//       Sigma_sub.diag() += lambda;
      
//       // Solve for beta2 and reshape by row
//       arma::vec beta2_vec = arma::solve(Sigma_sub, Sigma_row, arma::solve_opts::fast);
//       // Reshape by row (transpose, reshape, transpose back)
//       arma::mat beta2_temp = reshape(beta2_vec, rank, p).t();
//       beta2 = beta2_temp;
      
//       // Compute final A matrix for objective
//       arma::mat beta_prod = beta1 * beta2.t();
//       arma::mat A_final(2, 1 + p * p, fill::zeros);
//       A_final(0, 0) = 1.0;
//       A_final(1, 0) = 0.0;
      
//       // Vectorize beta_prod and place it in the second row
//       arma::vec beta_prod_vec = vectorise(beta_prod);
//       for (int i = 0; i < beta_prod_vec.n_elem; i++) {
//         A_final(1, i + 1) = beta_prod_vec(i);
//       }
      
//       if (has_V) {
//         A_final = A_final * V;
//       }
      
//       Sigma = A_final * Sigma_hat * A_final.t();
      
//       // Compute objective
//       objective(iter) = Sigma(0, 0) - 2 * Sigma(0, 1) + Sigma(1, 1) + 
//         lambda * (accu(beta1 % beta1) + accu(beta2 % beta2));
      
      
//       // Compute difference for convergence check
//       arma::mat diff_mat = beta1 * beta2.t() - beta1_old * beta2_old.t();
//       double norm_diff = norm(diff_mat, "fro");
//       double norm_curr = norm(beta1 * beta2.t(), "fro");
//       difference(iter) = norm_diff / (norm_curr + 1e-12);
      
//       // Check convergence
//       if (iter > 4 && difference(iter) < tolerance) {
//         break;
//       }
//     }
    
//     // Warning if didn't converge
//     if (iter == max_iter - 1) {
//       Rcpp::warning("Did not converge for initialization %d", init + 1);
//     }
    
//     double final_objective = objective(iter - 1);
    
//     // Check if this is the best solution so far
//     if (final_objective < best_objective) {
//       // Rcpp::Rcout << "Init " << init + 1 << ": " << final_objective << " < " << best_objective << std::endl;
//       best_objective = final_objective;
//       best_beta1 = beta1;
//       best_beta2 = beta2;
//       best_init = init + 1;
      
//       // Store only non-zero elements
//       int actual_iters = iter;
//       best_difference = difference.head(actual_iters);
//       best_objective_vec = objective.head(actual_iters);
//     }
//   }
  
//   // Return best solution
//   return List::create(
//     Named("beta1") = best_beta1,
//     Named("beta2") = best_beta2,
//     Named("difference") = best_difference,
//     Named("objective") = best_objective_vec,
//     Named("init_number") = best_init,
//     Named("final_objective") = best_objective
//   );
// }

// [[Rcpp::export]]
List matrix_regression_from_cov(const arma::mat& Sigma_hat,
                                    Nullable<arma::mat> V_nullable = R_NilValue,
                                    int rank = 1,
                                    double lambda = 0.1,
                                    int max_iter = 1000,
                                    double tolerance = 1e-6,
                                    int n_init = 1,
                                    Nullable<arma::mat> beta1_init_nullable = R_NilValue,
                                    Nullable<arma::mat> beta2_init_nullable = R_NilValue) {

  // Handle nullable V
  bool has_V = V_nullable.isNotNull();
  arma::mat V;
  if (has_V) {
    V = as<arma::mat>(V_nullable);
  }

  // Determine p
  int p;
  if (has_V) {
    p = std::round(std::sqrt((double)(V.n_rows - 1)));
  } else {
    p = std::round(std::sqrt((double)(Sigma_hat.n_rows - 1)));
  }

  // Initialize beta1 and beta2
  arma::mat beta1(p, rank, arma::fill::zeros);
  arma::mat beta2(p, rank, arma::fill::zeros);

  bool has_init = beta1_init_nullable.isNotNull() && beta2_init_nullable.isNotNull();
  if (has_init) {
    // Warm start: copy provided columns (may be fewer than rank)
    arma::mat b1_in = as<arma::mat>(beta1_init_nullable);
    arma::mat b2_in = as<arma::mat>(beta2_init_nullable);
    int prev_rank = std::min((int)b1_in.n_cols, rank);
    beta1.cols(0, prev_rank - 1) = b1_in.cols(0, prev_rank - 1);
    beta2.cols(0, prev_rank - 1) = b2_in.cols(0, prev_rank - 1);
  } else {
    // =====================================================================
    // DETERMINISTIC SVD INITIALIZATION
    // =====================================================================

    // Compute effective joint covariance matrix
    arma::mat Sigma_eff;
    if (has_V) {
      Sigma_eff = V * Sigma_hat * V.t();
    } else {
      Sigma_eff = Sigma_hat;
    }

    // Extract X covariance (Sigma_XX) and cross-covariance (Sigma_XY)
    // Index 0 is Y, Indices 1 to p^2 are vec(X)
    int p2 = p * p;
    arma::mat Sigma_XX = Sigma_eff.submat(1, 1, p2, p2);
    arma::vec Sigma_XY = Sigma_eff.submat(1, 0, p2, 0);

    // Solve unconstrained OLS in basis space
    Sigma_XX.diag() += lambda;
    arma::vec B_vec_0 = arma::solve(Sigma_XX, Sigma_XY, arma::solve_opts::fast);

    // Reshape to p x p matrix and force exact symmetry
    arma::mat B_mat_0 = arma::reshape(B_vec_0, p, p);
    arma::mat B_sym = 0.5 * (B_mat_0 + B_mat_0.t());

    // Symmetric Eigenvalue Decomposition
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, B_sym);

    // Sort eigenvalues by absolute magnitude (descending)
    arma::uvec sort_idx = arma::sort_index(arma::abs(eigval), "descend");
    int actual_rank = std::min(rank, p);
    arma::uvec top_idx = sort_idx.head(actual_rank);

    for(int r = 0; r < actual_rank; r++) {
      double ev = eigval(top_idx(r));
      double sqrt_ev = std::sqrt(std::abs(ev));
      beta1.col(r) = eigvec.col(top_idx(r)) * sqrt_ev;
      beta2.col(r) = eigvec.col(top_idx(r)) * (ev >= 0 ? sqrt_ev : -sqrt_ev);
    }
  }
  
  // =====================================================================
  // 2. ALS OPTIMIZATION (Runs ONLY ONCE now)
  // =====================================================================
  
  arma::vec objective(max_iter, fill::zeros);
  arma::vec difference(max_iter, fill::zeros);
  
  int iter;
  for (iter = 0; iter < max_iter; iter++) {
    arma::mat beta1_old = beta1;
    arma::mat beta2_old = beta2;
    
    // Update beta1
    arma::mat A, Sigma, temp;
    if (has_V) {
      A = compute_A1B_arma(beta2, V);
      Sigma = A * Sigma_hat * A.t();
    } else {
      temp = compute_A1B_arma(beta2, Sigma_hat);
      Sigma = compute_A1B_arma(beta2, temp.t()).t();
    }
    
    arma::mat Sigma_sub = Sigma.submat(1, 1, Sigma.n_rows-1, Sigma.n_cols-1);
    arma::vec Sigma_row = Sigma.submat(0, 1, 0, Sigma.n_cols-1).t();
    Sigma_sub.diag() += lambda;
    
    arma::vec beta1_vec = arma::solve(Sigma_sub, Sigma_row, arma::solve_opts::fast);
    beta1 = reshape(beta1_vec, p, rank);
    
    // Update beta2
    if (has_V) {
      A = compute_A2B_arma(beta1, V);
      Sigma = A * Sigma_hat * A.t();
    } else {
      temp = compute_A2B_arma(beta1, Sigma_hat);
      Sigma = compute_A2B_arma(beta1, temp.t()).t();
    }
    
    Sigma_sub = Sigma.submat(1, 1, Sigma.n_rows-1, Sigma.n_cols-1);
    Sigma_row = Sigma.submat(0, 1, 0, Sigma.n_cols-1).t();
    Sigma_sub.diag() += lambda;
    
    arma::vec beta2_vec = arma::solve(Sigma_sub, Sigma_row, arma::solve_opts::fast);
    arma::mat beta2_temp = reshape(beta2_vec, rank, p).t();
    beta2 = beta2_temp;
    
    // Compute final A matrix for objective
    arma::mat beta_prod = beta1 * beta2.t();
    arma::mat A_final(2, 1 + p * p, fill::zeros);
    A_final(0, 0) = 1.0;
    
    arma::vec beta_prod_vec = vectorise(beta_prod);
    for (int i = 0; i < beta_prod_vec.n_elem; i++) {
      A_final(1, i + 1) = beta_prod_vec(i);
    }
    
    if (has_V) {
      A_final = A_final * V;
    }
    
    Sigma = A_final * Sigma_hat * A_final.t();
    
    // Compute objective
    objective(iter) = Sigma(0, 0) - 2 * Sigma(0, 1) + Sigma(1, 1) + 
      lambda * (accu(beta1 % beta1) + accu(beta2 % beta2));
    
    // Compute difference for convergence check
    arma::mat diff_mat = beta1 * beta2.t() - beta1_old * beta2_old.t();
    double norm_diff = norm(diff_mat, "fro");
    double norm_curr = norm(beta1 * beta2.t(), "fro");
    difference(iter) = norm_diff / (norm_curr + 1e-12);
    
    // Check convergence
    if (iter > 4 && difference(iter) < tolerance) {
      break;
    }
  }
  
  // Warning if didn't converge
  if (iter == max_iter) {
    Rcpp::warning("ALS did not converge within max_iter");
    iter--; // Adjust for indexing
  }
  
  int actual_iters = iter + 1;
  arma::vec best_difference = difference.head(actual_iters);
  arma::vec best_objective_vec = objective.head(actual_iters);
  double final_objective = best_objective_vec(actual_iters - 1);
  
  // Return best (and only) solution
  return List::create(
    Named("beta1") = beta1,
    Named("beta2") = beta2,
    Named("difference") = best_difference,
    Named("objective") = best_objective_vec,
    Named("init_number") = 1,
    Named("final_objective") = final_objective
  );
}