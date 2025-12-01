// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// Soft-thresholding function
inline double soft_threshold(double z, double lambda) {
  if (z > lambda) return z - lambda;
  if (z < -lambda) return z + lambda;
  return 0.0;
}

// [[Rcpp::export]]
arma::mat lasso_from_cov(const arma::mat& Sigma_hat,
                                   const arma::vec& lambda_seq,
                                   const int max_iter = 1000,
                                   const double tolerance = 1e-6,
                                   const bool use_active_set = true,
                                   const double active_threshold = 1e-8) {
  
  // Extract XtX and XtY from Sigma_hat
  arma::mat XtX = Sigma_hat.submat(1, 1, Sigma_hat.n_rows - 1, Sigma_hat.n_cols - 1);
  arma::vec XtY = Sigma_hat.submat(1, 0, Sigma_hat.n_rows - 1, 0);
  
  int p = XtX.n_cols;
  int n_lambda = lambda_seq.n_elem;
  
  arma::mat beta_path(p, n_lambda, fill::zeros);
  arma::vec beta(p, fill::zeros);
  arma::vec beta_old(p);
  
  // Pre-compute diagonal elements
  arma::vec XtX_diag = XtX.diag();
  
  // Pre-compute residuals: r = XtY - XtX * beta
  arma::vec residuals = XtY;
  
  // Active set for tracking non-zero coefficients
  std::vector<bool> is_active(p, false);
  std::vector<int> active_indices;
  
  // Variables for convergence checking
  double max_change;
  int convergence_check_freq = std::min(10, p / 100 + 1); // Check every few iterations
  
  for (int lambda_idx = 0; lambda_idx < n_lambda; ++lambda_idx) {
    double lambda = lambda_seq[lambda_idx];
    
    // Warm start: use previous solution
    if (lambda_idx > 0) {
      beta = beta_path.col(lambda_idx - 1);
      residuals = XtY - XtX * beta;
    }
    
    // Update active set based on current beta
    if (use_active_set) {
      active_indices.clear();
      for (int j = 0; j < p; ++j) {
        is_active[j] = (std::abs(beta(j)) > active_threshold);
        if (is_active[j]) {
          active_indices.push_back(j);
        }
      }
    }
    
    for (int iter = 0; iter < max_iter; ++iter) {
      beta_old = beta;
      max_change = 0.0;
      
      // First pass: update active variables
      if (use_active_set && !active_indices.empty()) {
        for (int j : active_indices) {
          double beta_j_old = beta(j);
          double r_j = residuals(j) + XtX_diag(j) * beta_j_old;
          beta(j) = soft_threshold(r_j, lambda) / XtX_diag(j);
          
          if (beta(j) != beta_j_old) {
            double delta = beta(j) - beta_j_old;
            // Update residuals efficiently: r -= XtX(:,j) * delta
            residuals -= delta * XtX.col(j);
            max_change = std::max(max_change, std::abs(delta));
          }
        }
      }
      
      // Second pass: check inactive variables periodically or in final iterations
      bool check_inactive = !use_active_set || 
        (iter % 5 == 0) || 
        (iter > max_iter - 10);
      
      if (check_inactive) {
        for (int j = 0; j < p; ++j) {
          if (use_active_set && is_active[j]) continue; // Skip active variables
          
          double beta_j_old = beta(j);
          double r_j = residuals(j) + XtX_diag(j) * beta_j_old;
          beta(j) = soft_threshold(r_j, lambda) / XtX_diag(j);
          
          if (beta(j) != beta_j_old) {
            double delta = beta(j) - beta_j_old;
            residuals -= delta * XtX.col(j);
            max_change = std::max(max_change, std::abs(delta));
            
            // Add to active set if became non-zero
            if (use_active_set && !is_active[j] && std::abs(beta(j)) > active_threshold) {
              is_active[j] = true;
              active_indices.push_back(j);
            }
          }
        }
      }
      
      // Efficient convergence check
      if (iter % convergence_check_freq == 0) {
        if (max_change < tolerance) {
          break;
        }
      }
    }
    
    beta_path.col(lambda_idx) = beta;
  }
  
  return beta_path;
}


// // [[Rcpp::plugins(cpp11)]]
// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// 
// using namespace arma;
// using namespace Rcpp;
// 
// // Soft-thresholding function
// inline double soft_threshold(double z, double lambda) {
//   if (z > lambda) return z - lambda;
//   if (z < -lambda) return z + lambda;
//   return 0.0;
// }
// 
// // [[Rcpp::export]]
// arma::mat lasso_from_cov(const arma::mat& Sigma_hat,
//                          const arma::vec& lambda_seq,
//                          const int max_iter = 1000,
//                          const double tolerance = 1e-6) {
//   
//   // Extract XtX and XtY from Sigma_hat
//   arma::mat XtX = Sigma_hat.submat(1, 1, Sigma_hat.n_rows - 1, Sigma_hat.n_cols - 1);
//   arma::vec XtY = Sigma_hat.submat(1, 0, Sigma_hat.n_rows - 1, 0);
//   
//   int p = XtX.n_cols;
//   int n_lambda = lambda_seq.n_elem;
//   
//   arma::mat beta_path(p, n_lambda, fill::zeros);
//   arma::vec beta(p, fill::zeros);
//   arma::vec beta_old(p);
//   
//   arma::vec XtX_diag = XtX.diag();
//   
//   for (int lambda_idx = 0; lambda_idx < n_lambda; ++lambda_idx) {
//     double lambda = lambda_seq[lambda_idx];
//     
//     for (int iter = 0; iter < max_iter; ++iter) {
//       beta_old = beta;
//       
//       for (int j = 0; j < p; ++j) {
//         double r_j = XtY(j);
//         
//         if (j > 0) {
//           r_j -= dot(XtX.row(j).subvec(0, j - 1), beta.subvec(0, j - 1));
//         }
//         if (j < p - 1) {
//           r_j -= dot(XtX.row(j).subvec(j + 1, p - 1), beta.subvec(j + 1, p - 1));
//         }
//         
//         beta(j) = soft_threshold(r_j, lambda) / XtX_diag(j);
//       }
//       
//       if (max(abs(beta - beta_old)) < tolerance) {
//         break;
//       }
//     }
//     
//     beta_path.col(lambda_idx) = beta;
//   }
//   
//   return beta_path;
// }
