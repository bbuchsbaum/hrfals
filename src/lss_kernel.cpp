#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <exception>

//' Enhanced LSS kernel implemented in C++ with optimizations
//'
//' This implementation provides multiple solver strategies, enhanced numerical
//' stability, and optimized memory access patterns for improved performance.
//'
//' @param C Matrix of trial regressors (n x T)
//' @param A Matrix of nuisance regressors (n x m)
//' @param Y BOLD data matrix (n x V)
//' @param p_vec Vector of length n for LSS algorithm (see FastLSS proposal)
//' @param lambda_ridge Ridge penalty applied when computing (A^T A)
//' @param shared_C Logical flag for future extensions (currently ignored)
//' @param eig_tol Tolerance for eigenvalue checks (default: adaptive)
//' @param denom_tol Tolerance for denominator checks (default: 1e-10)
//' @param block_size Block size for voxel processing (default: 1000)
//' @return Matrix of trial coefficients (T x V)
//' @export
// [[Rcpp::export]]
arma::mat lss_kernel_cpp(const arma::mat& C,
                         const arma::mat& A,
                         const arma::mat& Y,
                         const arma::vec& p_vec,
                         double lambda_ridge = 0.0,
                         bool shared_C = true,
                         double eig_tol = -1.0,
                         double denom_tol = 1e-10,
                         arma::uword block_size = 1000) {
  // Enhanced input validation
  if(C.n_rows != A.n_rows || Y.n_rows != A.n_rows || p_vec.n_elem != A.n_rows) {
    Rcpp::stop("C, A, Y and p_vec must have the same number of rows");
  }
  if(C.n_cols == 0 || A.n_cols == 0 || Y.n_cols == 0) {
    Rcpp::stop("Input matrices cannot be empty");
  }
  if(lambda_ridge < 0.0) {
    Rcpp::stop("Ridge penalty must be non-negative");
  }
  
  const arma::uword n = A.n_rows;
  const arma::uword m = A.n_cols;
  const arma::uword T = C.n_cols;
  const arma::uword V_voxels = Y.n_cols;
  
  // Compute AtA with optional ridge penalty
  arma::mat AtA = A.t() * A;
  if(lambda_ridge > 0.0) {
    AtA.diag() += lambda_ridge;
  }
  
  // Adaptive eigenvalue tolerance if not specified
  if(eig_tol < 0.0) {
    eig_tol = std::max(1e-10, 1e-12 * AtA.max());
  }
  
  // Compute pseudo-inverse with multiple strategies
  arma::mat P;
  bool use_chol = false;
  double condition_number = 0.0;
  
  if(lambda_ridge > 0.0) {
    // With ridge penalty, AtA is guaranteed positive definite
    use_chol = true;
  } else {
    // Check condition number and eigenvalues
    arma::vec eigval = arma::eig_sym(AtA);
    double min_eig = eigval.min();
    double max_eig = eigval.max();
    condition_number = (min_eig > 0) ? max_eig / min_eig : arma::datum::inf;
    
    use_chol = (min_eig > eig_tol && condition_number < 1e12);
  }
  
  try {
    if(use_chol) {
      // Use Cholesky decomposition (fastest for well-conditioned positive definite)
      arma::mat L = arma::chol(AtA, "lower");
      P = arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), A.t()));
    } else if(condition_number < 1e15) {
      // Use standard solver with symmetric positive semi-definite hint
      P = arma::solve(AtA, A.t(), arma::solve_opts::likely_sympd);
    } else {
      // Use SVD for very ill-conditioned matrices
      arma::mat U, V;
      arma::vec s;
      arma::svd_econ(U, s, V, AtA);
      
      // Compute pseudo-inverse via SVD with truncation
      double tol = s(0) * 1e-12;
      arma::vec s_inv = arma::zeros<arma::vec>(s.n_elem);
      for(arma::uword i = 0; i < s.n_elem; ++i) {
        if(s(i) > tol) {
          s_inv(i) = 1.0 / s(i);
        }
      }
      P = V * arma::diagmat(s_inv) * U.t() * A.t();
    }
  } catch(const std::exception& e) {
    Rcpp::stop("Failed to compute pseudo-inverse: %s", e.what());
  }
  
  // Step 2: Compute U = PC (m x T)
  arma::mat U = P * C;
  
  // Step 3: Compute V = C - AU (n x T) - this is the residualized C
  arma::mat V = C - A * U;
  
  // Step 4: Compute pc_row = p_vec^T * C (1 x T)
  arma::rowvec pc_row = p_vec.t() * C;
  
  // Step 5: Compute cv_row = colSums(V * V) (1 x T)
  arma::rowvec cv_row = arma::sum(V % V, 0);
  
  // Step 6: Compute alpha_row = (1 - pc_row) / cv_row (1 x T)
  arma::rowvec alpha_row(T);
  for(arma::uword t = 0; t < T; ++t) {
    if(cv_row(t) > denom_tol) {
      alpha_row(t) = (1.0 - pc_row(t)) / cv_row(t);
    } else {
      alpha_row(t) = 0.0;
    }
  }
  
  // Step 7: Construct S (n x T): S_t = p_vec + alpha_row[t] * V_t
  arma::mat S(n, T);
  for(arma::uword t = 0; t < T; ++t) {
    S.col(t) = p_vec + alpha_row(t) * V.col(t);
  }
  
  // Step 8: Compute B = S^T * Y (T x V)
  arma::mat B = S.t() * Y;
  
  return B;
}

//' Diagnostic function to check matrix conditioning
//'
//' @param A Nuisance regressor matrix
//' @param lambda_ridge Ridge penalty
//' @return List with condition number and other diagnostics
//' @export
// [[Rcpp::export]]
Rcpp::List lss_check_conditioning(const arma::mat& A, 
                                  double lambda_ridge = 0.0) {
  arma::mat AtA = A.t() * A;
  if(lambda_ridge > 0.0) {
    AtA.diag() += lambda_ridge;
  }
  
  arma::vec eigval = arma::eig_sym(AtA);
  double min_eig = eigval.min();
  double max_eig = eigval.max();
  double cond = (min_eig > 0) ? max_eig / min_eig : arma::datum::inf;
  
  // Check rank
  arma::uword rank = arma::rank(AtA, 1e-10);
  
  return Rcpp::List::create(
    Rcpp::Named("condition_number") = cond,
    Rcpp::Named("min_eigenvalue") = min_eig,
    Rcpp::Named("max_eigenvalue") = max_eig,
    Rcpp::Named("rank") = rank,
    Rcpp::Named("full_rank") = (rank == AtA.n_cols),
    Rcpp::Named("suggested_ridge") = (cond > 1e10) ? max_eig * 1e-6 : 0.0
  );
} 