#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

//' Core LSS kernel implemented in C++
//'
//' @param C Matrix of trial regressors (n x T)
//' @param A Matrix of nuisance regressors (n x m)
//' @param Y BOLD data matrix (n x V)
//' @param lambda_ridge Ridge penalty applied when computing (A^T A)
//' @param shared_C Logical flag for future extensions (currently ignored)
//' @return Matrix of trial coefficients (T x V)
//' @export
// [[Rcpp::export]]
arma::mat lss_kernel_cpp(const arma::mat& C,
                         const arma::mat& A,
                         const arma::mat& Y,
                         double lambda_ridge = 0.0,
                         bool shared_C = true) {
  if(C.n_rows != A.n_rows || Y.n_rows != A.n_rows) {
    Rcpp::stop("C, A and Y must have the same number of rows");
  }
  arma::mat AtA = A.t() * A;
  if(lambda_ridge != 0.0) {
    AtA.diag() += lambda_ridge;
  }
  arma::mat P = arma::solve(AtA, A.t());
  arma::mat Cres = C - A * (P * C);
  arma::mat Yres = Y - A * (P * Y);

  arma::uword Tt = Cres.n_cols;
  arma::uword V  = Yres.n_cols;
  arma::mat B(Tt, V, arma::fill::zeros);

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(arma::uword t = 0; t < Tt; ++t) {
    double denom = arma::dot(Cres.col(t), Cres.col(t));
    if(denom != 0.0) {
      arma::rowvec num = Cres.col(t).t() * Yres;
      B.row(t) = num / denom;
    }
  }

  return B;
}
