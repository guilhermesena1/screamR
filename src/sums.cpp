#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::Map;

//' Sum columns of a sparse matrix
//'
//' @export
// [[Rcpp::export]]
NumericVector ColSums(Eigen::SparseMatrix<double> m){
  int ncol = m.cols();
  NumericVector ans(ncol);
  for (int j = 0; j < ncol; ++j) {
    ans[j] = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator
        i_(m,j); i_; ++i_) {
      ans[j] += i_.value();
    }
  }
  return ans;
}

//' Sum rows of a sparse matrix
//'
//' @export
// [[Rcpp::export]]
NumericVector RowSums(Eigen::SparseMatrix<double> m){
  return ColSums(Eigen::SparseMatrix<double>(m.transpose()));
}

//' Row variances of sparse matrix
//'
//' @export
// [[Rcpp::export]]
NumericVector ColVars(Eigen::SparseMatrix<double> m){
  int ncol = m.cols();
  int nrow = m.rows();
  NumericVector col_means = ColSums(m) / static_cast<double>(nrow);
  NumericVector ans(ncol);

  double z;
  for (int j = 0; j < ncol; ++j) {
    size_t num_zero = nrow;
    for (Eigen::SparseMatrix<double>::InnerIterator
        i_(m,j); i_; ++i_) {
      --num_zero;
      z = (i_.value() - col_means[j]);
      ans[j] += z*z;
    }

    // sum contribution of zeros to variance
    ans[j] += num_zero * col_means[j] * col_means[j];
  }

  return ans / static_cast<double>(nrow - 1);
}




