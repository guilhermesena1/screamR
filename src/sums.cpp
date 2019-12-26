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

