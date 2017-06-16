#include <omp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericMatrix distanceMatrix_elementwise(NumericMatrix MAT_CORR) {
  int i, j;
  NumericMatrix distanceMatrix(MAT_CORR.nrow(), MAT_CORR.nrow());

#pragma omp parallel for collapse(2)
  for (i = 0; i < MAT_CORR.nrow(); i++) {
    for (j = 0; j < MAT_CORR.ncol(); j++) {
      distanceMatrix(i,j) = std::pow(0.5*(1-MAT_CORR(i,j)), 0.5);  //pow is ^ in R. thus the expression says 0.5*(1-pij)^0.5)
    }
  }
  return distanceMatrix;
}


// [[Rcpp::export]]
NumericMatrix distanceMatrix_rowwise(NumericMatrix MAT_CORR) {
  int i,j,k;
  double temp_SUM = 0;
  NumericMatrix distanceMatrix(MAT_CORR.nrow(), MAT_CORR.nrow());

#pragma omp parallel for private(temp_SUM, j, k)
  for (i = 1; i < MAT_CORR.nrow(); ++i) {
    for (j = 0; j < i; ++j) {
      temp_SUM = 0;
      for (k = 0; k < MAT_CORR.nrow(); k++) {
        temp_SUM += std::pow(MAT_CORR(k,i) - MAT_CORR(k,j), 2);  //pow is ^ in R. thus the expression says 0.5*(1-pij)^2)
      }
      temp_SUM = std::pow(temp_SUM, 0.5);
      distanceMatrix(i,j) = temp_SUM;
      distanceMatrix(j,i) = temp_SUM;
    }
  }
  return distanceMatrix;
}
