#include <omp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat flatCluster(int index, int threshold, NumericMatrix clusterMatrix) {
  arma::mat temp_index;
  arma::mat temp_index_a;
  arma::mat temp_index_b;

  if(index <= threshold) {
    temp_index.set_size(1,1);
    temp_index(0,0) = index;
    return temp_index;
  }

  temp_index_a = flatCluster(clusterMatrix(index - threshold - 1,1), threshold, clusterMatrix);
  temp_index_b = flatCluster(clusterMatrix(index - threshold - 1,0), threshold, clusterMatrix);
  temp_index = join_rows(temp_index_a, temp_index_b);

  return temp_index;
}

// [[Rcpp::export]]
arma::mat clusterIndex(NumericMatrix clusterMatrix, NumericMatrix MAT_COV) {

  int num_asset;
  int nrow_clusterMatrix;

  nrow_clusterMatrix = clusterMatrix.nrow();
  num_asset = MAT_COV.nrow();

  arma::mat assetIndex;

  assetIndex = join_rows(flatCluster(clusterMatrix(nrow_clusterMatrix - 1,1), num_asset, clusterMatrix),
                         flatCluster(clusterMatrix(nrow_clusterMatrix - 1,0), num_asset, clusterMatrix));

  return assetIndex;
}

// [[Rcpp::export]]
NumericMatrix quasiDiag(NumericMatrix MAT_COV, arma::mat assetIndex) {

  int num_asset;
  int index_asset;
  num_asset = MAT_COV.nrow();

  NumericMatrix interMatrix(num_asset, num_asset);
  NumericMatrix quasiDiagMatrix(num_asset, num_asset);

#pragma omp parallel for private(index_asset)
  for(int i = 0; i < num_asset; ++i) {
    index_asset = assetIndex(0,i) - 1;
    // printf("current %d, total %d\n", num_asset-1, index_asset);
    interMatrix(_, i) = MAT_COV(_, index_asset);
  }

#pragma omp parallel for private(index_asset)
  for(int i = 0; i < num_asset; ++i) {
    index_asset = assetIndex(0,i) - 1;
    quasiDiagMatrix(i, _) = interMatrix(index_asset, _);
  }

  return quasiDiagMatrix;
}
