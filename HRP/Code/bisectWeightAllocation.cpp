#include <omp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

void bisectWeightAllocation(arma::mat& weightMat, arma::mat& covMat, arma::uword idx_start, arma::uword idx_end) {
  arma::colvec wi_upper;
  arma::colvec wi_lower;

  arma::mat temp_covMat_upper;
  arma::mat temp_covMat_lower;

  arma::uword idx_mid;

  double temp_scale_upper;
  double temp_scale_lower;

  if (idx_start != idx_end) {
    idx_mid = (idx_start + idx_end)/2;

    temp_covMat_upper = covMat.submat(idx_start, idx_start, idx_mid, idx_mid);
    temp_covMat_lower = covMat.submat(idx_mid+1, idx_mid+1, idx_end, idx_end);

    wi_upper = temp_covMat_upper.diag();
    wi_lower = temp_covMat_lower.diag();

    temp_scale_upper = as_scalar(wi_upper.t() * temp_covMat_upper * wi_upper);
    temp_scale_lower = as_scalar(wi_lower.t() * temp_covMat_lower * wi_lower);

    weightMat.submat(0, idx_start, 0, idx_mid) = weightMat.submat(0, idx_start, 0, idx_mid) * (temp_scale_lower /(temp_scale_upper + temp_scale_lower));
    weightMat.submat(0, idx_mid+1, 0, idx_end) = weightMat.submat(0, idx_mid+1, 0, idx_end) * (temp_scale_upper /(temp_scale_upper + temp_scale_lower));

#pragma omp task shared(weightMat, covMat) firstprivate(idx_start, idx_mid)
{
  bisectWeightAllocation(weightMat, covMat, idx_start, idx_mid);
}


#pragma omp task shared(weightMat, covMat) firstprivate(idx_mid, idx_end)
{
  bisectWeightAllocation(weightMat, covMat, idx_mid+1, idx_end);
}

  }
}

// [[Rcpp::export]]
arma::mat weightAllocation(NumericMatrix quasiDiagMatrix, arma::mat assetIndex) {
  int num_asset = quasiDiagMatrix.nrow();

  arma::mat covMat(quasiDiagMatrix.begin(), num_asset, num_asset, false);
  arma::mat weightMat_temp(1, num_asset, arma::fill::ones);
  arma::mat weightMat(1, num_asset, arma::fill::ones);

  omp_set_nested(0);

#pragma omp parallel
{
#pragma omp single
{
  bisectWeightAllocation(weightMat_temp, covMat, 0, num_asset-1);
}
}

for (int i = 0; i < num_asset; i++) {
  weightMat[0,i] = weightMat_temp[0,assetIndex[0,i]-1];
}

return weightMat;
}
