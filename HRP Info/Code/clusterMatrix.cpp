#include <omp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericMatrix clusterMatrix(NumericMatrix MAT_CORR) {

  //Auxillary Variables
  int i = 0;
  int dim = MAT_CORR.nrow();
  double max_MAT_CORR = max(MAT_CORR);

  arma::mat    temp_MAT_CORR(MAT_CORR.begin(), dim, dim, false);

  arma::uword temp_idx_row = 0;
  arma::uword temp_idx_col = 0;

  double min_dist = 0.0;

  arma::mat    temp_cluster_mat;
  arma::colvec temp_cluster_vec;
  arma::rowvec temp_cluster_rvec;

  //result matrix
  NumericMatrix clusterMatrix(dim-1, 4);

  arma::colvec clusterIndex((dim-1)*2);

  //fill diagonal of corr matrix
#pragma omp parallel for
  for(i = 0; i < dim; ++i) {
    temp_MAT_CORR(i,i) = max_MAT_CORR*2;
  }

#pragma omp parallel for
  for(i = 0; i < (dim-1)*2; ++i) {
    clusterIndex(i) = i+1;
  }


  for(i = 0; i < dim-1; i++) {
    //calculate clustermatrix row
    min_dist = temp_MAT_CORR.min(temp_idx_row, temp_idx_col);

    clusterMatrix(i,0) = clusterIndex(temp_idx_row);
    clusterMatrix(i,1) = clusterIndex(temp_idx_col);
    clusterMatrix(i,2) = min_dist;
    clusterMatrix(i,3) =
      (clusterMatrix(i,0) <= dim ? 1 : 0) + (clusterMatrix(i,1) <= dim ? 1 : 0);


    //re-construct correlation matrix
    clusterIndex.shed_row(temp_idx_row);
    clusterIndex.shed_row(temp_idx_col);

    temp_cluster_mat = join_rows(temp_MAT_CORR.col(temp_idx_row),
                                 temp_MAT_CORR.col(temp_idx_col));

    temp_cluster_vec = min(temp_cluster_mat,1);
    temp_cluster_rvec = temp_cluster_vec.t();
    temp_cluster_rvec.insert_cols(temp_cluster_vec.n_elem, 1);
    temp_cluster_rvec(temp_cluster_rvec.n_elem-1) = max_MAT_CORR;

    temp_MAT_CORR = join_rows(temp_MAT_CORR, temp_cluster_vec);
    temp_MAT_CORR = join_cols(temp_MAT_CORR, temp_cluster_rvec);

    temp_MAT_CORR.shed_row(temp_idx_row);
    temp_MAT_CORR.shed_row(temp_idx_col);
    temp_MAT_CORR.shed_col(temp_idx_row);
    temp_MAT_CORR.shed_col(temp_idx_col);

  }
  return clusterMatrix;
}
