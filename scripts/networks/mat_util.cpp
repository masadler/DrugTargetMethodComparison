#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd stoch_col_norm_(const Eigen::MatrixXd& W)
{
    Eigen::MatrixXd res(W.rows(), W.cols());
    Eigen::VectorXd colsums      = W.colwise().sum();
    const double    empt_col_val = 1.0 / W.cols();
    const double    zero_col     = 0.00001;
    for (unsigned int i = 0; i < W.cols(); ++i)
    {
        if ((W.col(i)).sum() <= zero_col)
            res.col(i).fill(empt_col_val);
        else
            res.col(i) = W.col(i) / colsums(i);
    }

    return res;
}

// [[Rcpp::export]]
NumericMatrix add_weights_adjacency(NumericMatrix A, IntegerVector idx1, IntegerVector idx2, NumericVector weight)
{
    for (int i=0; i<idx1.size(); i++){
        A(idx1[i], idx2[i]) = weight[i];
        A(idx2[i], idx1[i]) = weight[i];
    }
    return A;
}