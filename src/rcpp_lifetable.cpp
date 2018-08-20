#include <Rcpp.h>
#include "fundem/lifetable.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector first_moment_survival(
        NumericVector mx, NumericVector ax, NumericVector nx)
{
    std::vector<double> mxv(mx.begin(), mx.end());
    std::vector<double> axv(ax.begin(), ax.end());
    std::vector<double> nxv(nx.begin(), nx.end());

    std::vector<double> survival(mx.size());
    fundem::FirstMomentSurvival(
            &mxv[0], &axv[0], &nxv[0], &survival[0], mx.size());

    return NumericVector(survival.begin(), survival.end());
}
