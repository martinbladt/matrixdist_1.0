#include <Rcpp.h>
using namespace Rcpp;

NumericVector mWeibullden(NumericVector x, NumericVector pi, NumericMatrix T, double beta);

NumericVector mWeibullcdf(NumericVector x, NumericVector pi, NumericMatrix T, double beta, bool lower_tail = true);