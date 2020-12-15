#include <Rcpp.h>
using namespace Rcpp;

NumericVector mweibullden(NumericVector x, NumericVector alpha, NumericMatrix S, double beta);

NumericVector mweibullcdf(NumericVector x, NumericVector alpha, NumericMatrix S, double beta, bool lower_tail = true);