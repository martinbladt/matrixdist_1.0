#include <Rcpp.h>
using namespace Rcpp;

NumericVector mweibullden(NumericVector x, NumericVector pi, NumericMatrix T, double beta);

NumericVector mweibullcdf(NumericVector x, NumericVector pi, NumericMatrix T, double beta, bool lower_tail = true);