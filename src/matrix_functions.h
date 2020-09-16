#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix matrix_product(NumericMatrix a, NumericMatrix b);

NumericMatrix matrix_inverse(NumericMatrix A);

NumericMatrix matrix_exponential(const NumericMatrix & A);

NumericMatrix matrix_power(int n, const NumericMatrix & A);