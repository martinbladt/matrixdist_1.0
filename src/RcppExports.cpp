// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// embeddedMC
NumericMatrix embeddedMC(NumericMatrix T, NumericVector t);
RcppExport SEXP _matrixdist_embeddedMC(SEXP TSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type T(TSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(embeddedMC(T, t));
    return rcpp_result_gen;
END_RCPP
}
// cumulateMatrix
NumericMatrix cumulateMatrix(NumericMatrix A);
RcppExport SEXP _matrixdist_cumulateMatrix(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(cumulateMatrix(A));
    return rcpp_result_gen;
END_RCPP
}
// cumulateVector
NumericVector cumulateVector(NumericVector A);
RcppExport SEXP _matrixdist_cumulateVector(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(cumulateVector(A));
    return rcpp_result_gen;
END_RCPP
}
// initialState
long initialState(NumericVector cumulatedPi, double u);
RcppExport SEXP _matrixdist_initialState(SEXP cumulatedPiSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type cumulatedPi(cumulatedPiSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(initialState(cumulatedPi, u));
    return rcpp_result_gen;
END_RCPP
}
// newState
long newState(long previousState, NumericMatrix cumulatedEmbeddedMC, double u);
RcppExport SEXP _matrixdist_newState(SEXP previousStateSEXP, SEXP cumulatedEmbeddedMCSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< long >::type previousState(previousStateSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cumulatedEmbeddedMC(cumulatedEmbeddedMCSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(newState(previousState, cumulatedEmbeddedMC, u));
    return rcpp_result_gen;
END_RCPP
}
// rphasetype
NumericVector rphasetype(int n, NumericVector pi, NumericMatrix T, NumericVector t);
RcppExport SEXP _matrixdist_rphasetype(SEXP nSEXP, SEXP piSEXP, SEXP TSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type T(TSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(rphasetype(n, pi, T, t));
    return rcpp_result_gen;
END_RCPP
}
// matrix_product
NumericMatrix matrix_product(NumericMatrix a, NumericMatrix b);
RcppExport SEXP _matrixdist_matrix_product(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_product(a, b));
    return rcpp_result_gen;
END_RCPP
}
// matrix_sum
NumericMatrix matrix_sum(const NumericMatrix& A, const NumericMatrix& B);
RcppExport SEXP _matrixdist_matrix_sum(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_sum(A, B));
    return rcpp_result_gen;
END_RCPP
}
// LInf_norm
double LInf_norm(const NumericMatrix& A);
RcppExport SEXP _matrixdist_LInf_norm(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(LInf_norm(A));
    return rcpp_result_gen;
END_RCPP
}
// solve_linear_system
NumericMatrix solve_linear_system(NumericMatrix A, const NumericMatrix& B);
RcppExport SEXP _matrixdist_solve_linear_system(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_linear_system(A, B));
    return rcpp_result_gen;
END_RCPP
}
// matrix_inverse
NumericMatrix matrix_inverse(NumericMatrix A);
RcppExport SEXP _matrixdist_matrix_inverse(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_inverse(A));
    return rcpp_result_gen;
END_RCPP
}
// matrix_exponential
NumericMatrix matrix_exponential(const NumericMatrix& A);
RcppExport SEXP _matrixdist_matrix_exponential(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_exponential(A));
    return rcpp_result_gen;
END_RCPP
}
// matrixMax
double matrixMax(const NumericMatrix& A);
RcppExport SEXP _matrixdist_matrixMax(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matrixMax(A));
    return rcpp_result_gen;
END_RCPP
}
// matrixMaxDiagonal
double matrixMaxDiagonal(const NumericMatrix& A);
RcppExport SEXP _matrixdist_matrixMaxDiagonal(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matrixMaxDiagonal(A));
    return rcpp_result_gen;
END_RCPP
}
// matrix_power
NumericMatrix matrix_power(int n, const NumericMatrix& A);
RcppExport SEXP _matrixdist_matrix_power(SEXP nSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_power(n, A));
    return rcpp_result_gen;
END_RCPP
}
// phdensity
NumericVector phdensity(NumericVector x, NumericVector pi, NumericMatrix T);
RcppExport SEXP _matrixdist_phdensity(SEXP xSEXP, SEXP piSEXP, SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi(piSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(phdensity(x, pi, T));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _matrixdist_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_matrixdist_embeddedMC", (DL_FUNC) &_matrixdist_embeddedMC, 2},
    {"_matrixdist_cumulateMatrix", (DL_FUNC) &_matrixdist_cumulateMatrix, 1},
    {"_matrixdist_cumulateVector", (DL_FUNC) &_matrixdist_cumulateVector, 1},
    {"_matrixdist_initialState", (DL_FUNC) &_matrixdist_initialState, 2},
    {"_matrixdist_newState", (DL_FUNC) &_matrixdist_newState, 3},
    {"_matrixdist_rphasetype", (DL_FUNC) &_matrixdist_rphasetype, 4},
    {"_matrixdist_matrix_product", (DL_FUNC) &_matrixdist_matrix_product, 2},
    {"_matrixdist_matrix_sum", (DL_FUNC) &_matrixdist_matrix_sum, 2},
    {"_matrixdist_LInf_norm", (DL_FUNC) &_matrixdist_LInf_norm, 1},
    {"_matrixdist_solve_linear_system", (DL_FUNC) &_matrixdist_solve_linear_system, 2},
    {"_matrixdist_matrix_inverse", (DL_FUNC) &_matrixdist_matrix_inverse, 1},
    {"_matrixdist_matrix_exponential", (DL_FUNC) &_matrixdist_matrix_exponential, 1},
    {"_matrixdist_matrixMax", (DL_FUNC) &_matrixdist_matrixMax, 1},
    {"_matrixdist_matrixMaxDiagonal", (DL_FUNC) &_matrixdist_matrixMaxDiagonal, 1},
    {"_matrixdist_matrix_power", (DL_FUNC) &_matrixdist_matrix_power, 2},
    {"_matrixdist_phdensity", (DL_FUNC) &_matrixdist_phdensity, 3},
    {"_matrixdist_rcpp_hello_world", (DL_FUNC) &_matrixdist_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_matrixdist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
