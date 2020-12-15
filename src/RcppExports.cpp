// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// default_step_length
double default_step_length(const NumericMatrix& S);
RcppExport SEXP _matrixdist_default_step_length(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(default_step_length(S));
    return rcpp_result_gen;
END_RCPP
}
// runge_kutta
void runge_kutta(NumericMatrix& avector, NumericMatrix& bvector, NumericMatrix& cmatrix, double dt, double h, const NumericMatrix& S, const NumericMatrix& t);
RcppExport SEXP _matrixdist_runge_kutta(SEXP avectorSEXP, SEXP bvectorSEXP, SEXP cmatrixSEXP, SEXP dtSEXP, SEXP hSEXP, SEXP SSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type avector(avectorSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type bvector(bvectorSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type cmatrix(cmatrixSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type t(tSEXP);
    runge_kutta(avector, bvector, cmatrix, dt, h, S, t);
    return R_NilValue;
END_RCPP
}
// EMstep_RK
void EMstep_RK(double h, NumericVector& alpha, NumericMatrix& S, const NumericVector& obs, const NumericVector& weight, const NumericVector& rcens, const NumericVector& rcweight);
RcppExport SEXP _matrixdist_EMstep_RK(SEXP hSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP obsSEXP, SEXP weightSEXP, SEXP rcensSEXP, SEXP rcweightSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcens(rcensSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcweight(rcweightSEXP);
    EMstep_RK(h, alpha, S, obs, weight, rcens, rcweight);
    return R_NilValue;
END_RCPP
}
// a_rungekutta
void a_rungekutta(NumericMatrix& avector, double dt, double h, const NumericMatrix& S);
RcppExport SEXP _matrixdist_a_rungekutta(SEXP avectorSEXP, SEXP dtSEXP, SEXP hSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type avector(avectorSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type S(SSEXP);
    a_rungekutta(avector, dt, h, S);
    return R_NilValue;
END_RCPP
}
// logLikelihoodPH_RK
double logLikelihoodPH_RK(double h, NumericVector& alpha, NumericMatrix& S, const NumericVector& obs, const NumericVector& weight, const NumericVector& rcens, const NumericVector& rcweight);
RcppExport SEXP _matrixdist_logLikelihoodPH_RK(SEXP hSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP obsSEXP, SEXP weightSEXP, SEXP rcensSEXP, SEXP rcweightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcens(rcensSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcweight(rcweightSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikelihoodPH_RK(h, alpha, S, obs, weight, rcens, rcweight));
    return rcpp_result_gen;
END_RCPP
}
// logLikelihoodMweibull_RK
double logLikelihoodMweibull_RK(double h, NumericVector& alpha, NumericMatrix& S, double beta, const NumericVector& obs, const NumericVector& weight, const NumericVector& rcens, const NumericVector& rcweight);
RcppExport SEXP _matrixdist_logLikelihoodMweibull_RK(SEXP hSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP, SEXP obsSEXP, SEXP weightSEXP, SEXP rcensSEXP, SEXP rcweightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcens(rcensSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcweight(rcweightSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikelihoodMweibull_RK(h, alpha, S, beta, obs, weight, rcens, rcweight));
    return rcpp_result_gen;
END_RCPP
}
// logLikelihoodMpareto_RK
double logLikelihoodMpareto_RK(double h, NumericVector& alpha, NumericMatrix& S, double beta, const NumericVector& obs, const NumericVector& weight, const NumericVector& rcens, const NumericVector& rcweight);
RcppExport SEXP _matrixdist_logLikelihoodMpareto_RK(SEXP hSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP, SEXP obsSEXP, SEXP weightSEXP, SEXP rcensSEXP, SEXP rcweightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcens(rcensSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcweight(rcweightSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikelihoodMpareto_RK(h, alpha, S, beta, obs, weight, rcens, rcweight));
    return rcpp_result_gen;
END_RCPP
}
// logLikelihoodMlognormal_RK
double logLikelihoodMlognormal_RK(double h, NumericVector& alpha, NumericMatrix& S, double beta, const NumericVector& obs, const NumericVector& weight, const NumericVector& rcens, const NumericVector& rcweight);
RcppExport SEXP _matrixdist_logLikelihoodMlognormal_RK(SEXP hSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP, SEXP obsSEXP, SEXP weightSEXP, SEXP rcensSEXP, SEXP rcweightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcens(rcensSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcweight(rcweightSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikelihoodMlognormal_RK(h, alpha, S, beta, obs, weight, rcens, rcweight));
    return rcpp_result_gen;
END_RCPP
}
// logLikelihoodMloglogistic_RK
double logLikelihoodMloglogistic_RK(double h, NumericVector& alpha, NumericMatrix& S, NumericVector beta, const NumericVector& obs, const NumericVector& weight, const NumericVector& rcens, const NumericVector& rcweight);
RcppExport SEXP _matrixdist_logLikelihoodMloglogistic_RK(SEXP hSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP, SEXP obsSEXP, SEXP weightSEXP, SEXP rcensSEXP, SEXP rcweightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcens(rcensSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcweight(rcweightSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikelihoodMloglogistic_RK(h, alpha, S, beta, obs, weight, rcens, rcweight));
    return rcpp_result_gen;
END_RCPP
}
// logLikelihoodMgompertz_RK
double logLikelihoodMgompertz_RK(double h, NumericVector& alpha, NumericMatrix& S, double beta, const NumericVector& obs, const NumericVector& weight, const NumericVector& rcens, const NumericVector& rcweight);
RcppExport SEXP _matrixdist_logLikelihoodMgompertz_RK(SEXP hSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP, SEXP obsSEXP, SEXP weightSEXP, SEXP rcensSEXP, SEXP rcweightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcens(rcensSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcweight(rcweightSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikelihoodMgompertz_RK(h, alpha, S, beta, obs, weight, rcens, rcweight));
    return rcpp_result_gen;
END_RCPP
}
// logLikelihoodMgev_RK
double logLikelihoodMgev_RK(double h, NumericVector& alpha, NumericMatrix& S, NumericVector beta, const NumericVector& obs, const NumericVector& weight, const NumericVector& rcens, const NumericVector& rcweight);
RcppExport SEXP _matrixdist_logLikelihoodMgev_RK(SEXP hSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP, SEXP obsSEXP, SEXP weightSEXP, SEXP rcensSEXP, SEXP rcweightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcens(rcensSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcweight(rcweightSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikelihoodMgev_RK(h, alpha, S, beta, obs, weight, rcens, rcweight));
    return rcpp_result_gen;
END_RCPP
}
// reversTransformData
List reversTransformData(const NumericVector& observations, const NumericVector& weights, const NumericVector& beta);
RcppExport SEXP _matrixdist_reversTransformData(SEXP observationsSEXP, SEXP weightsSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(reversTransformData(observations, weights, beta));
    return rcpp_result_gen;
END_RCPP
}
// derivativeMatrixweibull
double derivativeMatrixweibull(double h, const NumericVector& obs, const NumericVector& weight, const NumericVector& rcens, const NumericVector& rcweight, NumericVector& alpha, NumericMatrix& S, double beta);
RcppExport SEXP _matrixdist_derivativeMatrixweibull(SEXP hSEXP, SEXP obsSEXP, SEXP weightSEXP, SEXP rcensSEXP, SEXP rcweightSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcens(rcensSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rcweight(rcweightSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(derivativeMatrixweibull(h, obs, weight, rcens, rcweight, alpha, S, beta));
    return rcpp_result_gen;
END_RCPP
}
// embeddedMC
NumericMatrix embeddedMC(NumericMatrix S);
RcppExport SEXP _matrixdist_embeddedMC(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(embeddedMC(S));
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
NumericVector rphasetype(int n, NumericVector alpha, NumericMatrix S);
RcppExport SEXP _matrixdist_rphasetype(SEXP nSEXP, SEXP alphaSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(rphasetype(n, alpha, S));
    return rcpp_result_gen;
END_RCPP
}
// riph
NumericVector riph(int n, String dist_type, NumericVector alpha, NumericMatrix S, NumericVector beta);
RcppExport SEXP _matrixdist_riph(SEXP nSEXP, SEXP dist_typeSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< String >::type dist_type(dist_typeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(riph(n, dist_type, alpha, S, beta));
    return rcpp_result_gen;
END_RCPP
}
// rmatrixgev
NumericVector rmatrixgev(int n, NumericVector alpha, NumericMatrix S, double mu, double sigma, double xi);
RcppExport SEXP _matrixdist_rmatrixgev(SEXP nSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP xiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    rcpp_result_gen = Rcpp::wrap(rmatrixgev(n, alpha, S, mu, sigma, xi));
    return rcpp_result_gen;
END_RCPP
}
// phdensity
NumericVector phdensity(NumericVector x, NumericVector alpha, NumericMatrix S);
RcppExport SEXP _matrixdist_phdensity(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(phdensity(x, alpha, S));
    return rcpp_result_gen;
END_RCPP
}
// phcdf
NumericVector phcdf(NumericVector x, NumericVector alpha, NumericMatrix S, bool lower_tail);
RcppExport SEXP _matrixdist_phcdf(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(phcdf(x, alpha, S, lower_tail));
    return rcpp_result_gen;
END_RCPP
}
// mweibullden
NumericVector mweibullden(NumericVector x, NumericVector alpha, NumericMatrix S, double beta);
RcppExport SEXP _matrixdist_mweibullden(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(mweibullden(x, alpha, S, beta));
    return rcpp_result_gen;
END_RCPP
}
// mweibullcdf
NumericVector mweibullcdf(NumericVector x, NumericVector alpha, NumericMatrix S, double beta, bool lower_tail);
RcppExport SEXP _matrixdist_mweibullcdf(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(mweibullcdf(x, alpha, S, beta, lower_tail));
    return rcpp_result_gen;
END_RCPP
}
// RunFunction
NumericVector RunFunction(NumericVector a, Function func);
RcppExport SEXP _matrixdist_RunFunction(SEXP aSEXP, SEXP funcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< Function >::type func(funcSEXP);
    rcpp_result_gen = Rcpp::wrap(RunFunction(a, func));
    return rcpp_result_gen;
END_RCPP
}
// mparetoden
NumericVector mparetoden(NumericVector x, NumericVector alpha, NumericMatrix S, double beta);
RcppExport SEXP _matrixdist_mparetoden(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(mparetoden(x, alpha, S, beta));
    return rcpp_result_gen;
END_RCPP
}
// mparetocdf
NumericVector mparetocdf(NumericVector x, NumericVector alpha, NumericMatrix S, double beta, bool lower_tail);
RcppExport SEXP _matrixdist_mparetocdf(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(mparetocdf(x, alpha, S, beta, lower_tail));
    return rcpp_result_gen;
END_RCPP
}
// mlognormalden
NumericVector mlognormalden(NumericVector x, NumericVector alpha, NumericMatrix S, double beta);
RcppExport SEXP _matrixdist_mlognormalden(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(mlognormalden(x, alpha, S, beta));
    return rcpp_result_gen;
END_RCPP
}
// mlognormalcdf
NumericVector mlognormalcdf(NumericVector x, NumericVector alpha, NumericMatrix S, double beta, bool lower_tail);
RcppExport SEXP _matrixdist_mlognormalcdf(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(mlognormalcdf(x, alpha, S, beta, lower_tail));
    return rcpp_result_gen;
END_RCPP
}
// mloglogisticden
NumericVector mloglogisticden(NumericVector x, NumericVector alpha, NumericMatrix S, NumericVector beta);
RcppExport SEXP _matrixdist_mloglogisticden(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(mloglogisticden(x, alpha, S, beta));
    return rcpp_result_gen;
END_RCPP
}
// mloglogisticcdf
NumericVector mloglogisticcdf(NumericVector x, NumericVector alpha, NumericMatrix S, NumericVector beta, bool lower_tail);
RcppExport SEXP _matrixdist_mloglogisticcdf(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(mloglogisticcdf(x, alpha, S, beta, lower_tail));
    return rcpp_result_gen;
END_RCPP
}
// mgompertzden
NumericVector mgompertzden(NumericVector x, NumericVector alpha, NumericMatrix S, double beta);
RcppExport SEXP _matrixdist_mgompertzden(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(mgompertzden(x, alpha, S, beta));
    return rcpp_result_gen;
END_RCPP
}
// mgompertzcdf
NumericVector mgompertzcdf(NumericVector x, NumericVector alpha, NumericMatrix S, double beta, bool lower_tail);
RcppExport SEXP _matrixdist_mgompertzcdf(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP betaSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(mgompertzcdf(x, alpha, S, beta, lower_tail));
    return rcpp_result_gen;
END_RCPP
}
// mgevden
NumericVector mgevden(NumericVector x, NumericVector alpha, NumericMatrix S, double mu, double sigma, double xi);
RcppExport SEXP _matrixdist_mgevden(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP xiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    rcpp_result_gen = Rcpp::wrap(mgevden(x, alpha, S, mu, sigma, xi));
    return rcpp_result_gen;
END_RCPP
}
// mgevcdf
NumericVector mgevcdf(NumericVector x, NumericVector alpha, NumericMatrix S, double mu, double sigma, double xi, bool lower_tail);
RcppExport SEXP _matrixdist_mgevcdf(SEXP xSEXP, SEXP alphaSEXP, SEXP SSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP xiSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(mgevcdf(x, alpha, S, mu, sigma, xi, lower_tail));
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
NumericMatrix solve_linear_system(NumericMatrix A1, const NumericMatrix& B);
RcppExport SEXP _matrixdist_solve_linear_system(SEXP A1SEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_linear_system(A1, B));
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
// clone_vector
NumericVector clone_vector(NumericVector v);
RcppExport SEXP _matrixdist_clone_vector(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(clone_vector(v));
    return rcpp_result_gen;
END_RCPP
}
// clone_matrix
NumericMatrix clone_matrix(NumericMatrix m);
RcppExport SEXP _matrixdist_clone_matrix(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(clone_matrix(m));
    return rcpp_result_gen;
END_RCPP
}
// matrix_VanLoan
NumericMatrix matrix_VanLoan(const NumericMatrix& A1, const NumericMatrix& A2, const NumericMatrix& B1);
RcppExport SEXP _matrixdist_matrix_VanLoan(SEXP A1SEXP, SEXP A2SEXP, SEXP B1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A2(A2SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type B1(B1SEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_VanLoan(A1, A2, B1));
    return rcpp_result_gen;
END_RCPP
}
// diagonal_vector
NumericMatrix diagonal_vector(const NumericVector& vec);
RcppExport SEXP _matrixdist_diagonal_vector(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(diagonal_vector(vec));
    return rcpp_result_gen;
END_RCPP
}
// sumPH
List sumPH(NumericVector alpha1, NumericMatrix S1, NumericVector alpha2, NumericMatrix S2);
RcppExport SEXP _matrixdist_sumPH(SEXP alpha1SEXP, SEXP S1SEXP, SEXP alpha2SEXP, SEXP S2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type alpha1(alpha1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type S2(S2SEXP);
    rcpp_result_gen = Rcpp::wrap(sumPH(alpha1, S1, alpha2, S2));
    return rcpp_result_gen;
END_RCPP
}
// random_structure
List random_structure(int p, String structure, double scale_factor);
RcppExport SEXP _matrixdist_random_structure(SEXP pSEXP, SEXP structureSEXP, SEXP scale_factorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< String >::type structure(structureSEXP);
    Rcpp::traits::input_parameter< double >::type scale_factor(scale_factorSEXP);
    rcpp_result_gen = Rcpp::wrap(random_structure(p, structure, scale_factor));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_matrixdist_default_step_length", (DL_FUNC) &_matrixdist_default_step_length, 1},
    {"_matrixdist_runge_kutta", (DL_FUNC) &_matrixdist_runge_kutta, 7},
    {"_matrixdist_EMstep_RK", (DL_FUNC) &_matrixdist_EMstep_RK, 7},
    {"_matrixdist_a_rungekutta", (DL_FUNC) &_matrixdist_a_rungekutta, 4},
    {"_matrixdist_logLikelihoodPH_RK", (DL_FUNC) &_matrixdist_logLikelihoodPH_RK, 7},
    {"_matrixdist_logLikelihoodMweibull_RK", (DL_FUNC) &_matrixdist_logLikelihoodMweibull_RK, 8},
    {"_matrixdist_logLikelihoodMpareto_RK", (DL_FUNC) &_matrixdist_logLikelihoodMpareto_RK, 8},
    {"_matrixdist_logLikelihoodMlognormal_RK", (DL_FUNC) &_matrixdist_logLikelihoodMlognormal_RK, 8},
    {"_matrixdist_logLikelihoodMloglogistic_RK", (DL_FUNC) &_matrixdist_logLikelihoodMloglogistic_RK, 8},
    {"_matrixdist_logLikelihoodMgompertz_RK", (DL_FUNC) &_matrixdist_logLikelihoodMgompertz_RK, 8},
    {"_matrixdist_logLikelihoodMgev_RK", (DL_FUNC) &_matrixdist_logLikelihoodMgev_RK, 8},
    {"_matrixdist_reversTransformData", (DL_FUNC) &_matrixdist_reversTransformData, 3},
    {"_matrixdist_derivativeMatrixweibull", (DL_FUNC) &_matrixdist_derivativeMatrixweibull, 8},
    {"_matrixdist_embeddedMC", (DL_FUNC) &_matrixdist_embeddedMC, 1},
    {"_matrixdist_cumulateMatrix", (DL_FUNC) &_matrixdist_cumulateMatrix, 1},
    {"_matrixdist_cumulateVector", (DL_FUNC) &_matrixdist_cumulateVector, 1},
    {"_matrixdist_initialState", (DL_FUNC) &_matrixdist_initialState, 2},
    {"_matrixdist_newState", (DL_FUNC) &_matrixdist_newState, 3},
    {"_matrixdist_rphasetype", (DL_FUNC) &_matrixdist_rphasetype, 3},
    {"_matrixdist_riph", (DL_FUNC) &_matrixdist_riph, 5},
    {"_matrixdist_rmatrixgev", (DL_FUNC) &_matrixdist_rmatrixgev, 6},
    {"_matrixdist_phdensity", (DL_FUNC) &_matrixdist_phdensity, 3},
    {"_matrixdist_phcdf", (DL_FUNC) &_matrixdist_phcdf, 4},
    {"_matrixdist_mweibullden", (DL_FUNC) &_matrixdist_mweibullden, 4},
    {"_matrixdist_mweibullcdf", (DL_FUNC) &_matrixdist_mweibullcdf, 5},
    {"_matrixdist_RunFunction", (DL_FUNC) &_matrixdist_RunFunction, 2},
    {"_matrixdist_mparetoden", (DL_FUNC) &_matrixdist_mparetoden, 4},
    {"_matrixdist_mparetocdf", (DL_FUNC) &_matrixdist_mparetocdf, 5},
    {"_matrixdist_mlognormalden", (DL_FUNC) &_matrixdist_mlognormalden, 4},
    {"_matrixdist_mlognormalcdf", (DL_FUNC) &_matrixdist_mlognormalcdf, 5},
    {"_matrixdist_mloglogisticden", (DL_FUNC) &_matrixdist_mloglogisticden, 4},
    {"_matrixdist_mloglogisticcdf", (DL_FUNC) &_matrixdist_mloglogisticcdf, 5},
    {"_matrixdist_mgompertzden", (DL_FUNC) &_matrixdist_mgompertzden, 4},
    {"_matrixdist_mgompertzcdf", (DL_FUNC) &_matrixdist_mgompertzcdf, 5},
    {"_matrixdist_mgevden", (DL_FUNC) &_matrixdist_mgevden, 6},
    {"_matrixdist_mgevcdf", (DL_FUNC) &_matrixdist_mgevcdf, 7},
    {"_matrixdist_matrix_product", (DL_FUNC) &_matrixdist_matrix_product, 2},
    {"_matrixdist_matrix_sum", (DL_FUNC) &_matrixdist_matrix_sum, 2},
    {"_matrixdist_LInf_norm", (DL_FUNC) &_matrixdist_LInf_norm, 1},
    {"_matrixdist_solve_linear_system", (DL_FUNC) &_matrixdist_solve_linear_system, 2},
    {"_matrixdist_matrix_inverse", (DL_FUNC) &_matrixdist_matrix_inverse, 1},
    {"_matrixdist_matrix_exponential", (DL_FUNC) &_matrixdist_matrix_exponential, 1},
    {"_matrixdist_matrixMax", (DL_FUNC) &_matrixdist_matrixMax, 1},
    {"_matrixdist_matrixMaxDiagonal", (DL_FUNC) &_matrixdist_matrixMaxDiagonal, 1},
    {"_matrixdist_matrix_power", (DL_FUNC) &_matrixdist_matrix_power, 2},
    {"_matrixdist_clone_vector", (DL_FUNC) &_matrixdist_clone_vector, 1},
    {"_matrixdist_clone_matrix", (DL_FUNC) &_matrixdist_clone_matrix, 1},
    {"_matrixdist_matrix_VanLoan", (DL_FUNC) &_matrixdist_matrix_VanLoan, 3},
    {"_matrixdist_diagonal_vector", (DL_FUNC) &_matrixdist_diagonal_vector, 1},
    {"_matrixdist_sumPH", (DL_FUNC) &_matrixdist_sumPH, 4},
    {"_matrixdist_random_structure", (DL_FUNC) &_matrixdist_random_structure, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_matrixdist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
