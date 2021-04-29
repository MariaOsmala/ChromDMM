#ifndef DMN_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define DMN_H

#include <RcppArmadillo.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>

using namespace Rcpp;

double neg_log_evidence_i(IntegerVector dataRow, NumericVector Lambda,
                          NumericVector LnGammaLambda0);

void disable_gsl_error_handler();


#endif