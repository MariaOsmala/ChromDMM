#ifndef REGUL_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define REGUL_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

double REG_H_LAMBDA(NumericVector lambda);
  
double REG_DERIV_LAMBDA_G(NumericVector lambda, int index);

#endif