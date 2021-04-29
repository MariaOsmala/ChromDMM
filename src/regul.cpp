// computes the regularization term using lambda, alpha= exp(lambda)
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

double REG_H_LAMBDA(Rcpp::NumericVector lambda) {
  int i;
  double hk = 0.0;
  
  //from second element to the last
  for (i = 1; i < lambda.length(); i++) {
    double diff = exp(lambda[i]) - exp(lambda[i-1]);
    hk += diff * diff;
  }
  return hk;
}



//computes the value of function g for j=index, input lambda values
double REG_DERIV_LAMBDA_G(Rcpp::NumericVector lambda, int index) {
  int prev = index - 1, next = index + 1;
  if (index == 0) prev = 0;
  if (index == lambda.length()-1) next = index;
  return 2*(2*exp(lambda[index]) - exp(lambda[prev]) - exp(lambda[next]));
}


