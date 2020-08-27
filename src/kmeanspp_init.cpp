#include <Rcpp.h>
using namespace Rcpp;

double sqdist(NumericVector ind1, NumericVector ind2) {
  int n = ind1.size();
  double out = 0.0;

  for(int i = 0; i < n; ++i) {
    out += pow(ind1[i] - ind2[i], 2.0);
  }
  return out;

}
//nearest_center( mat(noncenters[j], _), mat, centers )
double nearest_center(NumericVector x, NumericMatrix mat, IntegerVector ind) {
  int rows = mat.nrow(); //this is not used?
  double out = 1e20;

  for (int i = 0; i < ind.size(); ++i) {    //loop through all centers
    double tmp = sqdist(x, mat(ind[i], _)); //distance(Euclidean) between a sample and the center
    out = tmp < out ? tmp : out;
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector kmeanspp_initialize(NumericMatrix mat, int K) { //mat is NxS matrix
  RNGScope scope;  	// ensure RNG gets set/reset

  IntegerVector noncenters = seq(0, mat.nrow()-1); // 0:999
  int S = mat.ncol();                             // number of cols/bins
  Function sample("sample");  //Using the Function class, you can call R functions from Rcpp

  //take one center c_1, chosen uniformly at random from all samples
  IntegerVector centers = sample(noncenters, 1); //sample one from noncenters
  noncenters = setdiff(noncenters, centers);

  for (int k = 1; k < K; k++) {
    //probabilities of choosing the next center
    NumericVector probs(noncenters.size());
    for (int j = 0; j < noncenters.size(); j++) {

      //compute the shortest distance from a data point to the closest center already chosen
      probs[j] = nearest_center(mat(noncenters[j], _), mat, centers);
    }
    int newcenter = as<int>(sample(noncenters, 1, false, probs)); //replace is false
    centers.push_back(newcenter);
    noncenters = setdiff(noncenters, centers);
  }
  return centers;
}