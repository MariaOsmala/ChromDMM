#include <RcppArmadillo.h> // we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>

#include "regul.h"
#include "dmn.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]


#define BIG_DBL 1.0e9

/*The function below computes the value of the lower bound to be optimized
 * i.e. the whole expected log posterior Q with shifting and flipping
 * Computes the term only depending on \alpha_{k}^m or \lambda_l^m
 * . Does not consider \pi. Computes the value for one cluster k and datatype m
 * lambda are the current values, lparams contains the current parameter values (E(z_i))
 * List of 8
 $ pi       : Z array(dim=c(S,2,N) )These are z_k for all i, S and 2 flip states
 $ data     : int [1:1000, 1:L_x 0 1 0 1 0 0 0 0 1 1 ...
 $ nu       : num 1
 $ etah     : num 1
 $ nuh      : num 10
 $ hkm      :List of 1001
 ..$ : num 0.000227
 ..$ : num 0.592
 $ hkm_index: int 42
 * 
 */

// [[Rcpp::export]]
double neg_log_evidence_lambda_pi_shift_flip(Rcpp::NumericVector lambda, Rcpp::List lparams)
{
  
  
  int i, s, f, j;
  
  Rcpp::IntegerMatrix aanX = as<Rcpp::IntegerMatrix>(lparams["data"]); // N x L_x
  arma::Cube<double> adPi = as<arma::cube>(lparams["pi"]); //size Sx2xN
  // adPi.n_cols 2 (flips)
  // adPi.n_rows S
  // adPi.n_slices N
  // arma::mat A = adPi.row(0); //2xN
  // arma::mat B = adPi.col(0); //SxN
  // arma::mat C = adPi.slice(0); //Sx2
  // adPi(s,f,i)
 
  // Rcpp::NumericVector adPi = as<Rcpp::NumericVector>(lparams["pi"]);   // size N, these are z_{km}
  double GAMMA_ITA = as<double>(lparams["eta"]);
  double GAMMA_NU = as<double>(lparams["nu"]);
  double GAMMA_ITA_H = as<double>(lparams["etah"]);
  double GAMMA_NU_H = as<double>(lparams["nuh"]);
  
  const int Lx = aanX.ncol(); // L_x
  const int N = aanX.nrow();
  const int S = adPi.n_rows;
  const int La = Lx + S -1;
  
  /*dLogE collects the terms \sum_{n=1}^N \sum_{j=1}^S E[z_i] \log \gamma (x_ij+ alpha_j)
   and \sum_{n=1}^N E[z_i]* lng( \sum_j(x_ij+alpha-j) )*/
  double dLogE = 0.0;
  

  double dLogEAlpha_final = 0.0;
  double dSumAlpha = 0.0; // sum of alpha \sum_{j=1}^S \alpha_{j}
  double dSumLambda = 0.0; // sum of lambda?
  double dHk = 0.0; //this is not used
  // double dWeight = 0.0; // \sum_{n=1}^N E(z_i)
  Rcpp::NumericVector dWeight(S);
  Rcpp::NumericVector dAlpha(lambda.length());
  Rcpp::NumericMatrix dWeight_f(S,N);
  // Rcpp::NumericVector adSumAlphaN(N); // \sum_{j}^S \alpha_{j}+x_{ij}
  Rcpp::NumericMatrix adSumAlphaN(S,N);
  Rcpp::NumericVector dSumAlpha_s(S);
  Rcpp::NumericVector dLogEAlpha(S); // \log \gamma ( \sum_{j=1}^S  \alpha_{j} )
  
  for(s = 0; s < S; s++){
    dSumAlpha_s[s]=0.0;
    double tmp=0;
    for (i = 0; i < N; i++) {
      adSumAlphaN[s,i]=0.0;
      double tmp2=0.0;
      for(f = 0; f < 2; f++){
          tmp += adPi(s,f,i); //sum over i and f
          tmp2 += adPi(s,f,i);
      }
      dWeight_f[s,i]=tmp2;
    }
    dWeight[s]=tmp;
    // Rcpp::Rcout << dWeight(s) << "";
  }
  


  for (j = 0; j < La; j++) {
    dAlpha[j] = exp(lambda[j]);
    dSumLambda += lambda[j];
    dSumAlpha += dAlpha[j];
  }
  
  
  
  for(s=0; s<S; s++){
 
    for (j = 0; j < Lx; j++) {
      
      dSumAlpha_s[s] += dAlpha[j+s];
      /* compute the logarithm of the Gamma function,
      dAlpha can not be zero or negative.
      Function computed using the real Lanczos method */
      
      // dLogEAlpha1_[s]
      dLogEAlpha[s] += gsl_sf_lngamma(dAlpha[j+s]); // \sum_{j=1}^S \log \gamma (\alpha_{j} )
     
      const double lngammaAlpha0_1 = gsl_sf_lngamma(dAlpha[j+s]); //\log \gamma (\alpha_{js} ) f=1
      const double lngammaAlpha0_2 = gsl_sf_lngamma(dAlpha[La-S+s-j]); //\log \gamma (\alpha_{js} ) f=2
      for (i = 0; i < N; i++) {
        const double dN = aanX(i, j); // x_{ij}
        const double dAlphaN_1 = dAlpha[j+s] + dN; // \alpha_{j}+x_{ij} f=1
        const double dAlphaN_2 = dAlpha[La-S+s-j] + dN; // f=2
        
        const double lngammaAlphaN_1 = dN ? gsl_sf_lngamma(dAlphaN_1) : lngammaAlpha0_1; //if dN exists or is non-zero, compute lnGamma(dAlphaN), else compute lnGamma(DaLpha)
        const double lngammaAlphaN_2 = dN ? gsl_sf_lngamma(dAlphaN_2) : lngammaAlpha0_2; //if dN exists or is non-zero, compute lnGamma(dAlphaN), else compute lnGamma(DaLpha)
        
        adSumAlphaN[s,i] += dAlphaN_1; // \sum_{j}^Lx \alpha_{j}+x_{ij}
        // adPi(s,f,i)
        // dLogE -= adPi[i] * lngammaAlphaN; // -\sum_i E[z_i] *\sum_j lngamma(x_{ij}+alpha_j)
        dLogE -= adPi(s,0,i) * lngammaAlphaN_1; // -\sum_i \sum_s E[z_isf] *\sum_j lngamma(x_{ij}+alpha_j), f=1
        dLogE -= adPi(s,1,i) * lngammaAlphaN_2; // -\sum_i \sum_s E[z_isf] *\sum_j lngamma(x_{ij}+alpha_j), f=2
      }
    }
    // dLogEAlpha1_[s]
    // dLogEAlpha[s] += gsl_sf_lngamma(dAlpha[j+s]); // \sum_{j=1}^S \log \gamma (\alpha_{j} )  
   dLogEAlpha_final += dLogEAlpha[s] * dWeight[s]; //dLogEAlpha1_[s] // \sum_{j=1}^S \log \gamma (\alpha_{j} )
   dLogEAlpha_final -= gsl_sf_lngamma(dSumAlpha_s[s])* dWeight[s];  //dLogEAlpha_2[s] -\log \gamma ( \sum_{j=1}^S  \alpha_{j} )
 }
  

  for(i = 0; i < N; i++){
    for(s = 0; s < S; s++){
      dLogE += dWeight_f[s,i] * gsl_sf_lngamma(adSumAlphaN[s,i]); // \sum_{n=1}^N E[z_i]* lngamma( \sum_j(x_ij+alpha_j) )
    }
  }
  double reg_term;

  if ((GAMMA_ITA_H==0) && (GAMMA_NU_H==0)) {
    reg_term = 0.0;
  } else {
    double hk = REG_H_LAMBDA(lambda); // computes the value of the regularization term
    reg_term = GAMMA_NU_H*hk - (GAMMA_ITA_H-1)*log(hk); //This was nuh *hg-etah*log (hk), it is now corrected to (GAMMA_ITA_H-1)!!!
  }



  //should it be (GAMMA_ITA-1)*dSumLambda???
  return dLogE + dLogEAlpha_final + // complete data likelihood term
    GAMMA_NU*dSumAlpha - (GAMMA_ITA - 1) * dSumLambda + reg_term; //prior term, ITA corrected to ITA-1 !!!
}


/* A function to return the gradient for the BFGS, derivative of the 
 * expected negative log posterior wrt lambda
 * lambda are the current values, lparams contains the current parameter values (E(z_i))
 * 
 */

// [[Rcpp::export]]
Rcpp::NumericVector neg_log_derive_evidence_lambda_pi_shift_flip(Rcpp::NumericVector ptLambda,
                                                Rcpp::List lparams)
{
  Rcpp::IntegerMatrix aanX = as<Rcpp::IntegerMatrix>(lparams["data"]); // N x S
  arma::Cube<double> adPi = as<arma::cube>(lparams["pi"]); // Sx2xN
  double GAMMA_ITA = as<double>(lparams["eta"]);
  double GAMMA_NU = as<double>(lparams["nu"]);
  double GAMMA_ITA_H = as<double>(lparams["etah"]);
  double GAMMA_NU_H = as<double>(lparams["nuh"]);
  Rcpp::List hkm = as<List>(lparams["hkm"]);
  Rcpp::IntegerVector hkm_index = as<Rcpp::IntegerVector>(lparams["hkm_index"]);
  
  int i, j, s, f, k; // k is j'
  
  const int Lx = aanX.ncol();
  const int N = aanX.nrow();
  const int S = adPi.n_rows;
  const int La = Lx + S - 1;
  
  Rcpp::NumericVector g(La);
  Rcpp::NumericVector adDeriv(La); //derivative for each j
  Rcpp::NumericMatrix adStore(S,N); // \sum_j^S x_ij + \sum_j^S \alpha_j
  Rcpp::NumericVector adAlpha(La);
  
  Rcpp::NumericVector dSumStore(La); // \sum_n^N E[z_i]* psi( \sum_j^S x_ij + \sum_j^S \alpha_j )
  Rcpp::NumericVector dStore(S); // sum of alpha over k
  Rcpp::NumericVector xStore(N); // sum of x_i over k
  Rcpp::NumericVector dWeight_j(La); // dWeight_j
  Rcpp::NumericVector dWeight(S); // dWeight
  Rcpp::NumericMatrix dWeight_f(S,N);
  

  for(s = 0; s < S; s++){
    dStore[s] = 0.0;
    double tmp=0.0;
    for (i = 0; i < N; i++) {
      adStore[s,i] = 0.0;
      double tmp2 = 0.0;
      for(f = 0; f < 2; f++){
        tmp += adPi(s,f,i); //sum over i and f
        tmp2 += adPi(s,f,i);
      }
      dWeight_f[s,i]=tmp2;
    }
    dWeight[s]=tmp;
   
  }
  
  for(j = 0; j < La; j++){
    adAlpha[j] = exp(ptLambda[j]);
    dSumStore[j] = 0.0;
    // Rcpp::Rcout << "j: " << j <<" ";
    for(i = 0; i < N; i++){
      for(f = 0; f < 2; f++){
        for(s =( j - Lx + 1); s < (j+1) ; s++){
          if(s > -1 && s < S ){
            if(i==0 && f==0){
              // Rcpp::Rcout << "s: " << s <<" ";
              }
           
            dWeight_j[j] += adPi(s,f,i);
          }
        }
      }
    }
    // Rcpp::Rcout  <<"\n";
  }
  
  
  
  for(s = 0; s < S; s++){
    // Rcpp::Rcout << "s: " << s <<" ";
    for(k = 0; k < Lx; k++ ){ // Why K, should be j, does not really matter here
      // Rcpp::Rcout  << k+s <<" ";
      dStore[s] += adAlpha[k+s];
      }
    // Rcpp::Rcout <<"\n";
  }

  
  for(i = 0; i < N; i++){
    xStore[i] = 0.0;
    for(k = 0; k < Lx; k++ ){
      xStore[i] += aanX(i, k );
    }
    for(s = 0; s < S; s++){
      adStore[s,i]=xStore[i]+dStore[s];
    }
  }
  
  
  
  for (j = 0; j < La; j++) {
    // Rcpp::Rcout << "j: " << j << " ";
    double alphaS0 = gsl_sf_psi(adAlpha[j]);
    adDeriv[j] = dWeight_j[j] * alphaS0; // adDeriv_1[j]
    // Rcpp::Rcout << "s_start:" << (j-Lx+1) << " s_end:"<< j << " ";
    for(s = (j - Lx + 1); s < (j + 1); s++){
      if(s > -1 && s < S){
        // Rcpp::Rcout << j-s << " ";
        // Rcpp::Rcout << Lx-1 -j+s << " ";
        for (i = 0; i < N; i++) {
          int dN_1 = aanX[i, j - s ]; //c++ indexing starts from 0 -> -1
          int dN_2 = aanX[i, Lx - 1 - j + s]; //c++ indexing starts from 0 -> -1
          
          double dAlphaN_1 = adAlpha[j] + dN_1;
          double dAlphaN_2 = adAlpha[j] + dN_2;
          
          double psiAlphaN_1 = dN_1 ? gsl_sf_psi(dAlphaN_1) : alphaS0;
          double psiAlphaN_2 = dN_2 ? gsl_sf_psi(dAlphaN_2) : alphaS0;
          
          adDeriv[j] -= adPi(s,0,i) * psiAlphaN_1; //adDeriv_2[j], f=1
          adDeriv[j] -= adPi(s,1,i) * psiAlphaN_2; //adDeriv_2[j], f=2
          
        } //i
        
        
      } //if true s
    } //s
    // Rcpp::Rcout <<"\n";
  } //j
  
  for(j = 0; j < La; j++){
    // Rcpp::Rcout << "j: " << j << " ";
    // Rcpp::Rcout << "s_start:" << (j-Lx+1) << " s_end:"<< j << " ";
    for(s = (j - Lx + 1); s < (j + 1); s++){
      
      if(s > -1 && s < S ){
        // Rcpp::Rcout << s << " ";
       for (i = 0; i < N; i++){
          dSumStore[j] += dWeight_f[s,i] * gsl_sf_psi(adStore[s,i]);
        }
      }
    }
    // Rcpp::Rcout <<  "\n";
 }
  
  
  Rcpp::NumericVector dStore_final(La);
  for(j = 0; j < La; j++){
    dStore_final[j] = 0.0;
    for(s = (j - Lx + 1); s < (j + 1); s++){
        if(s > -1 && s < S){
          dStore_final[j] += dWeight[s] * gsl_sf_psi(dStore[s]); //dStore_2
          }
        
    }
  }
 

  double hk = REG_H_LAMBDA(ptLambda);

  if (hkm_index[0] < hkm.size()) hkm[hkm_index[0]] = hk;
  else Rprintf("hkm_index exceeds hkm vector length!!! hkm_index: %i, hkm.length: %i\n",
               hkm_index[0], hkm.size());

  hkm_index[0] += 1;

  for (j = 0; j < La; j++) {

    double reg_term;

    if ((GAMMA_ITA_H==0) && (GAMMA_NU_H==0)) {
      reg_term = 0.0;
    } else {
      double gjk = REG_DERIV_LAMBDA_G(ptLambda, j); //deriv. of h_kj^(m) wrt \alpha_kj^(m)
      reg_term = GAMMA_NU_H*gjk - (GAMMA_ITA_H-1)*gjk/hk;
    }
    //This is now corrected, GAMMA_ITA converted to GAMMA_ITA-1
    double value = adAlpha[j] *
      (GAMMA_NU + adDeriv[j] - dStore_final[j] + dSumStore[j] + reg_term) - ( GAMMA_ITA -1 ); // should be (GAMMA_ITA -1)

    g[j] = value;
  }
  return g;
}


/*The function below computes the value of the lower bound to be optimized
 * i.e. the whole expected log posterior Q with only shifting, no flipping. 
 * Computes the term only depending on \alpha_{k}^m or \lambda_l^m
 * . Does not consider \pi. Computes the value for one cluster k and datatype m
 * lambda are the current values, lparams contains the current parameter values (E(z_i))
 * List of 8
 $ pi       : Z array(dim=c(S,N) )These are z_k for all i and s 
 $ data     : int [1:1000, 1:L_x]
 $ nu       : num 1
 $ etah     : num 1
 $ nuh      : num 10
 $ hkm      :List 
 ..$ : num 
 ..$ : num 
 $ hkm_index: int 
 * 
 */

// [[Rcpp::export]]
double neg_log_evidence_lambda_pi_shift(Rcpp::NumericVector lambda, Rcpp::List lparams)
{
  
  int i, s, j;
  
  Rcpp::IntegerMatrix aanX = as<Rcpp::IntegerMatrix>(lparams["data"]); // N x L_x
  Rcpp::NumericMatrix adPi = as<Rcpp::NumericMatrix>(lparams["pi"]); //S x N
  
  double GAMMA_ITA = as<double>(lparams["eta"]);
  double GAMMA_NU = as<double>(lparams["nu"]);
  double GAMMA_ITA_H = as<double>(lparams["etah"]);
  double GAMMA_NU_H = as<double>(lparams["nuh"]);
  
  const int Lx = aanX.ncol(); // L_x
  const int N = aanX.nrow();
  const int S = adPi.nrow();
  const int La = Lx + S - 1; //Wx=Wa-(S-1)*B
  
  double dLogE = 0.0;
  
  
  double dLogEAlpha_final = 0.0;
  double dSumAlpha = 0.0; // sum of alpha \sum_{j=1}^S \alpha_{j}
  double dSumLambda = 0.0; // sum of lambda?

 
  Rcpp::NumericVector dWeight(S);
  Rcpp::NumericVector dAlpha(lambda.length());

  Rcpp::NumericMatrix adSumAlphaN(S,N);
  Rcpp::NumericVector dSumAlpha_s(S);
  Rcpp::NumericVector dLogEAlpha(S); 
  
  //Compute dWeight_[s]
  for(s = 0; s < S; s++){
    dSumAlpha_s[s]=0.0;
    double tmp=0.0;
    for (i = 0; i < N; i++) {
      adSumAlphaN[s,i]=0.0;
      tmp += adPi[s,i]; //sum over i 
    }
    dWeight[s]=tmp;
    //Rcpp::Rcout << dWeight(s) << "";
  }
  
  
  
  for (j = 0; j < La; j++) {
    dAlpha[j] = exp(lambda[j]);
    dSumLambda += lambda[j];
    dSumAlpha += dAlpha[j];
  }
  
  
  
  for(s = 0; s < S; s++){
    for (j = 0; j < Lx; j++) {
      
      dSumAlpha_s[s] += dAlpha[j+s];
      
      //dLogEAlpha1_[s]
      dLogEAlpha[s] += gsl_sf_lngamma(dAlpha[j+s]); // \sum_{j=1}^S \log \gamma (\alpha_{j} )
      
      const double lngammaAlpha0_1 = gsl_sf_lngamma(dAlpha[j+s]); //\log \gamma (\alpha_{js} ) 
      
      for (i = 0; i < N; i++) {
        
        const double dN = aanX[i, j]; // x_{ij}
        const double dAlphaN_1 = dAlpha[j+s] + dN; // \alpha_{j}+x_{ij} f=1
        const double lngammaAlphaN_1 = dN ? gsl_sf_lngamma(dAlphaN_1) : lngammaAlpha0_1; //if dN exists or is non-zero, compute lnGamma(dAlphaN), else compute lnGamma(DaLpha)
       
        adSumAlphaN[s,i] += dAlphaN_1; // \sum_{j}^Lx \alpha_{j}+x_{ij}
        dLogE -= adPi[s,i] * lngammaAlphaN_1; // -\sum_i \sum_s E[z_isf] *\sum_j lngamma(x_{ij}+alpha_j), 
        
      } //i
    } //j
    //dLogEAlpha1_[s]
    
    dLogEAlpha_final += dLogEAlpha[s] * dWeight[s]; //dLogEAlpha1_[s] // \sum_{j=1}^S \log \gamma (\alpha_{j} )
    dLogEAlpha_final -= gsl_sf_lngamma(dSumAlpha_s[s]) * dWeight[s];  //dLogEAlpha_2[s] -\log \gamma ( \sum_{j=1}^S  \alpha_{j} )
  } //s
  
  
  for(i = 0; i < N; i++){
    for(s = 0; s < S; s++){
      dLogE += adPi[s,i] * gsl_sf_lngamma(adSumAlphaN[s,i]); // \sum_{n=1}^N E[z_i]* lngamma( \sum_j(x_ij+alpha_j) )
    }
  }
  
  double reg_term;
  
  if ((GAMMA_ITA_H==0) && (GAMMA_NU_H==0)) {
    reg_term = 0.0;
  } else {
    double hk = REG_H_LAMBDA(lambda); // computes the value of the regularization term
    reg_term = GAMMA_NU_H*hk - (GAMMA_ITA_H-1)*log(hk); //This was nuh *hg-etah*log (hk), it is now corrected to (GAMMA_ITA_H-1)!!!
  }
  
  //should it be (GAMMA_ITA-1)*dSumLambda???
  return dLogE + dLogEAlpha_final + // complete data likelihood term
    GAMMA_NU*dSumAlpha - (GAMMA_ITA - 1) * dSumLambda + reg_term; //prior term
}


/* With shift only. A function to return the gradient for the BFGS, derivative of the 
 * expected negative log posterior wrt lambda
 * lambda are the current values, lparams contains the current parameter values (E(z_i))
 * 
 */

// [[Rcpp::export]]
Rcpp::NumericVector neg_log_derive_evidence_lambda_pi_shift(Rcpp::NumericVector ptLambda,
                                                                 Rcpp::List lparams)
{
  Rcpp::IntegerMatrix aanX = as<Rcpp::IntegerMatrix>(lparams["data"]); // N x Lx
  Rcpp::NumericMatrix adPi = as<Rcpp::NumericMatrix>(lparams["pi"]); // S x N
  double GAMMA_ITA = as<double>(lparams["eta"]);
  double GAMMA_NU = as<double>(lparams["nu"]);
  double GAMMA_ITA_H = as<double>(lparams["etah"]);
  double GAMMA_NU_H = as<double>(lparams["nuh"]);
  Rcpp::List hkm = as<List>(lparams["hkm"]);
  Rcpp::IntegerVector hkm_index = as<Rcpp::IntegerVector>(lparams["hkm_index"]);
  
  int i, j, s, k; // k is j'
  
  const int Lx = aanX.ncol();
  const int N = aanX.nrow();
  const int S = adPi.nrow();
  const int La = Lx + S - 1;
  
  Rcpp::NumericVector g(La); //gradient
  Rcpp::NumericVector adDeriv(La); //derivative for each j
  Rcpp::NumericMatrix adStore(S,N); // \sum_j^S x_ij + \sum_j^S \alpha_j
  Rcpp::NumericVector adAlpha(La);
  
  Rcpp::NumericVector dSumStore(La); // \sum_n^N E[z_i]* psi( \sum_j^S x_ij + \sum_j^S \alpha_j )
  Rcpp::NumericVector dStore(S); // sum of alpha over k
  Rcpp::NumericVector xStore(N); // sum of x_i over k
  Rcpp::NumericVector dWeight_j(La); // dWeight_j
  Rcpp::NumericVector dWeight(S); // dWeight
  Rcpp::NumericVector dStore_final(La);
  
  
  //Compute dWeight_[s]
  for(s = 0; s < S; s++){
    dStore[s] = 0.0;
    double tmp = 0.0;
    for (i = 0; i < N; i++) {
      adStore[s,i] = 0.0;
      
      tmp += adPi[s,i]; //sum over i a
    }
    dWeight[s]=tmp;
  }
  
  for(j = 0; j < La; j++){
    adAlpha[j] = exp(ptLambda[j]);
    dSumStore[j] = 0.0;
    dWeight_j[j] = 0.0;
    //Rcpp::Rcout << "j: " << j <<" ";
    for(i = 0; i < N; i++){
        for(s =( j - Lx + 1); s < (j+1) ; s++){
          if(s > -1 && s < S ){
            //if(i==0 && f==0){
              // Rcpp::Rcout << "s: " << s <<" ";
            //}
            dWeight_j[j] += adPi[s,i];
          }
        } //s
      } //i
    // Rcpp::Rcout  <<"\n";
  } //j
  
  for(s = 0; s < S; s++){
    // Rcpp::Rcout << "s: " << s <<" ";
    for(k = 0; k < Lx; k++ ){
      // Rcpp::Rcout  << k+s <<" ";
      dStore[s] += adAlpha[k+s];
    }
    // Rcpp::Rcout <<"\n";
  }
  
  
  for(i = 0; i < N; i++){
    xStore[i] = 0.0;
    for(k = 0; k < Lx; k++ ){
      xStore[i] += aanX[i, k ];
    }
    for(s = 0; s < S; s++){
      adStore[s,i]=xStore[i]+dStore[s];
    }
  }
  
  
  
  for (j = 0; j < La; j++) {
    // Rcpp::Rcout << "j: " << j << " ";
    double alphaS0 = gsl_sf_psi(adAlpha[j]);
    adDeriv[j] = dWeight_j[j] * alphaS0; // adDeriv_1[j]
    // Rcpp::Rcout << "s_start:" << (j-Lx+1) << " s_end:"<< j << " ";
    for(s = (j - Lx + 1); s < (j + 1); s++){
      if(s > -1 && s < S){
        // Rcpp::Rcout << j-s << " ";
        // Rcpp::Rcout << Lx-1 -j+s << " ";
        for (i = 0; i < N; i++) {
          int dN_1 = aanX[i, j - s ]; //c++ indexing starts from 0 -> -1
          
          double dAlphaN_1 = adAlpha[j] + dN_1;
      
          
          double psiAlphaN_1 = dN_1 ? gsl_sf_psi(dAlphaN_1) : alphaS0;
          
          adDeriv[j] -= adPi[s,i] * psiAlphaN_1; //adDeriv_2[j], f=1
          
          
        } //i
      } //if true s
    } //s
    // Rcpp::Rcout <<"\n";
  } //j
  
  for(j = 0; j < La; j++){
    // Rcpp::Rcout << "j: " << j << " ";
    // Rcpp::Rcout << "s_start:" << (j-Lx+1) << " s_end:"<< j << " ";
    for(s = (j - Lx + 1); s < (j + 1); s++){
      
      if(s > -1 && s < S ){
        // Rcpp::Rcout << s << " ";
        for (i = 0; i < N; i++){
          dSumStore[j] += adPi[s,i] * gsl_sf_psi(adStore[s,i]);
        }
      }
    }
    // Rcpp::Rcout <<  "\n";
  }
  
  
  
  for(j = 0; j < La; j++){
    dStore_final[j] = 0.0;
    for(s = (j - Lx + 1); s < (j + 1); s++){
      if(s > -1 && s < S){
        dStore_final[j] += dWeight[s] * gsl_sf_psi(dStore[s]); //dStore_2
      }
    }
  }
  
  
  double hk = REG_H_LAMBDA(ptLambda);
  
  if (hkm_index[0] < hkm.size()) hkm[hkm_index[0]] = hk;
  else Rprintf("hkm_index exceeds hkm vector length!!! hkm_index: %i, hkm.length: %i\n",
               hkm_index[0], hkm.size());
  
  hkm_index[0] += 1;
  
  for (j = 0; j < La; j++) {
    
    double reg_term;
    
    if ((GAMMA_ITA_H==0) && (GAMMA_NU_H==0)) {
      reg_term = 0.0;
    } else {
      double gjk = REG_DERIV_LAMBDA_G(ptLambda, j); //deriv. of h_kj^(m) wrt \alpha_kj^(m)
      reg_term = GAMMA_NU_H*gjk - (GAMMA_ITA_H-1)*gjk/hk;
    }
   
    double value = adAlpha[j] *
      (GAMMA_NU + adDeriv[j] - dStore_final[j] + dSumStore[j] + reg_term) - ( GAMMA_ITA -1 ); // should be (GAMMA_ITA -1)
    
    g[j] = value;
  }
  return g;
}



// /* DONE Computes the E-step, i.e. the posterior probabilities of the cluster labels
//  * Z=Ez K x S x N matrix (from the previous E-step, these are not used?)
//  * data=binned.data list of M matrixes of size N x Lx
//  * W=weights vector of length K, \pi_k^{old}
//  * lambda=lambda list of M, each K times L matrix
//  *
//  * returns K x S x N matrix Z
//  *
//  * Why offset is subtracted, for numerical reasons, but the offset not added later?
//  */

// [[Rcpp::export]]
arma::Cube<double> calc_z_shift(arma::Cube<double> Z, Rcpp::List data,
                     Rcpp::NumericVector W, Rcpp::NumericVector xi, Rcpp::List Lambda)
{
  // Z is K x S x N
  int i, j, k, m, s, l;
  Rcpp::IntegerMatrix temp = as<Rcpp::IntegerMatrix>(data[0]);
  const int N = temp.nrow();
  const int K = W.length(); // number of clusters
  const int M = data.size();
  const int S = xi.length();
  
  arma::Cube<double> evidence_matrix(M, K, S); 
  // Rcpp::NumericMatrix evidence_matrix(M, K); // will contain the log p(\mathbf{x}_{i}^m| \theta) for sample i

  Rcpp::NumericVector Lx(M);
  Rcpp::NumericVector La(M);
  for (m = 0; m < M ; m++) {
    temp = as<Rcpp::IntegerMatrix>(data[m]); //N x L matrix
    Lx[m] = temp.ncol(); //L
    La[m] = Lx[m] + S - 1; //assume all Lx are the same
  }

  Rcpp::List LngammaLambda0(M); // lngamma ( \alpha_{jk}^{(m)}  )

  // Compute lngammaalpha_jkms
  for (m = 0; m < M; m++) {
    // Rcpp::Rcout << "m: "<<m << "\n";
    Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x La matrix
    arma::Cube<double> LngammaLambda0_matrix(K, Lx[m],S); // lngamma ( \alpha_{jk}  )
    for(k = 0; k < K; k++){
      // Rcpp::Rcout << "k: "<<k << "\n";
      for(s = 0; s < S; s++){
        // Rcpp::Rcout << "s: "<<s << "\n";
        for(j = 0; j < Lx[m]; j++){
          // Rcpp::Rcout <<"j: "<<j << ", j+s: "<<j+s << "\n";
          const double dAlpha = exp(Lambda_matrix(k, j+s));
          LngammaLambda0_matrix(k, j, s) = gsl_sf_lngamma(dAlpha);
        } //j
      } //s
    } //k
    LngammaLambda0(m) = LngammaLambda0_matrix;
  } //m

  for (i = 0; i < N; i ++) {
   
    double dSum = 0.0;
    Rcpp::NumericVector offset(M); //save the smallest negLogEvidence for each cluster(largest absolute value)
    // Compute the evidence matrix, the DirichletMultinomial pdf for given i, m and k
    for (m = 0; m < M; m++) {

      offset[m] = BIG_DBL; //1.0e9

      Rcpp::IntegerMatrix data_matrix = as<Rcpp::IntegerMatrix>(data[m]); // N x Lx matrix
      Rcpp::NumericMatrix lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x La matrix
      // Rcpp::NumericMatrix Ln = as<Rcpp::NumericMatrix>(LngammaLambda0[m]); // K x L lngamma ( \alpha_{jk}^{(m)}  )
      arma::Cube<double> LngammaLambda0_matrix = LngammaLambda0(m); //K x Lx x S
      Rcpp::IntegerVector data_row = data_matrix(i, _); // one row of X of length Lx
      for (k = 0; k < K; k++) {
      
        //Computes the logarithm of -p(\mathbf{x}_i|z_{ik}=1,\bm{\theta}) but of only those
        // terms that depend on \alpha
        // I.e. computes the negative logarithm of the unnormalized DirichletMultinomial pdf value
        //logarithm due to numerical purposes
        // computes the - log DirMulti(x_i^m | alpha_k^m), terms not depending on alpha excluded
        for(s = 0; s < S; s++){
        
          Rcpp::NumericVector lambda_ks(Lx[m]);
          for(l = 0; l < Lx[m]; l++){
           
            lambda_ks[l]=lambda_matrix[k, s+l]; //lambda_{ks}^m
            }
          // LngammaLambda0_matrix K x Lx x S
          double dNegLogEviI = neg_log_evidence_i(data_row, lambda_ks,
                                               as<Rcpp::NumericVector>( Rcpp::wrap(LngammaLambda0_matrix.slice(s).row(k)) ) );
         
          if (dNegLogEviI < offset[m]) //1.0e9 Can this explode?
            offset[m] = dNegLogEviI; // offset[m] is the smallest dNegLogEviI over k
          evidence_matrix(m, k, s) = dNegLogEviI; //M x K x S
        
        } //s
        //data_row x_i^m 1 x L
        // lambda_matrix (k, _) log alpha_k^m 1 x L
        // LngammaLambda0_matrix(k, _) 1 x L lngamma ( \alpha_{jk}^{(m)}  )


        
      } //over K
    } //over M

    for (k = 0; k < K; k++) {
      //Rcpp::Rcout << "k: "<<k << "\n";
      for(s = 0; s < S; s++){
        //Rcpp::Rcout << "s: "<<s << "\n";
        Z(k, s, i) = 0; //K x S x N
        for (m = 0; m < M; m++) {
          //Rcpp::Rcout << "m: "<<m << "\n";
          Z(k, s, i) += (-(evidence_matrix(m, k, s) - offset[m])); //why offset is substracted?? For numerical reasons?
        }
        
        Z(k, s, i) = W[k]*xi[s] * exp(Z(k, s, i)); //back from logaritm, multiply by the shift probs i.e. \xi_k
        dSum += Z(k, s, i); //normalization constant
      } //s
    } // k
    
    for (k = 0; k < K; k++){
      for(s = 0; s < S; s++){
        Z(k, s, i) /= dSum;
        }
    }
  } // i
  return Z; //K x N matrix
}

/* DONE Function whose convergence is checked in each EM iteration
 * Computes log p(\theta,Z| X) \propto log p(X,Z|theta) + log p(theta)
 * the terms in log p(theta) that do not depend on alpha do not affect the convergence, 
 * why they are computed? In addition, one constant term is missing?
 * (\eta-1) term needs to be corrected !!??
 * W weights 1xK, these are \pi_k values
 * lambda list of length M, these are K x La matrices
 * xi are the shift state probabilities
 * binned.data list of length M, NxLx matrices
 */

// [[Rcpp::export]]
double neg_log_likelihood_shift(Rcpp::NumericVector W, Rcpp::NumericVector xi, Rcpp::List Lambda,
                          Rcpp::List data, double eta, double nu,
                          double etah, double nuh, int S)
{
   
  Rcpp::IntegerMatrix temp = as<Rcpp::IntegerMatrix>(data[0]); // N x Lx
  const int N = temp.nrow();
    // Rcpp::Rcout << "N: "<<N << "\n";
  const int K = W.length();
    // Rcpp::Rcout << "K: "<<K << "\n";
  const int M = data.size();
    // Rcpp::Rcout << "M: "<<M << "\n";
  int i, j, k, m, s;
  
  Rcpp::NumericVector Lx(M);
  Rcpp::NumericVector La(M);
  for (m = 0; m < M; m++) {
    temp = as<Rcpp::IntegerMatrix>(data[m]);
    Lx[m] = temp.ncol();
    La[m] = Lx[m] + S - 1; //assume all Lx are the same
  }
  
  

  Rcpp::List LngammaLambda0(M);
  
  // return -dRet - dL5 - dL6 - dL7 - dL8 - regterm1 - regterm2;
  /* -log \sum_k^K(\pi_k*p(x|\theta_k))+M*L*K*lng(eta)- eta*K*M*L*log(nu)
   * + nu* \sum_j^L \alpha_kj^m -eta * \sum_j^L \lambda_kj^m
   -M * K * (etah * log(nuh) - gsl_sf_lngamma(etah)) - \sum_k^K ((etah - 1) * log(hkm) - nuh*hkm) */
  double dRet = 0.0; //dRet of log \sum_k^K(\pi_k*p(x|\theta_k)) without Gamma (J_i^m-1) term
  double dL5 = 0.0; // -M*L*K*lng(eta) from prior
  double dL6 = 0.0; // eta*K*M*L*log(nu) from prior
  double dL7 = 0.0; //-nu* \sum_j^L \alpha_kj^m from prior
  double dL8 = 0.0; // eta * \sum_j^L \lambda_kj^m from prior
  double regterm1 = 0.0; //M * K * (etah * log(nuh) - gsl_sf_lngamma(etah)); from prior
  double regterm2 = 0.0; // \sum_k^K ((etah - 1) * log(hkm) - nuh*hkm);  form prior
  
  //const int S = adPi.n_rows;
  arma::Cube<double> LogBAlpha(M,K,S); 
  //NumericMatrix LogBAlpha(M, K); // // \sum_j^L lng (alpha_jk) -lngamma( \sum_j^L \alpha_jk)
  
  Rcpp::NumericVector dSumAlphaK(S);
  
  for (m = 0; m < M; m++) {
    arma::Cube<double> LngammaLambda0_matrix(K, Lx[m],S); // lng (alpha_jk)
    Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x La matrices
    
    for (k = 0; k < K; k++){
      for(s = 0; s < S; s++){
        dSumAlphaK[s] = 0.0; // separate sum for each k,m,s
        LogBAlpha(m, k, s) = 0.0;
        for (j = 0; j < Lx[m]; j++){
          double dAlpha = exp(Lambda_matrix(k, j + s )); // \ alpha_jks
          double lngammaAlpha = gsl_sf_lngamma(dAlpha); // lng (alpha_jk)
          LngammaLambda0_matrix(k, j, s) = lngammaAlpha; // // lng (alpha_jk)
          
          dSumAlphaK[s] += dAlpha; // \sum_j^L \alpha_jk
          LogBAlpha(m, k, s) += lngammaAlpha; // LogBAlpha_1
        } //over j
        LogBAlpha(m, k,s) -= gsl_sf_lngamma(dSumAlphaK[s]); // -LogBAlpha_2
      } //over s
    } // over k
    LngammaLambda0(m) = LngammaLambda0_matrix; // list of length M
  } // over m
  
    // Rcpp::Rcout << "LogBAlpha: "<< LogBAlpha(0,0,10) << "\n";
  
  for (i = 0; i < N; i++) {
    double dProb = 0.0;
   
    arma::Cube<double> LogStore(M, K, S);
    Rcpp::NumericVector offset(M);
    // Computes the log p(x_i^m|, z_ik=1, theta) for k and m, results in matrix LogStore
    for (m = 0; m < M; m++) {
      offset[m] = -BIG_DBL;
      double dSum = 0.0; // \sum_j^L x_ij
      double dFactor = 0.0; // \sum_j^L lng( x_ij +1) - lngamma( \sum_j^L x_ij +1)
    
      Rcpp::IntegerMatrix data_matrix = as<Rcpp::IntegerMatrix>(data[m]);
      Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]);
      arma::Cube<double> LngammaLambda0_matrix = LngammaLambda0(m); //size K x Lx x S
      // NumericMatrix LngammaLambda0_matrix = as<NumericMatrix>(LngammaLambda0[m]); // lng (alpha_jk) K x L matrix
      
      for (j = 0; j < Lx[m]; j++) {
        dSum += data_matrix(i, j); // dSum_m, this is J_i^m
        // dFactor does not depend on parameters, why computed?
        // dFactor_1
        dFactor += gsl_sf_lngamma(data_matrix(i, j) + 1.0); // \sum_j^L lng( x_ij +1),
      }
      // dFactor does not depend on parameters, why computed?
      dFactor -= gsl_sf_lngamma(dSum + 1.0);  // \sum_j^L lng( x_ij +1) - lngamma( \sum_j^L x_ij +1)
      
      for (k = 0; k < K; k++) {
        Rcpp::NumericVector dSumAlphaKN(S); // \sum_j^L \alpha_jk + x_ij
        Rcpp::NumericVector dLogBAlphaN(S); // lng( \alpha_jk + x_ij ) - lng( \sum_j^L \alpha_jk + x_ij )
        for(s = 0; s < S; s++){
          dSumAlphaKN[s]=0.0;
          dLogBAlphaN[s]=0.0;
          for (j = 0; j < Lx[m]; j++) {
            int countN = data_matrix(i, j); // x_ij
            double dAlphaN = exp(Lambda_matrix(k, j+s)) + countN; // \alpha_jk + x_ij
            dSumAlphaKN[s] += dAlphaN; // \sum_j^L \alpha_jk + x_ij
            dLogBAlphaN[s] += countN ? gsl_sf_lngamma(dAlphaN) : LngammaLambda0_matrix(k, j, s); // dLogBAlphaN_1[m,k,s]
          }
          dLogBAlphaN[s] -= gsl_sf_lngamma(dSumAlphaKN[s]); // // dLogBAlphaN_2[m,k,s]
          
          // LogStore (m,k) seems to be p(\mathbf{x}_i|\theta) 
          // LogStore= lng( \alpha_jk + x_ij ) - lng( \sum_j^L \alpha_jk + x_ij )
          // -\sum_j^L lng (alpha_jk) + lngamma( \sum_j^L \alpha_jk) - \sum_j^L lng( x_ij +1)
          LogStore(m, k, s) = dLogBAlphaN[s] - LogBAlpha(m, k,s) - dFactor; // the positive? log marginal likelihood of sample x_i^m for cluster k
          if (LogStore(m, k,s) > offset(m))
            offset(m) = LogStore(m, k,s); //offset will be the largest of LogStore(m,:,:)
        } //over S
        
        // Rcpp::Rcout << "LogStore: "<< LogStore(0,0,10) << "\n";
      } //over K
    } //over M
    
    Rcpp::NumericVector dProb_S(K);
    for(k = 0; k < K; k++){
      dProb_S[k]=0.0;
      for (s = 0; s < S; s++) {
        
        // dProb_S[k] += xi[s]*exp( sum( Rcpp::wrap( LogStore(arma::span::all,arma::span(k),arma::span(s) ))) - sum(offset) ); //sum of column k over m rows
        dProb_S[k] += xi[s]*exp( sum( as<Rcpp::NumericVector>(Rcpp::wrap( LogStore.slice(s).col(k) ) )) - sum(offset) ); //sum of column k over m rows
        // Rcpp::Rcout << "dProb_S: "<< dProb_S[k] << "\n";
        
      }
      // normalize mixing weights
      // Rcpp::Rcout << "W[k]: "<< W[k] << "\n";
      
      double piK = W[k]/N;
      
      // Rcpp::Rcout << "piK: "<< piK << "\n";
      dProb += piK*dProb_S[k]; //sum of column k over m rows
      
      // Rcpp::Rcout << "dProb: "<< dProb << "\n";
    } // k
    
    
    dRet += log(dProb)+Rcpp::sum(offset); //logarithm of the normalization term, sum over all n
    // Rcpp::Rcout << "dRet: "<< dRet << "\n";
    
  } //over N
  
  
  
  
  dL5 = -sum(La) * K * gsl_sf_lngamma(eta); // -M*L*K*lng(eta) sum(La)=M*La
  dL6 = eta * K * sum(La) * log(nu);  //eta*K*M*L*log(nu)
  
  if ((etah!=0) || (nuh!=0)) {
    regterm1 = M * K * (etah * log(nuh) - gsl_sf_lngamma(etah));
  }
  
  for (m = 0; m < M; m++) {
    Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x La matrix
    
    for (k = 0; k < K; k++) {
      
      if ((etah!=0) || (nuh!=0)) {
        // exp(Lambda_matrix(i, _) are alpha_k, all L elements
        // diff function in c++ computes the difference of elements of a vector
        double hkm = sum( diff( exp(Lambda_matrix(k, _)) ) * diff(exp(Lambda_matrix(k, _))));
        regterm2 += (etah - 1) * log(hkm) - nuh*hkm; // \sum_k^K (etah - 1) * log(hkm) - nuh*hkm;
      }
      
      dL7 += sum(exp(Lambda_matrix(k, _))); // \sum_j^L \alpha_kj^m
      dL8 += sum(Lambda_matrix(k, _));  //  \sum_j^L \lambda_kj^m
    } // over K
  }  // over M
  dL7 *= -nu; //-nu* \sum_j^L \alpha_kj^m
  dL8 *= (eta-1); // eta * \sum_j^L \lambda_kj^m SHOULD BE (eta-1)
  return -dRet - dL5 - dL6 - dL7 - dL8 - regterm1 - regterm2;
  //
  
}
