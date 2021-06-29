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
  
  double GAMMA_ITA = as<double>(lparams["eta"]);
  double GAMMA_NU = as<double>(lparams["nu"]);
  double GAMMA_ITA_H = as<double>(lparams["etah"]);
  double GAMMA_NU_H = as<double>(lparams["nuh"]);
  
  const int Lx = aanX.ncol(); // L_x
  const int N = aanX.nrow();
  const int S = adPi.n_rows;
  const int La = Lx + S - 1;

  /*dLogE collects the terms \sum_{n=1}^N \sum_{j=1}^S E[z_i] \log \gamma (x_ij+ alpha_j)
   and \sum_{n=1}^N E[z_i]* lng( \sum_j(x_ij+alpha-j) )*/
  double dLogE = 0.0;
  double dLogEAlpha_final = 0.0;
  double dSumAlpha = 0.0; // sum of alpha \sum_{j=1}^S \alpha_{j}
  double dSumLambda = 0.0; // sum of lambda?

  // double dWeight = 0.0; // \sum_{n=1}^N E(z_i)
  Rcpp::NumericVector dWeight(S); //
  Rcpp::NumericVector dAlpha(La); //
  Rcpp::NumericMatrix dWeight_f(S,N); //
  // Rcpp::NumericVector adSumAlphaN(N); // \sum_{j}^S \alpha_{j}+x_{ij}
  Rcpp::NumericMatrix adSumAlphaN(S,N); //
  Rcpp::NumericVector dSumAlpha_s(S); //
  Rcpp::NumericVector dLogEAlpha(S); // \log \gamma ( \sum_{j=1}^S  \alpha_{j} )

  for(s = 0; s < S; s++){
    dSumAlpha_s[s]=0.0;
    dLogEAlpha[s]=0.0;
    dWeight[s]=0.0;
    for (i = 0; i < N; i++) {
      adSumAlphaN(s,i)=0.0;
      dWeight_f(s,i)=0.0;
      for(f = 0; f < 2; f++){
          dWeight[s] += adPi(s,f,i); //sum over i and f
          dWeight_f(s,i) += adPi(s,f,i);
      }
    }
  }

  for (j = 0; j < La; j++) {
    dAlpha[j] = exp(lambda[j]);
    dSumLambda += lambda[j];
    dSumAlpha += dAlpha[j];
  }

  for( s=0; s < S; s++ ){

    for ( j = 0; j < Lx; j++ ) {

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

        adSumAlphaN(s,i) += dAlphaN_1; // \sum_{j}^Lx \alpha_{j}+x_{ij}
        // adPi(s,f,i)
        // dLogE -= adPi[i] * lngammaAlphaN; // -\sum_i E[z_i] *\sum_j lngamma(x_{ij}+alpha_j)
        dLogE -= adPi(s,0,i) * lngammaAlphaN_1; // -\sum_i \sum_s E[z_isf] *\sum_j lngamma(x_{ij}+alpha_j), f=1
        dLogE -= adPi(s,1,i) * lngammaAlphaN_2; // -\sum_i \sum_s E[z_isf] *\sum_j lngamma(x_{ij}+alpha_j), f=2
      } //i
    }  // j
    // dLogEAlpha1_[s]
    // dLogEAlpha[s] += gsl_sf_lngamma(dAlpha[j+s]); // \sum_{j=1}^S \log \gamma (\alpha_{j} )
   dLogEAlpha_final += dLogEAlpha[s] * dWeight[s]; //dLogEAlpha1_[s] // \sum_{j=1}^S \log \gamma (\alpha_{j} )
   dLogEAlpha_final -= gsl_sf_lngamma(dSumAlpha_s[s])* dWeight[s];  //dLogEAlpha_2[s] -\log \gamma ( \sum_{j=1}^S  \alpha_{j} )
 } //s


  for(i = 0; i < N; i++){
    for(s = 0; s < S; s++){
      dLogE += dWeight_f(s,i) * gsl_sf_lngamma(adSumAlphaN(s,i)); // \sum_{n=1}^N E[z_i]* lngamma( \sum_j(x_ij+alpha_j) )
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
Rcpp::NumericVector neg_log_derive_evidence_lambda_pi_shift_flip(Rcpp::NumericVector ptLambda, Rcpp::List lparams)
{
  Rcpp::IntegerMatrix aanX = as<Rcpp::IntegerMatrix>(lparams["data"]); // N x S
  arma::Cube<double> adPi = as<arma::cube>(lparams["pi"]); // Sx2xN
  double GAMMA_ITA = as<double>(lparams["eta"]);
  double GAMMA_NU = as<double>(lparams["nu"]);
  double GAMMA_ITA_H = as<double>(lparams["etah"]);
  double GAMMA_NU_H = as<double>(lparams["nuh"]);
  Rcpp::List hkm = as<List>(lparams["hkm"]);
  Rcpp::IntegerVector hkm_index = as<Rcpp::IntegerVector>(lparams["hkm_index"]);
  Rcpp::List gradient = as<Rcpp::List>(lparams["gradient"]);
  Rcpp::IntegerVector gradient_index = as<Rcpp::IntegerVector>(lparams["gradient_index"]);

  int i, j, s, f, k; // k is j'

  const int Lx = aanX.ncol();
  const int N = aanX.nrow();
  const int S = adPi.n_rows;
  const int La = Lx + S - 1;

  Rcpp::NumericVector g(La); //INIT
  Rcpp::NumericVector adDeriv(La); //INIT derivative for each j
  Rcpp::NumericMatrix adStore(S,N); //INIT \sum_j^S x_ij + \sum_j^S \alpha_j
  Rcpp::NumericVector adAlpha(La); //INIT

  Rcpp::NumericVector dSumStore(La); //INIT \sum_n^N E[z_i]* psi( \sum_j^S x_ij + \sum_j^S \alpha_j )
  Rcpp::NumericVector dStore(S); //INIT sum of alpha over k
  Rcpp::NumericVector xStore(N); //INIT sum of x_i over k
  Rcpp::NumericVector dWeight_j(La); //INIT dWeight_j
  Rcpp::NumericVector dWeight(S); //INIT dWeight
  Rcpp::NumericMatrix dWeight_f(S,N); //INIT


  for(s = 0; s < S; s++){
    dStore[s] = 0.0;
    dWeight[s]=0.0;
    for (i = 0; i < N; i++) {
      adStore(s,i) = 0.0;
      dWeight_f(s,i) = 0.0;
      for(f = 0; f < 2; f++){
        dWeight[s] += adPi(s,f,i); //sum over i and f
        dWeight_f(s,i) += adPi(s,f,i);
      }
    }
  }

  for(j = 0; j < La; j++){
    adAlpha[j] = exp(ptLambda[j]);
    dSumStore[j] = 0.0;
    dWeight_j[j] = 0.0;
    for(i = 0; i < N; i++){
      for(f = 0; f < 2; f++){
        for(s = ( j - Lx + 1); s < ( j + 1 ); s++){
          if(s > -1 && s < S ){
            dWeight_j[j] += adPi(s,f,i);
          }
        }
      }
    }
  }



  for(s = 0; s < S; s++){
    for(k = 0; k < Lx; k++ ){ 
      dStore[s] += adAlpha[k+s];
    }
  }


  for(i = 0; i < N; i++){
    xStore[i] = 0.0;
    for(k = 0; k < Lx; k++ ){
      xStore[i] += aanX(i, k );
    }
    for(s = 0; s < S; s++){
      adStore(s,i)=xStore[i]+dStore[s];
    }
  }



  for (j = 0; j < La; j++) {
    double alphaS0 = gsl_sf_psi(adAlpha[j]);
    adDeriv[j] = dWeight_j[j] * alphaS0; // adDeriv_1[j]
    for(s = (j - Lx + 1); s < (j + 1); s++){
      if(s > -1 && s < S){
        for (i = 0; i < N; i++) {
          int dN_1 = aanX(i, j - s ); //c++ indexing starts from 0 -> -1
          int dN_2 = aanX(i, Lx - 1 - j + s); //c++ indexing starts from 0 -> -1

          double dAlphaN_1 = adAlpha[j] + dN_1;
          double dAlphaN_2 = adAlpha[j] + dN_2;

          double psiAlphaN_1 = dN_1 ? gsl_sf_psi(dAlphaN_1) : alphaS0;
          double psiAlphaN_2 = dN_2 ? gsl_sf_psi(dAlphaN_2) : alphaS0;

          adDeriv[j] -= adPi(s,0,i) * psiAlphaN_1; //adDeriv_2[j], f=1
          adDeriv[j] -= adPi(s,1,i) * psiAlphaN_2; //adDeriv_2[j], f=2

        } //i
      } //if true s
    } //s
  } //j

  for(j = 0; j < La; j++){
    for(s = (j - Lx + 1); s < (j + 1); s++){
      if(s > -1 && s < S ){
        for (i = 0; i < N; i++){
          dSumStore[j] += dWeight_f(s,i) * gsl_sf_psi(adStore(s,i));
        }
      }
    }
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
  
  double gnorm=0.0;
  
  for (j = 0; j < La; j++) {
    gnorm += pow(g[j], 2.0);
  }
  gnorm=sqrt(gnorm);
  // Rprintf("Gradient norm: %f\n", gnorm);
  
  if (gradient_index[0] < gradient.size()) {
    gradient[gradient_index[0]] = gnorm;
  }
  else{
    Rprintf("gradient_index exceeds gradient vector length!!! gradient_index: %i, gradient.length: %i\n",
            gradient_index[0], gradient.size());
  }
  gradient_index[0] += 1;
  
  
  return g;
}

//Z is now length 2 list of arrays/arma::cubes of size K x S x N

// [[Rcpp::export]]
Rcpp::List calc_z_shift_flip(Rcpp::List Z, Rcpp::List data,
                                Rcpp::NumericVector W, Rcpp::NumericVector xi, Rcpp::NumericVector zeta, Rcpp::List Lambda)
{
  // Z is K x S x N
  int i, j, k, m, s, f, l;
  Rcpp::IntegerMatrix temp = as<Rcpp::IntegerMatrix>(data[0]);
  const int N = temp.nrow();
  const int K = W.length(); // number of clusters
  const int M = data.size();
  const int S = xi.length();
  const int F = zeta.length();
  
  arma::Cube<double> Z_1=as<arma::cube>(Z[0]);
  arma::Cube<double> Z_2=as<arma::cube>(Z[1]);
  
  Rcpp::List evidence_matrix(F);  
  arma::Cube<double> evidence_matrix_1(M, K, S);
  arma::Cube<double> evidence_matrix_2(M, K, S);
  
  Rcpp::NumericVector Lx(M);
  Rcpp::NumericVector La(M);
  for (m = 0; m < M ; m++) {
    temp = as<Rcpp::IntegerMatrix>(data[m]); //N x L matrix
    Lx[m] = temp.ncol(); //L
    La[m] = Lx[m] + S - 1; //assume all Lx are the same
  }
  
  Rcpp::List LngammaLambda0_1(M); // for f=1 no flip
  Rcpp::List LngammaLambda0_2(M); // for f=2 with flip
  
  // Compute lngammaalpha_jkmsf
  for (m = 0; m < M; m++) {
    Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x La matrix
    arma::Cube<double> LngammaLambda0_matrix_1(K, Lx[m],S); // for f=1 (no flip)
    arma::Cube<double> LngammaLambda0_matrix_2(K, Lx[m],S); // for f=2 (with flip)
    for(k = 0; k < K; k++){
       for(s = 0; s < S; s++){
         for(j = 0; j < Lx[m]; j++){
          LngammaLambda0_matrix_1(k, j, s) = gsl_sf_lngamma( exp( Lambda_matrix(k, j+s ) ) );
          LngammaLambda0_matrix_2(k, j, s) = gsl_sf_lngamma( exp( Lambda_matrix(k, La[m]-S+s-j ) ) );
        } //j
      } //s
    } //k
    LngammaLambda0_1(m) = LngammaLambda0_matrix_1;
    LngammaLambda0_2(m) = LngammaLambda0_matrix_2;
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
      arma::Cube<double> LngammaLambda0_matrix_1 = LngammaLambda0_1(m); //K x Lx x S
      arma::Cube<double> LngammaLambda0_matrix_2 = LngammaLambda0_2(m); //K x Lx x S
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
            
            lambda_ks[l]=lambda_matrix(k, s+l); //lambda_{ks}^m
          }
          // LngammaLambda0_matrix K x Lx x S
          double dNegLogEviI = neg_log_evidence_i(data_row, lambda_ks,
                                                  as<Rcpp::NumericVector>( Rcpp::wrap(LngammaLambda0_matrix_1.slice(s).row(k)) ) );
          
          if (dNegLogEviI < offset[m]) //1.0e9 Can this explode?
            offset[m] = dNegLogEviI; // offset[m] is the smallest dNegLogEviI over k and s
          evidence_matrix_1(m, k, s) = dNegLogEviI; //M x K x S
          
          dNegLogEviI = neg_log_evidence_i(data_row, Rcpp::rev( lambda_ks ), //alpha_{La-S+s-j+1}
                                           as<Rcpp::NumericVector>( Rcpp::wrap(LngammaLambda0_matrix_2.slice(s).row(k)) ) );
          if (dNegLogEviI < offset[m]){ //1.0e9 Can this explode?
            offset[m] = dNegLogEviI; // offset[m] is the smallest dNegLogEviI over k and f
          }
          evidence_matrix_2(m, k, s) = dNegLogEviI; //M x K x 2
          
          
          
        } //s
      } //over K
    } //over M
    
    for (k = 0; k < K; k++) {
      //Rcpp::Rcout << "k: "<<k << "\n";
      for(s = 0; s < S; s++){
        //Rcpp::Rcout << "s: "<<s << "\n";
        Z_1(k, s, i) = 0.0; //K x S x N
        Z_2(k, s, i) = 0.0; //K x S x N
        for (m = 0; m < M; m++) {
          //Rcpp::Rcout << "m: "<<m << "\n";
          Z_1(k, s, i) += (-(evidence_matrix_1(m, k, s) - offset[m])); //why offset is substracted?? For numerical reasons?
          Z_2(k, s, i) += (-(evidence_matrix_2(m, k, s) - offset[m])); //why offset is substracted?? For numerical reasons?
        }
        
        Z_1(k, s, i) = W[k]*xi[s] *zeta[0]* exp(Z_1(k, s, i)); //back from logaritm, multiply by the shift probs i.e. \xi_k
        Z_2(k, s, i) = W[k]*xi[s] *zeta[1]* exp(Z_2 (k, s, i)); //back from logaritm, multiply by the shift probs i.e. \xi_k
        dSum += Z_1(k, s, i); //normalization constant
        dSum += Z_2(k, s, i); //normalization constant
      } //s
    } // k
    
    for (k = 0; k < K; k++){
      for(s = 0; s < S; s++){
        Z_1(k, s, i) /= dSum;
        Z_2(k, s, i) /= dSum;
      }
    }
  } // i
  
  Z[0] = Z_1;
  Z[1] = Z_2;
  return Z; //K x N matrix
}

// [[Rcpp::export]]
double neg_log_likelihood_shift_flip(Rcpp::NumericVector W, Rcpp::NumericVector xi, Rcpp::NumericVector zeta, Rcpp::List Lambda,
                                Rcpp::List data, double eta, double nu,
                                double etah, double nuh, int S)
{
  
  Rcpp::IntegerMatrix temp = as<Rcpp::IntegerMatrix>(data[0]); // N x Lx
  const int N = temp.nrow();
  const int K = W.length();
  const int M = data.size();
  int i, j, k, m, s, f;
  
  Rcpp::NumericVector Lx(M);
  Rcpp::NumericVector La(M);
  for (m = 0; m < M; m++) {
    temp = as<Rcpp::IntegerMatrix>(data[m]);
    Lx[m] = temp.ncol();
    La[m] = Lx[m] + S - 1; //assume all Lx are the same
  }
  
  
  Rcpp::List LngammaLambda0_1(M); // for f=1 no flip
  Rcpp::List LngammaLambda0_2(M); // for f=2 with flip
  
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
    arma::Cube<double> LngammaLambda0_matrix_1(K, Lx[m],S); // lng (alpha_jk)
    arma::Cube<double> LngammaLambda0_matrix_2(K, Lx[m],S); 
    Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x La matrices
    for (k = 0; k < K; k++){
       for(s = 0; s < S; s++){
         dSumAlphaK[s] = 0.0; // separate sum for each k,m,s
         LogBAlpha(m, k, s) = 0.0;
         for (j = 0; j < Lx[m]; j++){
            double dAlpha_1 = exp(Lambda_matrix(k, j + s )); // \ alpha_jks
            double lngammaAlpha_1 = gsl_sf_lngamma(dAlpha_1); // lng (alpha_jk)
            LngammaLambda0_matrix_1(k, j, s) = lngammaAlpha_1; // // lng (alpha_jk)
            LngammaLambda0_matrix_2(k, j, s) = gsl_sf_lngamma( exp( Lambda_matrix(k, La[m]-S+s-j ) ) );
            dSumAlphaK[s] += dAlpha_1; // \sum_j^L \alpha_jk
            LogBAlpha(m, k, s) += lngammaAlpha_1; // LogBAlpha_1
         } //over j
         LogBAlpha(m, k,s) -= gsl_sf_lngamma(dSumAlphaK[s]); // -LogBAlpha_2
       } //over s
    } // over k
    LngammaLambda0_1(m) = LngammaLambda0_matrix_1; // list of length M
    LngammaLambda0_2(m) = LngammaLambda0_matrix_2; // list of length M
  } // over m
  
  for (i = 0; i < N; i++) {
    double FinaldProb = 0.0;
    arma::Cube<double> LogStore_1(M, K, S);
    arma::Cube<double> LogStore_2(M, K, S);
    Rcpp::NumericVector offset(M);
    // Computes the log p(x_i^m|, z_ik=1, theta) for k and m, results in matrix LogStore
    for (m = 0; m < M; m++){
      offset[m] = -BIG_DBL;
      double dSum = 0.0; // \sum_j^L x_ij
      double dFactor = 0.0; // \sum_j^L lng( x_ij +1) - lngamma( \sum_j^L x_ij +1)
      
      Rcpp::IntegerMatrix data_matrix = as<Rcpp::IntegerMatrix>(data[m]);
      Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]);
      arma::Cube<double> LngammaLambda0_matrix_1 = LngammaLambda0_1(m); //size K x Lx x S
      arma::Cube<double> LngammaLambda0_matrix_2 = LngammaLambda0_2(m); //size K x Lx x S
      
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
        Rcpp::NumericVector dLogBAlphaN_1(S); // lng( \alpha_jk + x_ij ) - lng( \sum_j^L \alpha_jk + x_ij )
        Rcpp::NumericVector dLogBAlphaN_2(S);
        for(s = 0; s < S; s++){
          dSumAlphaKN[s]=0.0;
          dLogBAlphaN_1[s]=0.0;
          dLogBAlphaN_2[s]=0.0;
          for (j = 0; j < Lx[m]; j++) {
            int countN = data_matrix(i, j); // x_ij
            double dAlphaN_1 = exp(Lambda_matrix(k, j+s)) + countN; // \alpha_jk + x_ij
            double dAlphaN_2 = exp(Lambda_matrix(k, La[m] - S + s -j )) + countN;
            dSumAlphaKN[s] += dAlphaN_1; // \sum_j^L \alpha_jk + x_ij
            dLogBAlphaN_1[s] += countN ? gsl_sf_lngamma(dAlphaN_1) : LngammaLambda0_matrix_1(k, j, s); // dLogBAlphaN_1[m,k,s]
            dLogBAlphaN_2[s] += countN ? gsl_sf_lngamma(dAlphaN_2) : LngammaLambda0_matrix_2(k, j , s);
          }
          dLogBAlphaN_1[s] -= gsl_sf_lngamma(dSumAlphaKN[s]); // // dLogBAlphaN_2[m,k,s]
          dLogBAlphaN_2[s] -= gsl_sf_lngamma(dSumAlphaKN[s]);
          // LogStore (m,k) seems to be p(\mathbf{x}_i|\theta) 
          // LogStore= lng( \alpha_jk + x_ij ) - lng( \sum_j^L \alpha_jk + x_ij )
          // -\sum_j^L lng (alpha_jk) + lngamma( \sum_j^L \alpha_jk) - \sum_j^L lng( x_ij +1)
          LogStore_1(m, k, s) = dLogBAlphaN_1[s] - LogBAlpha(m, k,s) - dFactor; // the positive? log marginal likelihood of sample x_i^m for cluster k
          if (LogStore_1(m, k,s) > offset(m))
            offset(m) = LogStore_1(m, k,s); //offset will be the largest of LogStore(m,:,:), or should it be largest of LogStore(m,:,k)??
          LogStore_2(m, k, s) = dLogBAlphaN_2[s] - LogBAlpha(m, k,s) - dFactor; // the positive? log marginal likelihood of sample x_i^m for cluster k
          if (LogStore_2(m, k,s) > offset(m))
            offset(m) = LogStore_2(m, k,s); //offset will be the largest of LogStore(m,:,:), or should it be largest of LogStore(m,:,k)??
          
        } //over S
      } //over K
    } //over M
    
    Rcpp::NumericVector dProb(K);
    Rcpp::NumericMatrix dProb_F(K, S);
    for(k = 0; k < K; k++){
      dProb[k]=0.0;
      for (s = 0; s < S; s++) {
        dProb_F(k,s)=0.0;
        
        dProb_F(k,s) +=zeta[0]*exp( Rcpp::sum( as<Rcpp::NumericVector>(Rcpp::wrap( LogStore_1.slice(s).col(k) ) )) -sum(offset) );
        dProb_F(k,s) +=zeta[1]*exp( Rcpp::sum( as<Rcpp::NumericVector>(Rcpp::wrap( LogStore_2.slice(s).col(k) ) )) -sum(offset) );
        
        dProb[k] += xi[s]*dProb_F(k,s);
        
      }
      
      // normalize mixing weights
            
      double piK = W[k]/Rcpp::sum(W);
      FinaldProb += piK*dProb[k]; //sum of column k over m rows
    } // k
    
    dRet += log(FinaldProb)+Rcpp::sum(offset); //logarithm of the normalization term, sum over all n
    
    
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



/*Model with flipping only */

// [[Rcpp::export]]
double neg_log_evidence_lambda_pi_flip(Rcpp::NumericVector lambda, Rcpp::List lparams){
  int i, j, f;
  
  Rcpp::IntegerMatrix aanX = as<Rcpp::IntegerMatrix>(lparams["data"]); // N x L_x
  Rcpp::NumericMatrix adPi = as<Rcpp::NumericMatrix>(lparams["pi"]);   // 2xN, these are z_{km}
  double GAMMA_ITA = as<double>(lparams["eta"]);
  double GAMMA_NU = as<double>(lparams["nu"]);
  double GAMMA_ITA_H = as<double>(lparams["etah"]);
  double GAMMA_NU_H = as<double>(lparams["nuh"]);
  
  const int L = aanX.ncol(); // L_x
  const int N = aanX.nrow();
  const int F = adPi.nrow();
  
  /*dLogE collects the terms \sum_{n=1}^N \sum_{j=1}^L E[z_i] \log \gamma (x_ij+ alpha_j)
   and \sum_{n=1}^N E[z_i]* lng( \sum_j(x_ij+alpha-j) )*/
  double dLogE = 0.0;
  
  double dLogEAlpha = 0.0; // \log \gamma ( \sum_{j=1}^L  \alpha_{j} )
  double dSumAlpha = 0.0; // sum of alpha \sum_{j=1}^L \alpha_{j}
  double dSumLambda = 0.0; // sum of lambda?

  double dWeight = 0.0; // \sum_{n=1}^N E(z_i)
  
  Rcpp::NumericVector dWeight_f(N);
  
  Rcpp::NumericVector adSumAlphaN(N); // \sum_{j}^L \alpha_{j}+x_{ij}
  
  for (i = 0; i < N; i++) {
    adSumAlphaN[i] = 0.0;
    dWeight_f[i] = 0.0;
    for (f = 0; f < F; f++) {
      dWeight += adPi(f,i);
      dWeight_f[i] += adPi(f,i);
    }
    
  }
  
  
  for (j = 0; j < L; j++) {
    const double dLambda = lambda[j];
    const double dAlpha = exp(dLambda);
    /* compute the logarithm of the Gamma function,
     dAlpha can not be zero or negative.
     Function computed using the real Lanczos method */
    dLogEAlpha += gsl_sf_lngamma(dAlpha); //dLogEAlpha_1 \sum_{j=1}^L \log \gamma (\alpha_{j} )
    dSumLambda += dLambda;
    dSumAlpha += dAlpha;
    const double lngammaAlpha0_1 = gsl_sf_lngamma(dAlpha); //lngamma of \alpha_{j}
    const double lngammaAlpha0_2 = gsl_sf_lngamma( exp( lambda[ L - 1 - j ] )  ); //lngamma of \alpha_{L-j+1}
    for (i = 0; i < N; i++) {
      const double dN = aanX(i, j); // x_{ij}
      const double dAlphaN_1 = dAlpha + dN; // \alpha_{j}+x_{ij}
      const double dAlphaN_2 = exp( lambda[ L - 1 - j ] ) + dN; // \alpha_{L-j+1}+x_{ij}
      const double lngammaAlphaN_1 = dN ? gsl_sf_lngamma(dAlphaN_1) : lngammaAlpha0_1; //if dN exists or is non-zero, compute lnGamma(dAlphaN), else compute lnGamma(DaLpha)
      const double lngammaAlphaN_2 = dN ? gsl_sf_lngamma(dAlphaN_2) : lngammaAlpha0_2; 
      adSumAlphaN[i] += dAlphaN_1; // \sum_{j}^L \alpha_{j}+x_{ij} , weight by pi
      
      dLogE -= adPi(0,i) * lngammaAlphaN_1; //dlogE_1 weight by pi, -\sum_i E[z_i] *\sum_j lngamma(x_{ij}+alpha_j)
      dLogE -= adPi(1,i) * lngammaAlphaN_2; 
    } //i
  } //j
  dLogEAlpha -= gsl_sf_lngamma(dSumAlpha);// dLogEAlpha_2 \sum_{j=1}^L \log \gamma (\alpha_{j} ) -\log \gamma ( \sum_{j=1}^L  \alpha_{j} )
  
  for(i = 0; i < N; i++){
    dLogE += dWeight_f[i] * gsl_sf_lngamma(adSumAlphaN[i]); //dLogE_2 \sum_{n=1}^N E[z_i]* lngamma( \sum_j(x_ij+alpha_j) )
  }
    
  
  double reg_term;
  
  if ((GAMMA_ITA_H==0) && (GAMMA_NU_H==0)) {
    reg_term = 0.0;
  } else {
    double hk = REG_H_LAMBDA(lambda); // computes the value of the regularization term
    
    reg_term = GAMMA_NU_H*hk - (GAMMA_ITA_H-1)*log(hk); //This was nuh *hg-etah*log (hk), it is now corrected to (GAMMA_ITA_H-1)!!!
  }
  
  return dLogE + dWeight*dLogEAlpha + // complete data likelihood term
    GAMMA_NU*dSumAlpha - (GAMMA_ITA - 1) * dSumLambda + reg_term;//prior term, ITA corrected to ITA-1 !!!
  

}


// [[Rcpp::export]]
Rcpp::NumericVector neg_log_derive_evidence_lambda_pi_flip(Rcpp::NumericVector ptLambda,
                                                      Rcpp::List lparams)
{
  Rcpp::IntegerMatrix aanX = as<Rcpp::IntegerMatrix>(lparams["data"]); // N x L
  Rcpp::NumericMatrix adPi = as<Rcpp::NumericMatrix>(lparams["pi"]); // 2xN
  double GAMMA_ITA = as<double>(lparams["eta"]);
  double GAMMA_NU = as<double>(lparams["nu"]);
  double GAMMA_ITA_H = as<double>(lparams["etah"]);
  double GAMMA_NU_H = as<double>(lparams["nuh"]);
  Rcpp::List hkm = as<Rcpp::List>(lparams["hkm"]);
  Rcpp::IntegerVector hkm_index = as<Rcpp::IntegerVector>(lparams["hkm_index"]);
  
  Rcpp::List gradient = as<Rcpp::List>(lparams["gradient"]);
  Rcpp::IntegerVector gradient_index = as<Rcpp::IntegerVector>(lparams["gradient_index"]);
  
  int i, j, f;
  const int L = aanX.ncol();
  const int N = aanX.nrow();
  const int F = adPi.nrow();
  
  Rcpp::NumericVector g(L);
  Rcpp::NumericVector adDeriv(L); //derivative for each j
  Rcpp::NumericVector adStore(N); // \sum_j^L x_ij + \sum_j^L \alpha_j
  Rcpp::NumericVector adAlpha(L);
  double dSumStore = 0.0; // \sum_n^N E[z_i]* psi( \sum_j^L x_ij + \sum_j^L \alpha_j )
  double dStore = 0.0; // sum of alpha over j
  double dWeight = 0.0; // sum of z_i over i/N
  Rcpp::NumericVector dWeight_f(N);
  
  for (i = 0; i < N; i++) {  
    adStore[i] = 0.0;
    dWeight_f[i] = 0.0;
    for (f = 0; f < F; f++) {
      dWeight += adPi(f,i);
      dWeight_f[i] += adPi(f,i);
    }
  }
  
  // There are some unnecessary terms here
  for (j = 0; j < L; j++) {
    adAlpha[j] = exp(ptLambda[j]);
    dStore += adAlpha[j];
    adDeriv[j] = dWeight* gsl_sf_psi(adAlpha[j]); //adDeriv_1[j]
    double alphaS0 = gsl_sf_psi(adAlpha[j]);
   
    for (i = 0; i < N; i++) {
      int dN_1 = aanX(i, j);
      int dN_2 = aanX(i, L - 1 - j);
      double dAlphaN_1 = adAlpha[j] + dN_1;
      double dAlphaN_2 = adAlpha[j] + dN_2;
      double psiAlphaN_1 = dN_1 ? gsl_sf_psi(dAlphaN_1) : alphaS0;
      double psiAlphaN_2 = dN_2 ? gsl_sf_psi(dAlphaN_2) : alphaS0;
      adDeriv[j] -= adPi(0,i)*psiAlphaN_1; //adDeriv_2[j] f=1
      adDeriv[j] -= adPi(1,i)*psiAlphaN_2; //adDeriv_2[j] f=2
      adStore[i] += dAlphaN_1; //  \sum_j^L x_ij + \sum_j^L \alpha_j
    }
  }
  
  for (i = 0; i < N; i++){
    dSumStore += dWeight_f[i] * gsl_sf_psi(adStore[i]);
  }
  dStore = dWeight * gsl_sf_psi(dStore); //dStore_2
  
  double hk = REG_H_LAMBDA(ptLambda);
  
  if (hkm_index[0] < hkm.size()) {
    hkm[hkm_index[0]] = hk;
  }
  else{
    Rprintf("hkm_index exceeds hkm vector length!!! hkm_index: %i, hkm.length: %i\n",
            hkm_index[0], hkm.size());
  }
  hkm_index[0] += 1;
  
  for (j = 0; j < L; j++) {
    
    double reg_term;
    
    if ((GAMMA_ITA_H==0) && (GAMMA_NU_H==0)) {
      reg_term = 0.0;
    } else {
      double gjk = REG_DERIV_LAMBDA_G(ptLambda, j); //deriv. of h_kj^(m) wrt \alpha_kj^(m)
      reg_term = GAMMA_NU_H*gjk - (GAMMA_ITA_H-1)*gjk/hk;
    }
    //This is now corrected, GAMMA_ITA converted to GAMMA_ITA-1
    double value = adAlpha[j] *
      (GAMMA_NU + adDeriv[j] - dStore + dSumStore + reg_term) - ( GAMMA_ITA -1 ); // should be (GAMMA_ITA -1)
    
    g[j] = value;
  }
  
  double gnorm=0.0;
  
  for (j = 0; j < L; j++) {
    gnorm += pow(g[j], 2.0);
  }
  gnorm=sqrt(gnorm);
  
  //Rprintf("Gradient norm: %f\n", gnorm);
  
  if (gradient_index[0] < gradient.size()) {
    gradient[gradient_index[0]] = gnorm;
  }
  else{
    Rprintf("gradient_index exceeds gradient vector length!!! gradient_index: %i, gradient.length: %i\n",
            gradient_index[0], gradient.size());
  }
  gradient_index[0] += 1;
  
  
  return g;
}

// /* DONE Computes the E-step, i.e. the posterior probabilities of the cluster labels
//  * Z=Ez K x 2 x N matrix (from the previous E-step, these are not used?)
//  * data=binned.data list of M matrixes of size N x Lx
//  * W=weights vector of length K, \pi_k^{old}
//  * lambda=lambda list of M, each K times L matrix
//  *
//  * returns K x S x N matrix Z
//  *
//  * 
//  */

// [[Rcpp::export]]
arma::Cube<double> calc_z_flip(arma::Cube<double> Z, Rcpp::List data,
                                Rcpp::NumericVector W, Rcpp::NumericVector zeta, Rcpp::List Lambda)
{
  // Z is K x 2 x N
  int i, j, k, m, f;
  Rcpp::IntegerMatrix temp = as<Rcpp::IntegerMatrix>(data[0]);
  const int N = temp.nrow();
  const int K = W.length(); // number of clusters
  const int M = data.size();
  const int F = zeta.length();
  
  arma::Cube<double> evidence_matrix(M, K, F); 
  
  Rcpp::NumericVector Lx(M);
  for (m = 0; m < M ; m++) {
    temp = as<Rcpp::IntegerMatrix>(data[m]); //N x L matrix
    Lx[m] = temp.ncol(); //L
  }
  
  Rcpp::List LngammaLambda0(M); // lngamma ( \alpha_{jk}^{(m)}  )
  
  // Compute lngammaalpha_jkms
  for (m = 0; m < M; m++) {
    Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x La matrix
    arma::Cube<double> LngammaLambda0_matrix(K, Lx[m],F); // lngamma ( \alpha_{jk}  )
    for(k = 0; k < K; k++){
      for(j = 0; j < Lx[m]; j++){
        const double dAlpha_1 = exp( Lambda_matrix(k, j) );
        const double dAlpha_2 = exp( Lambda_matrix(k, Lx[m] - 1 -j ) );
        LngammaLambda0_matrix(k, j, 0) = gsl_sf_lngamma(dAlpha_1);
        LngammaLambda0_matrix(k, j, 1) = gsl_sf_lngamma(dAlpha_2);
      } //j
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
      arma::Cube<double> LngammaLambda0_matrix = LngammaLambda0(m); //K x Lx x 2
      Rcpp::IntegerVector data_row = data_matrix(i, _); // one row of X of length Lx
      for (k = 0; k < K; k++) {

        //Computes the logarithm of -p(\mathbf{x}_i|z_{ik}=1,\bm{\theta}) but of only those
        // terms that depend on \alpha
        // I.e. computes the negative logarithm of the unnormalized DirichletMultinomial pdf value
        //logarithm due to numerical purposes
        // computes the - log DirMulti(x_i^m | alpha_k^m), terms not depending on alpha excluded
        
        // LngammaLambda0_matrix K x Lx x 2
        double dNegLogEviI = neg_log_evidence_i(data_row, lambda_matrix(k,_),
                                                  as<Rcpp::NumericVector>( Rcpp::wrap(LngammaLambda0_matrix.slice(0).row(k)) ) );
        
        if (dNegLogEviI < offset[m]){ //1.0e9 Can this explode?
            offset[m] = dNegLogEviI; // offset[m] is the smallest dNegLogEviI over k and f
        }
        evidence_matrix(m, k, 0) = dNegLogEviI; //M x K x 2
      
        dNegLogEviI = neg_log_evidence_i(data_row, Rcpp::rev( lambda_matrix(k,_) ),
                                                as<Rcpp::NumericVector>( Rcpp::wrap(LngammaLambda0_matrix.slice(1).row(k)) ) );
        if (dNegLogEviI < offset[m]){ //1.0e9 Can this explode?
          offset[m] = dNegLogEviI; // offset[m] is the smallest dNegLogEviI over k and f
        }
        evidence_matrix(m, k, 1) = dNegLogEviI; //M x K x 2
      
              
      } //over K
    } //over M
    
    for (k = 0; k < K; k++) {
        Z(k, 0, i) = 0.0; //K x 2 x N
        Z(k, 1, i) = 0.0; //K x 2 x N
        for (m = 0; m < M; m++) {
          
          Z(k, 0, i) += (-(evidence_matrix(m, k, 0) - offset[m])); 
          Z(k, 1, i) += (-(evidence_matrix(m, k, 1) - offset[m])); 
        }
        Z(k, 0, i) = W[k]*zeta[0] * exp(Z(k, 0, i)); 
        Z(k, 1, i) = W[k]*zeta[1] * exp(Z(k, 1, i));
        dSum += Z(k, 0, i); //normalization constant
        dSum += Z(k, 1, i); //normalization constant
     
    } // k
    
    for (k = 0; k < K; k++){
      
        Z(k, 0, i) /= dSum;
        Z(k, 1, i) /= dSum;
    }
  } // i
  return Z; //K x 2 x N matrix
}

// [[Rcpp::export]]
double neg_log_likelihood_flip(Rcpp::NumericVector W, Rcpp::NumericVector zeta, Rcpp::List Lambda,
                          Rcpp::List data, double eta, double nu,
                          double etah, double nuh)
{
  Rcpp::IntegerMatrix temp = as<Rcpp::IntegerMatrix>(data[0]); // N x Lx
  const int N = temp.nrow();
  const int K = W.length();
  const int M = data.size();
  const int F = zeta.length();
  int i, j, k, m, f;
  
  Rcpp::NumericVector L(M);
  for (m = 0; m < M; m++) {
    temp = as<Rcpp::IntegerMatrix>(data[m]);
    L[m] = temp.ncol();
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
  
  Rcpp::NumericMatrix LogBAlpha(M, K); // // \sum_j^L lng (alpha_jk) -lngamma( \sum_j^L \alpha_jk)
  
  for (m = 0; m < M; m++) {
    Rcpp::NumericMatrix LngammaLambda0_matrix(K, L[m]); // lng (alpha_jk)
    Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x La matrices
    
    for (k = 0; k < K; k++){
      double dSumAlphaK = 0.0;
      LogBAlpha(m, k) = 0.0;
      for (j = 0; j < L[m]; j++){
        double dAlpha = exp(Lambda_matrix(k, j)); // \ alpha_jk
        double lngammaAlpha = gsl_sf_lngamma(dAlpha); // lng (alpha_jk)
        LngammaLambda0_matrix(k, j) = lngammaAlpha; // // lng (alpha_jk)
        
        dSumAlphaK += dAlpha; // \sum_j^L \alpha_jk
        LogBAlpha(m, k) += lngammaAlpha; // LogBAlpha_1 \sum_j^L lng (alpha_jk)
      } //over j
      LogBAlpha(m, k) -= gsl_sf_lngamma(dSumAlphaK); // LogBAlpha_2 \sum_j^L lng (alpha_jk) -lngamma( \sum_j^L \alpha_jk)
    } // over k
    LngammaLambda0[m] = LngammaLambda0_matrix; // list of length M
  } // over m
  for (i = 0; i < N; i++) {
    double dProb = 0.0;
    arma::Cube<double> LogStore(M, K, 2);
    Rcpp::NumericVector offset(M);
    // Computes the log p(x_i^m|, z_ik=1, theta) for k and m, results in matrix LogStore
    for (m = 0; m < M; m++) {
      double dSum = 0.0; // \sum_j^L x_ij
      double dFactor = 0.0; // \sum_j^L lng( x_ij +1) - lngamma( \sum_j^L x_ij +1)
      offset(m) = -BIG_DBL;
      
      Rcpp::IntegerMatrix data_matrix = as<Rcpp::IntegerMatrix>(data[m]);
      Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]);
      Rcpp::NumericMatrix LngammaLambda0_matrix = as<Rcpp::NumericMatrix>(LngammaLambda0[m]); // lng (alpha_jk) K x L matrix
      
      for (j = 0; j < L[m]; j++) {
        dSum += data_matrix(i, j); // dSum_m, this is J_i^m
        // dFactor does not depend on parameters, why computed?
        // dFactor_1
        dFactor += gsl_sf_lngamma(data_matrix(i, j) + 1.0); // dFactor_1 \sum_j^L lng( x_ij +1),
      }
      // dFactor does not depend on parameters, why computed?
      dFactor -= gsl_sf_lngamma(dSum + 1.0);  //dFactor_2 \sum_j^L lng( x_ij +1) - lngamma( \sum_j^L x_ij +1)
      
      for (k = 0; k < K; k++) {
        double dSumAlphaKN = 0.0; // \sum_j^L \alpha_jk + x_ij
        double dLogBAlphaN_1 = 0.0; // lng( \alpha_jk + x_ij ) - lng( \sum_j^L \alpha_jk + x_ij )
        double dLogBAlphaN_2 = 0.0; 
        for (j = 0; j < L[m]; j++) {
          int countN = data_matrix(i, j); // x_ij
          double dAlphaN_1 = exp(Lambda_matrix(k, j)) + countN; // \alpha_jk + x_ij
          double dAlphaN_2 = exp(Lambda_matrix(k, L[m] - 1 - j) ) + countN; // \alpha_jk + x_ij
          dSumAlphaKN += dAlphaN_1; // \sum_j^L \alpha_jk + x_ij
          dLogBAlphaN_1 += countN ? gsl_sf_lngamma(dAlphaN_1) : LngammaLambda0_matrix(k, j); //dLogBAlphaN_1_1 lng( \alpha_jk + x_ij )
          dLogBAlphaN_2 += countN ? gsl_sf_lngamma(dAlphaN_2) : LngammaLambda0_matrix(k, L[m] - 1 - j ); //dLogBAlphaN_2_1 lng( \alpha_jk + x_ij )
        }
        dLogBAlphaN_1 -= gsl_sf_lngamma(dSumAlphaKN); // dLogBAlphaN_1_2 lng( \alpha_jk + x_ij ) - lng( \sum_j^L \alpha_jk + x_ij )
        dLogBAlphaN_2 -= gsl_sf_lngamma(dSumAlphaKN); 
        //LogStore (m,k) seems to be p(\mathbf{x}_i|\theta) but without term \Gamma(J_i^{m}+1)
        // LogStore= lng( \alpha_jk + x_ij ) - lng( \sum_j^L \alpha_jk + x_ij )
        // -\sum_j^L lng (alpha_jk) + lngamma( \sum_j^L \alpha_jk) - \sum_j^L lng( x_ij +1)
        
        LogStore(m, k,0) = dLogBAlphaN_1 - LogBAlpha(m, k) - dFactor; // the positive? log marginal likelihood of sample x_i^m for cluster k
        if (LogStore(m, k,0 ) > offset(m))
          offset(m) = LogStore(m, k, 0); //offset will be the largest of LogStore(m,1:K, 1:2) 
        
        LogStore(m, k,1) = dLogBAlphaN_2 - LogBAlpha(m, k) - dFactor; // the positive? log marginal likelihood of sample x_i^m for cluster k
        if (LogStore(m, k,1 ) > offset(m))
          offset(m) = LogStore(m, k, 1); //offset will be the largest of LogStore(m,1:K, 1:2) 
        
        
      } //over K
    } //over M
    
    //offset?
    Rcpp::NumericVector dProb_F(K);
    for (k = 0; k < K; k++) {
      dProb_F[k]=0.0;
      dProb_F[k] += zeta[0]*exp( Rcpp::sum( as<Rcpp::NumericVector>(Rcpp::wrap( LogStore.slice(0).col(k) ) )) - sum(offset) );
      dProb_F[k] += zeta[1]*exp( Rcpp::sum( as<Rcpp::NumericVector>(Rcpp::wrap( LogStore.slice(1).col(k) ) )) - sum(offset) );
      // normalize mixing weights
      double piK = W[k]/Rcpp::sum(W); 
      dProb += piK*dProb_F[k];
    }
    dRet += log(dProb)+sum(offset); //logarithm of the normalization term, sum over all n
  } //over N
  
  dL5 = -sum(L) * K * gsl_sf_lngamma(eta); // -M*L*K*lng(eta)
  dL6 = eta * K * sum(L) * log(nu);  //eta*K*M*L*log(nu)
  
  if ((etah!=0) || (nuh!=0)) {
    regterm1 = M * K * (etah * log(nuh) - gsl_sf_lngamma(etah));
  }
  
  for (m = 0; m < M; m++) {
    Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x L matrix
    
    for (k = 0; k < K; k++) {
      
      if ((etah!=0) || (nuh!=0)) {
        //exp(Lambda_matrix(i, _) are alpha_k, all L elements
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
double neg_log_evidence_lambda_pi_shift(Rcpp::NumericVector lambda, Rcpp::List lparams){
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
  Rcpp::NumericVector dAlpha(La); 
  Rcpp::NumericMatrix adSumAlphaN(S,N); 
  Rcpp::NumericVector dSumAlpha_s(S); 
  Rcpp::NumericVector dLogEAlpha(S); 
  //Compute dWeight_[s]
  for(s = 0; s < S; s++){
    dSumAlpha_s[s]=0.0;
    dWeight[s]=0.0;
    dLogEAlpha[s]=0.0;
    for (i = 0; i < N; i++) {
      adSumAlphaN(s,i)=0.0;
      dWeight[s] += adPi(s,i); //sum over i 
    }
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
        
        const double dN = aanX(i, j); // x_{ij}
        const double dAlphaN_1 = dAlpha[j+s] + dN; // depends on i, j and s
        const double lngammaAlphaN_1 = dN ? gsl_sf_lngamma(dAlphaN_1) : lngammaAlpha0_1; //if dN exists or is non-zero, compute lnGamma(dAlphaN), else compute lnGamma(DaLpha)
       
        adSumAlphaN(s,i) += dAlphaN_1; // \sum_{j}^Lx \alpha_{j}+x_{ij}
        dLogE -= adPi(s,i) * lngammaAlphaN_1; //dLogE_1 -\sum_i \sum_s E[z_isf] *\sum_j lngamma(x_{ij}+alpha_j), 
        
      } //i
    } //j
    
    dLogEAlpha_final += dLogEAlpha[s] * dWeight[s]; //dLogEAlpha1_[s] // \sum_{j=1}^S \log \gamma (\alpha_{j} )
    dLogEAlpha_final -= gsl_sf_lngamma(dSumAlpha_s[s]) * dWeight[s];  //dLogEAlpha_2[s] -\log \gamma ( \sum_{j=1}^S  \alpha_{j} )
  } //s
  
  
  for(i = 0; i < N; i++){
    for(s = 0; s < S; s++){
      dLogE += adPi(s,i) * gsl_sf_lngamma(adSumAlphaN(s,i)); //dLogE_2 \sum_{n=1}^N E[z_i]* lngamma( \sum_j(x_ij+alpha_j) )
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
                                                                 Rcpp::List lparams){
  Rcpp::IntegerMatrix aanX = as<Rcpp::IntegerMatrix>(lparams["data"]); // N x Lx
  Rcpp::NumericMatrix adPi = as<Rcpp::NumericMatrix>(lparams["pi"]); // S x N
  double GAMMA_ITA = as<double>(lparams["eta"]);
  double GAMMA_NU = as<double>(lparams["nu"]);
  double GAMMA_ITA_H = as<double>(lparams["etah"]);
  double GAMMA_NU_H = as<double>(lparams["nuh"]);
  Rcpp::List hkm = as<List>(lparams["hkm"]);
  Rcpp::IntegerVector hkm_index = as<Rcpp::IntegerVector>(lparams["hkm_index"]);
  
  Rcpp::List gradient = as<Rcpp::List>(lparams["gradient"]);
  Rcpp::IntegerVector gradient_index = as<Rcpp::IntegerVector>(lparams["gradient_index"]);
  
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
    dWeight[s]=0.0;
    for (i = 0; i < N; i++) {
      adStore(s,i) = 0.0;
      
      dWeight[s] += adPi(s,i); //sum over i a
    }
   
  }
  
  for(j = 0; j < La; j++){
    adAlpha[j] = exp(ptLambda[j]);
    dSumStore[j] = 0.0;
    dWeight_j[j] = 0.0;
    for(i = 0; i < N; i++){
        for(s = (j - Lx + 1); s < (j + 1); s++){
          if(s > -1 && s < S ){
            dWeight_j[j] += adPi(s,i);
          } //if
          
        } //s
      } //i
  } //j
  
  for(s = 0; s < S; s++){
    for(k = 0; k < Lx; k++ ){
      dStore[s] += adAlpha[k+s];
    }
  }
  
  
  for(i = 0; i < N; i++){
    xStore[i] = 0.0;
    for(k = 0; k < Lx; k++ ){
      xStore[i] += aanX(i, k ); //N x Lx
    }
    for(s = 0; s < S; s++){
      adStore(s,i)=xStore[i]+dStore[s];
    }
  }
  
  
  
  for (j = 0; j < La; j++) {
    double alphaS0 = gsl_sf_psi(adAlpha[j]);
    adDeriv[j] = dWeight_j[j] * alphaS0; // adDeriv_1[j]
    for(s = (j - Lx + 1); s < (j + 1); s++){
      if(s > -1 && s < S){
        for (i = 0; i < N; i++) {
          int dN_1 = aanX(i, j - s ); //c++ indexing starts from 0 -> -1
          double dAlphaN_1 = adAlpha[j] + dN_1;
          double psiAlphaN_1 = dN_1 ? gsl_sf_psi(dAlphaN_1) : alphaS0;
          adDeriv[j] -= adPi(s,i) * psiAlphaN_1; //adDeriv_2[j]
        } //i
      } //if true s
    } //s
   } //j
  
  for(j = 0; j < La; j++){
    for(s = (j - Lx + 1); s < (j + 1); s++){
      if(s > -1 && s < S ){
        for (i = 0; i < N; i++){
          dSumStore[j] += adPi(s,i) * gsl_sf_psi(adStore(s,i));
        }
      }
    }
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
  
  if (hkm_index[0] < hkm.size()){

    hkm[hkm_index[0]] = hk;
  }else{
    Rprintf("hkm_index exceeds hkm vector length!!! hkm_index: %i, hkm.length: %i\n", hkm_index[0], hkm.size());
  } 
  
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
  
  double gnorm=0.0;
  
  for (j = 0; j < La; j++) {
    gnorm += pow(g[j], 2.0);
  }
  gnorm=sqrt(gnorm);
    
  if (gradient_index[0] < gradient.size()) {
    gradient[gradient_index[0]] = gnorm;
  }
  else{
    Rprintf("gradient_index exceeds gradient vector length!!! gradient_index: %i, gradient.length: %i\n",
            gradient_index[0], gradient.size());
  }
  gradient_index[0] += 1;
  
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
   Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x La matrix
    arma::Cube<double> LngammaLambda0_matrix(K, Lx[m],S); // lngamma ( \alpha_{jk}  )
    for(k = 0; k < K; k++){
      for(s = 0; s < S; s++){
        for(j = 0; j < Lx[m]; j++){
          const double dAlpha = exp(Lambda_matrix(k, j+s));
          LngammaLambda0_matrix(k, j, s) = gsl_sf_lngamma(dAlpha);
        } //j
      } //s
    } //k
    LngammaLambda0(m) = Rcpp::wrap(LngammaLambda0_matrix);
  } //m

  for (i = 0; i < N; i ++) {
   
    double dSum = 0.0;
    Rcpp::NumericVector offset(M); //save the smallest negLogEvidence for each cluster(largest absolute value)
    // Compute the evidence matrix, the DirichletMultinomial pdf for given i, m and k
    for (m = 0; m < M; m++) {

      offset[m] = BIG_DBL; //1.0e9

      Rcpp::IntegerMatrix data_matrix = as<Rcpp::IntegerMatrix>(data[m]); // N x Lx matrix
      Rcpp::NumericMatrix lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x La matrix
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
           
            lambda_ks[l]=lambda_matrix(k, s+l); //lambda_{ks}^m
            }
          // LngammaLambda0_matrix K x Lx x S
          double dNegLogEviI = neg_log_evidence_i(data_row, lambda_ks,
                                               as<Rcpp::NumericVector>( Rcpp::wrap(LngammaLambda0_matrix.slice(s).row(k)) ) );
         
          if (dNegLogEviI < offset[m]) //1.0e9 Can this explode?
            offset[m] = dNegLogEviI; // offset[m] is the smallest dNegLogEviI over k and s
          evidence_matrix(m, k, s) = dNegLogEviI; //M x K x S
        
        } //s
      } //over K
    } //over M

    for (k = 0; k < K; k++) {
      for(s = 0; s < S; s++){
        Z(k, s, i) = 0.0; //K x S x N
        for (m = 0; m < M; m++) {
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
  const int K = W.length();
  const int M = data.size();
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
  
  arma::Cube<double> LogBAlpha(M,K,S); 
    
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
            offset(m) = LogStore(m, k,s); //offset will be the largest of LogStore(m,:,:), or should it be largest of LogStore(m,:,k)??
        } //over S
      } //over K
    } //over M
    
    Rcpp::NumericVector dProb_S(K);
    for(k = 0; k < K; k++){
      dProb_S[k]=0.0;
      for (s = 0; s < S; s++) {
        
        dProb_S[k] += xi[s]*exp( sum( as<Rcpp::NumericVector>(Rcpp::wrap( LogStore.slice(s).col(k) ) )) - sum(offset) ); //sum of column k over m rows
      }
      // normalize mixing weights
      double piK = W[k]/Rcpp::sum(W);
      
      dProb += piK*dProb_S[k]; //sum of column k over m rows
      
    } // k
    
    
    dRet += log(dProb)+Rcpp::sum(offset); //logarithm of the normalization term, sum over all n
    
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


