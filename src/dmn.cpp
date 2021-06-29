#include <RcppArmadillo.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>


#include "regul.h"

#define BIG_DBL 1.0e9

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]


// [[Rcpp::export]]
void disable_gsl_error_handler() {
 gsl_set_error_handler_off();
}



/*Soft K-means clustering treats the cluster assignments as probability distributions
over the clusters. There is a connection between Euclidean distance and multivariate normal
models with a fixed covariance, soft K-means can be expressed as a multivariate
normal mixture model*/

// [[Rcpp::export]]
Rcpp::List soft_kmeans(Rcpp::NumericMatrix data,
                 int K,
                 bool verbose=false,
                 bool randomInit=true,
                 Rcpp::NumericMatrix centers=Rcpp::NumericMatrix(), //if not given, initialized randomly
                 double stiffness=50.0,
                 double threshold=1.0e-6,
                 int maxIt=1000,
                 bool rowNorm=true) {
    /*data=kmeans.binned.data, N x (M*S)
     centers=kmeanspp.centers, K x (M*S), can
     rowNorm=F */
    const int S = data.ncol(), N = data.nrow();
    Rcpp::RNGScope scope;

    int i, j, k, iter = 0;

    Rcpp::NumericVector W(K);
    Rcpp::NumericMatrix Z(K, N); //cluster membership probabilities
    Rcpp::List data_dimnames = data.attr("dimnames");
    Z.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector(), data_dimnames[0]); //sample names
    Rcpp::NumericMatrix Mu(K, S);
    Mu.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector(), data_dimnames[1]); //bin names

    double dMaxChange = BIG_DBL; //1.0e9

    if (verbose)
        Rprintf("  Soft kmeans\n");

    Rcpp::NumericMatrix Y(N, S);
    Rcpp::NumericVector Mu_row(S);

    if (rowNorm) {
      for (i = 0; i < N; i++) { //over rows
        int iTotal = 0;
        for (j = 0; j < S; j++) //over columns
          iTotal += data(i, j);  //rowsum

        if (iTotal == 0) iTotal = 1; //workaround for -nan error

        for (j = 0; j < S; j++)
          Y(i,j) = (data(i, j)) / (double)iTotal; //divide the row elements using the rowsum
      }
    } else {
      Y = data;
    }

    Rcpp::Function sample("sample"); //access to R's sample function

    /* initialise */
    if (centers.ncol() == 0) {
      //initialise randomly
      for (i = 0; i < N; i++) {
          if (K == 1)
              k = 0;
          else if (randomInit) {
              k = as<int>(sample(seq(0, K-1), 1)); //random integer btw 0 and K-1
          }
          else
              k = i % K;

          Z(k, i) = 1.0;
      }
    } else {
      //initialize Z matrix based on given centers
      Mu = centers;
      Mu.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector(), data_dimnames[1]);

      for (i = 0; i < N; i++) {
        double dNorm = 0.0, dOffset = BIG_DBL;
        for (k = 0; k < K; k++) {
            Z(k, i) = 0;
            for (j = 0; j < S; j++) {
                const double dDiff = (Mu(k, j) - Y(i, j));
                Z(k, i) += dDiff * dDiff; //Euclidean distance
            }
            Z(k, i) = -stiffness * sqrt(Z(k, i));
            if (Z(k,i) < dOffset)
                dOffset = Z(k, i);
        }
        for (k = 0; k < K; k++) {
            Z(k, i) = exp(Z(k, i) - dOffset);
            dNorm += Z(k, i);
        }
        for (k = 0; k < K; k++) {
            Z(k, i) = Z(k, i) / dNorm; //membership probabilities
        }
      }
    }
    /*If the sum of differences between the previous cluster centers and last cluster
    centers are neglible, stop*/
    while (dMaxChange > threshold && iter < maxIt) {
        /* update mu i.e. the cluster centers*/
        dMaxChange = 0.0;
        for (i = 0; i < K; i++){
            double dNormChange = 0.0;
            W[i] = 0.0;
            for (j = 0; j < N; j++)
                W[i] += Z(i,j);

            for (j = 0; j < S; j++) {
                Mu_row[j] = 0.0;
                for (k = 0; k < N; k++)
                    Mu_row[j] += Z(i, k)*Y(k, j);
            }

            for (j = 0; j < S; j++) {
                double dDiff = 0.0;
                Mu_row[j] /= W[i] ? W[i] : 1;

                dDiff = (Mu_row[j] - Mu(i, j));
                dNormChange += dDiff * dDiff;
                Mu(i,j) = Mu_row[j];
            }
            dNormChange = sqrt(dNormChange);
            if (dNormChange > dMaxChange)
                dMaxChange = dNormChange;
        }

        /* calc distances and update Z */
        for (i = 0; i < N; i++) {
            double dNorm = 0.0, dOffset = BIG_DBL;
            for (k = 0; k < K; k++) {
                Z(k, i) = 0;
                for (j = 0; j < S; j++) {
                    const double dDiff = (Mu(k, j) - Y(i, j));
                    Z(k, i) += dDiff * dDiff;
                }
                Z(k, i) = -stiffness * sqrt(Z(k, i));
                if (Z(k,i) < dOffset)
                    dOffset = Z(k, i);
            }
            for (k = 0; k < K; k++) {
                Z(k, i) = exp(Z(k, i) - dOffset);
                dNorm += Z(k, i);
            }
            for (k = 0; k < K; k++) {
                Z(k, i) = Z(k, i) / dNorm;
            }
        }
        iter++;

        if (verbose && (iter % 10 == 0))
            Rprintf("    iteration %d change %f\n", iter, dMaxChange);
    }

    return Rcpp::List::create(_["centers"] = Mu,
                        _["weights"] = W,
                        _["labels"] = Z); 
}

/*The function below computes the value of the lower bound to be optimized
 * i.e. the whole expected log posterior Q
 * Computes the term only depending on \alpha_{k}^m or \lambda_l^m
 * . Does not consider \pi. Computes the value for one cluster k and datatype m
 * lambda are the current values, lparams contains the current parameter values (E(z_i))
 * List of 8
 $ pi       : Named num [1:1000] 0.38 0.313 0.289 0.346 0.278 ... These are z_i
 $ data     : int [1:1000, 1:44] 0 1 0 1 0 0 0 0 1 1 ...
 $ nu       : num 1
 $ etah     : num 1
 $ nuh      : num 10
 $ hkm      :List of 1001
 ..$ : num 0.000227
 ..$ : num 0.592
 $ hkm_index: int 42
 * CHECKED, eta term converted to eta - 1
 */

// [[Rcpp::export]]
double neg_log_evidence_lambda_pi(Rcpp::NumericVector lambda,  Rcpp::List lparams)
{
    int i, j;

    Rcpp::IntegerMatrix aanX = as<Rcpp::IntegerMatrix>(lparams["data"]); // N x L_x
    Rcpp::NumericVector adPi = as<Rcpp::NumericVector>(lparams["pi"]);   // size N, these are z_{km}
    double GAMMA_ITA = as<double>(lparams["eta"]);
    double GAMMA_NU = as<double>(lparams["nu"]);
    double GAMMA_ITA_H = as<double>(lparams["etah"]);
    double GAMMA_NU_H = as<double>(lparams["nuh"]);
    

    const int L = aanX.ncol(); // L_x
    const int N = aanX.nrow();
    
    /*dLogE collects the terms \sum_{n=1}^N \sum_{j=1}^L E[z_i] \log \gamma (x_ij+ alpha_j)
    and \sum_{n=1}^N E[z_i]* lng( \sum_j(x_ij+alpha-j) )*/
    double dLogE = 0.0;
    
    double dLogEAlpha = 0.0; // \log \gamma ( \sum_{j=1}^L  \alpha_{j} )
    double dSumAlpha = 0.0; // sum of alpha \sum_{j=1}^L \alpha_{j}
    double dSumLambda = 0.0; // sum of lambda?
    double dHk = 0.0;
    double dWeight = 0.0; // \sum_{n=1}^N E(z_i)

    Rcpp::NumericVector adSumAlphaN(N); // \sum_{j}^L \alpha_{j}+x_{ij}

    for (i = 0; i < N; i++) {
        adSumAlphaN[i] = 0.0;
        dWeight += adPi[i]; //sum over z_i
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
        const double lngammaAlpha0 = gsl_sf_lngamma(dAlpha); //lngamma of \alpha_{j}
        for (i = 0; i < N; i++) {
            const double dN = aanX(i, j); // x_{ij}
            const double dAlphaN = dAlpha + dN; // \alpha_{j}+x_{ij}
            const double lngammaAlphaN = dN ? gsl_sf_lngamma(dAlphaN) : lngammaAlpha0; //if dN exists or is non-zero, compute lnGamma(dAlphaN), else compute lnGamma(DaLpha)
            adSumAlphaN[i] += dAlphaN; // \sum_{j}^L \alpha_{j}+x_{ij} , weight by pi
            dLogE -= adPi[i] * lngammaAlphaN; //dlogE_1 weight by pi, -\sum_i E[z_i] *\sum_j lngamma(x_{ij}+alpha_j)
        }
    }
    dLogEAlpha -= gsl_sf_lngamma(dSumAlpha);// dLogEAlpha_2 \sum_{j=1}^L \log \gamma (\alpha_{j} ) -\log \gamma ( \sum_{j=1}^L  \alpha_{j} )

    for(i = 0; i < N; i++)
        dLogE += adPi[i] * gsl_sf_lngamma(adSumAlphaN[i]); //dLogE_2 \sum_{n=1}^N E[z_i]* lngamma( \sum_j(x_ij+alpha_j) )

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

/* A function to return the gradient for the BFGS, derivative of the 
 * expected negative log posterior wrt lambda
 * lambda are the current values, lparams contains the current parameter values (E(z_i))
 * There might be again an error with GAMMA_ITA vs (GAMMA_ITA -1), now corrected
 * GAMMA_ITA_H for the regularization seems to as expected
 */

// [[Rcpp::export]]
Rcpp::NumericVector neg_log_derive_evidence_lambda_pi(Rcpp::NumericVector ptLambda,
                                                      Rcpp::List lparams)
{
    Rcpp::IntegerMatrix aanX = as<Rcpp::IntegerMatrix>(lparams["data"]); // N x L
    Rcpp::NumericVector adPi = as<Rcpp::NumericVector>(lparams["pi"]); // N
    double GAMMA_ITA = as<double>(lparams["eta"]);
    double GAMMA_NU = as<double>(lparams["nu"]);
    double GAMMA_ITA_H = as<double>(lparams["etah"]);
    double GAMMA_NU_H = as<double>(lparams["nuh"]);
    Rcpp::List hkm = as<Rcpp::List>(lparams["hkm"]);
    Rcpp::IntegerVector hkm_index = as<Rcpp::IntegerVector>(lparams["hkm_index"]);
    
    Rcpp::List gradient = as<Rcpp::List>(lparams["gradient"]);
    Rcpp::IntegerVector gradient_index = as<Rcpp::IntegerVector>(lparams["gradient_index"]);

    int i, j;
    const int L = aanX.ncol();
    const int N = aanX.nrow();

    Rcpp::NumericVector g(L);
    Rcpp::NumericVector adDeriv(L); //derivative for each j
    Rcpp::NumericVector adStore(N); // \sum_j^L x_ij + \sum_j^L \alpha_j
    Rcpp::NumericVector adAlpha(L);
    double dSumStore = 0.0; // \sum_n^N E[z_i]* psi( \sum_j^L x_ij + \sum_j^L \alpha_j )
    double dStore = 0.0; // sum of alpha over j
    double dWeight = 0.0; // sum of z_i over i/N

    for (i = 0; i < N; i++) {  
        adStore[i] = 0.0;
        dWeight += adPi[i];
    }
    for (j = 0; j < L; j++) {
        adAlpha[j] = exp(ptLambda[j]);
        dStore += adAlpha[j];
        adDeriv[j] = dWeight* gsl_sf_psi(adAlpha[j]); //adDeriv_1[j]
        double alphaS0 = gsl_sf_psi(adAlpha[j]);
        for (i = 0; i < N; i++) {
            int dN = aanX(i, j);
            double dAlphaN = adAlpha[j] + dN;
            double psiAlphaN = dN ? gsl_sf_psi(dAlphaN) : alphaS0;
            adDeriv[j] -= adPi[i]*psiAlphaN; //adDeriv_2[j]
            adStore[i] += dAlphaN; //  \sum_j^L x_ij + \sum_j^L \alpha_j
        }
    }

    for (i = 0; i < N; i++){
        dSumStore += adPi[i] * gsl_sf_psi(adStore[i]);
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





/* calc_z calls this
 * Computes the logarithm of -p(x_i^m|z_{ik}=1, \theta) for given i, k and m, 
 * but only those terms that depend on \alpha
 * I.e. computes the negative logarithm of the unnormalized DirichletMultinomial pdf value 
 * neg_log_evidence_i(data_row, lambda_matrix(k, _),LngammaLambda0_matrix(k, _));
 * dataRow=data_row 1 x S
 * Lambda=lambda_matrix (k, _) 1 x S
 * LnGammaLambda0=LngammaLambda0_matrix(k, _) 1 x S lngamma ( \alpha_{jk}^{(m)}  )
 *
*/

// [[Rcpp::export]]
double neg_log_evidence_i(Rcpp::IntegerVector dataRow, Rcpp::NumericVector Lambda,
                          Rcpp::NumericVector LnGammaLambda0)
{
    int j;
    const int L = dataRow.length();
    double dLogE = 0.0; // \sum_j^L lngamma ( x_{ij} + alpha_j  )
    double dLogEAlpha = 0.0; // \sum_j^L lngamma ( \alpha_{jk}^{(m)}  )
    double dSumAlpha = 0.0;  // \sum_j^L \alpha_{j}
    double dSumAlphaN = 0.0; // \sum_j^L (x_{ij} + \alpha_{j})

    for (j = 0; j < L; j++) {
        const double n = dataRow[j]; // x_{ij}
        const double dAlpha = exp(Lambda[j]); //alpha_j
        const double dAlphaN = n + dAlpha; // x_{ij} + alpha_j
        // dLogEAlpha_1
        dLogEAlpha += LnGammaLambda0[j]; // \sum_j^L lngamma ( \alpha_{jk}^{(m)}  )
        dSumAlpha += dAlpha; // \sum_j^L \alpha_{j}
        dSumAlphaN += dAlphaN; // \sum_j^L (x_{ij} + \alpha_{j})
        dLogE -= n ? gsl_sf_lngamma(dAlphaN) : LnGammaLambda0[j] ;
    }
    // dLogEAlpha_2 
    dLogEAlpha -= gsl_sf_lngamma(dSumAlpha); // \sum_j^L lngamma ( \alpha_{jk}^{(m)}  ) - lngamma(  \sum_j^L (x_{ij} + \alpha_{j}) )
    dLogE += gsl_sf_lngamma(dSumAlphaN); // \sum_j^L (x_{ij} + \alpha_{j})
    return dLogE + dLogEAlpha;
}

/* Computes the E-step, i.e. the posterior probabilities of the cluster labels
 * Z=Ez K x N matrix (from the previous E-step, these are not used?)
 * data=binned.data list of M matrixes of size N x L
 * W=weights vector of length K, \pi_k^{old}
 * lambda=lambda list of M, each K times L matrix
 *
 * returns K x N matrix Z
 *
 * Why offset is subtracted?
 */

// [[Rcpp::export]]
Rcpp::NumericMatrix calc_z(Rcpp::NumericMatrix Z, Rcpp::List data,
                           Rcpp::NumericVector W, Rcpp::List Lambda)
{
    int i, j, k, m;
    Rcpp::IntegerMatrix temp = as<Rcpp::IntegerMatrix>(data[0]);
    const int N = temp.nrow();
    const int K = W.length();
    const int M = data.size();
    Rcpp::NumericMatrix evidence_matrix(M, K); // will contain the log p(\mathbf{x}_{i}^m| \theta) for sample i

    Rcpp::NumericVector L(M);
    for (i=0; i<M; i++) {
      temp = as<Rcpp::IntegerMatrix>(data[i]); //N x L matrix
      L[i] = temp.ncol(); //L
    }

    Rcpp::List LngammaLambda0(M); // lngamma ( \alpha_{jk}^{(m)}  )

    // Compute lngammaalpha_jkm
    for (m = 0; m < M; m++) {
      Rcpp::NumericMatrix Lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x L matrix
      Rcpp::NumericMatrix LngammaLambda0_matrix(K, L[m]); // lngamma ( \alpha_{jk}  )
      for(k = 0; k < K; k++){
          for(j = 0; j < L[m]; j++){
              const double dAlpha = exp(Lambda_matrix(k, j));
              LngammaLambda0_matrix(k, j) = gsl_sf_lngamma(dAlpha);
          } //j
      } //k
      LngammaLambda0[m] = LngammaLambda0_matrix;
    } //m

    for (i = 0; i < N; i ++) {
        double dSum = 0.0;
        Rcpp::NumericVector offset(M); //save the smallest negLogEvidence for each cluster(largest absolute value)
        // Compute the evidence matrix, the DirichletMultinomial pdf for given i, m and k
        for (m = 0; m < M; m++) {
          offset[m] = BIG_DBL; //1.0e9

          Rcpp::IntegerMatrix data_matrix = as<Rcpp::IntegerMatrix>(data[m]); // N x L matrix
          Rcpp::NumericMatrix lambda_matrix = as<Rcpp::NumericMatrix>(Lambda[m]); //K x L matrix
          Rcpp::NumericMatrix LngammaLambda0_matrix = as<Rcpp::NumericMatrix>(LngammaLambda0[m]); // K x L lngamma ( \alpha_{jk}^{(m)}  )
          Rcpp::IntegerVector data_row = data_matrix(i, _); // one row of X of length L
          for (k = 0; k < K; k++) {
              //Computes the logarithm of -p(\mathbf{x}_i|z_{ik}=1,\bm{\theta}) but of only those
              // terms that depend on \alpha
              // I.e. computes the negative logarithm of the unnormalized DirichletMultinomial pdf value
              //logarithm due to numerical purposes
              // computes the - log DirMulti(x_i^m | alpha_k^m), terms not depending on alpha excluded
              double dNegLogEviI =neg_log_evidence_i(data_row, lambda_matrix(k, _),
                                     LngammaLambda0_matrix(k, _));
                  //data_row x_i^m 1 x L
                  // lambda_matrix (k, _) log alpha_k^m 1 x L
                  // LngammaLambda0_matrix(k, _) 1 x L lngamma ( \alpha_{jk}^{(m)}  )


              if (dNegLogEviI < offset[m]) //1.0e9 Can this explode?
                  offset[m] = dNegLogEviI; // offset[m] is the smallest dNegLogEviI over k
              evidence_matrix(m, k) = dNegLogEviI;
          } //over K
        } //over M
        
        for (k = 0; k < K; k++) {
          Z(k, i) = 0.0;
            for (m = 0; m < M; m++) {
              Z(k, i) += (-(evidence_matrix(m, k) - offset[m])); //why offset is substracted?? For numerical reasons?
          }
          Z(k, i) = W[k] * exp(Z(k, i)); //back from logaritm, multiply by the weights i.e. \pi_k
          dSum += Z(k, i); //normalization constant
        } // k
        for (k = 0; k < K; k++)
            Z(k, i) /= dSum;
    } // i
  return Z; //K x N matrix
}

/* Function whose convergence is checked in each EM iteration
 * Computes log p(\theta,Z| X) \propto log p(X,Z|theta) + log p(theta)
 * the terms in log p(theta) that do not depend on alpha do not affect the convergence, 
 * why they are computed? In addition, one constant term is missing?
 * (\eta-1) term needs to be corrected !!??
 * W weights 1xK, these are \pi_k values
 * lambda list of length M, these are K x La matrices
 * binned.data list of length M, NxLx matrices
*/

// [[Rcpp::export]]
double neg_log_likelihood(Rcpp::NumericVector W, Rcpp::List Lambda,
                          Rcpp::List data, double eta, double nu,
                          double etah, double nuh)
{
    Rcpp::IntegerMatrix temp = as<Rcpp::IntegerMatrix>(data[0]); // N x Lx
    const int N = temp.nrow();
    const int K = W.length();
    const int M = data.size();
    int i, j, k, m;

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
              LogBAlpha(m, k) += lngammaAlpha; // \sum_j^L lng (alpha_jk)
          } //over j
          LogBAlpha(m, k) -= gsl_sf_lngamma(dSumAlphaK); // \sum_j^L lng (alpha_jk) -lngamma( \sum_j^L \alpha_jk)
      } // over k
      LngammaLambda0[m] = LngammaLambda0_matrix; // list of length M
    } // over m
    for (i = 0; i < N; i++) {
        double dProb = 0.0;
        Rcpp::NumericMatrix LogStore(M, K);
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
              dFactor += gsl_sf_lngamma(data_matrix(i, j) + 1.0); // \sum_j^L lng( x_ij +1),
          }
          // dFactor does not depend on parameters, why computed?
          dFactor -= gsl_sf_lngamma(dSum + 1.0);  // \sum_j^L lng( x_ij +1) - lngamma( \sum_j^L x_ij +1)

          for (k = 0; k < K; k++) {
              double dSumAlphaKN = 0.0; // \sum_j^L \alpha_jk + x_ij
              double dLogBAlphaN = 0.0; // lng( \alpha_jk + x_ij ) - lng( \sum_j^L \alpha_jk + x_ij )
              for (j = 0; j < L[m]; j++) {
                  int countN = data_matrix(i, j); // x_ij
                  double dAlphaN = exp(Lambda_matrix(k, j)) + countN; // \alpha_jk + x_ij
                  dSumAlphaKN += dAlphaN; // \sum_j^L \alpha_jk + x_ij
                  dLogBAlphaN += countN ? gsl_sf_lngamma(dAlphaN) : LngammaLambda0_matrix(k, j); // lng( \alpha_jk + x_ij )
              }
              dLogBAlphaN -= gsl_sf_lngamma(dSumAlphaKN); // lng( \alpha_jk + x_ij ) - lng( \sum_j^L \alpha_jk + x_ij )

              //LogStore (m,k) seems to be p(\mathbf{x}_i|\theta) but without term \Gamma(J_i^{m}+1)
              // LogStore= lng( \alpha_jk + x_ij ) - lng( \sum_j^L \alpha_jk + x_ij )
              // -\sum_j^L lng (alpha_jk) + lngamma( \sum_j^L \alpha_jk) - \sum_j^L lng( x_ij +1)
              LogStore(m, k) = dLogBAlphaN - LogBAlpha(m, k) - dFactor; // the positive? log marginal likelihood of sample x_i^m for cluster k
              if (LogStore(m, k) > offset(m))
                  offset(m) = LogStore(m, k); //offset will be the largest of LogStore(m,1:K) 
          } //over K
        } //over M
       
        //offset?
        for (k = 0; k < K; k++) {
            // normalize mixing weights
            double piK = W[k]/N; 
            dProb += piK*exp(sum(LogStore(_, k)) - sum(offset)); //sum of column k over m rows
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





/* Computes the Hessian matrix, some terms missing? Error in the order of calculation
 * lambda[[m]][k, ], 1xS
* Pi=Ez[k, ], 1x1000
* binned.data[[m]], 1000xS
* nu 1
 */

// [[Rcpp::export]]
NumericMatrix hessian(Rcpp::NumericVector Lambda, Rcpp::NumericVector Pi,
                      Rcpp::IntegerMatrix data, double nu)
{
    const int L = data.ncol();
    const int N = data.nrow();
    Rcpp::NumericMatrix Hessian(L, L);

    int i = 0, j = 0;
    Rcpp::NumericVector adAlpha(L);
    Rcpp::NumericVector adAJK(L);
    Rcpp::NumericVector adCJK(L);
    Rcpp::NumericVector adAJK0(L);
    Rcpp::NumericVector adCJK0(L);

    double dCK0 = 0.0;
    double dAK0;
    double dCSum;
    double dAlphaSum = 0.0;
    double dW = 0.0;
    double dCK = 0.0;
    double dAK;

    for (j = 0; j < L; j++) {
        adAlpha[j] = exp(Lambda[j]);
        dAlphaSum += adAlpha[j];
        adAJK0[j] = adAJK[j] = adCJK0[j] = adCJK[j] = 0.0;
        const double dPsiAlpha = gsl_sf_psi(adAlpha[j]);
        const double dPsi1Alpha = gsl_sf_psi_1(adAlpha[j]);
        for (i = 0; i < N; i++) {
            const int n = data(i, j);
            //Rcout << adCJK0[j] << " ";
            adCJK0[j] += Pi[i] * ( n ? gsl_sf_psi(adAlpha[j] + n) : dPsiAlpha ); //Pi=Ez[k,], Test that Pi[i]*n !=0 ?? or Pi[i] multiplied by the conditional
            //Rcout << adCJK0[j] << " " << Pi[i] << " "<< n << " " << Pi[i]*n << " " << gsl_sf_psi(adAlpha[j] + n) << " "<< dPsiAlpha << "\n";
            adAJK0[j] += Pi[i] * dPsiAlpha;
            adCJK[j] += Pi[i] * ( n ? gsl_sf_psi_1(adAlpha[j] + n): dPsi1Alpha );
            adAJK[j] += Pi[i] * dPsi1Alpha;
        }
    }

    for (i = 0; i < N; i++) { 
        dW += Pi[i];
        dCSum = 0.0;
        for (j = 0; j < L; j++)
            dCSum += adAlpha[j] + data(i, j);
        dCK  += Pi[i]*gsl_sf_psi_1(dCSum);
        dCK0 += Pi[i]*gsl_sf_psi(dCSum);
    }

    dAK = dW * gsl_sf_psi_1(dAlphaSum);
    dAK0 = dW * gsl_sf_psi(dAlphaSum);
    for (i = 0; i < L; i++) {
        for (j = 0; j < L; j++) {
            double dVal = 0.0;
            if (i == j) {// diagonal terms
                double dG1 = -adAlpha[i] *
                    (dAK0 - dCK0 + adCJK0[i] - adAJK0[i]);
                double dG2 = -adAlpha[i] *
                    adAlpha[i]*(dAK - dCK + adCJK[i] - adAJK[i]);
                double dG3 = adAlpha[i] * nu;
                dVal = dG1 + dG2 + dG3;
            } else
                dVal = -adAlpha[i] * adAlpha[j] * (dAK - dCK);
            Hessian(i, j) = dVal;
        }
    }
    return Hessian;
}


