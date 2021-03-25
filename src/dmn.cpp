#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>

#include <Rcpp.h>
using namespace Rcpp;

#define BIG_DBL 1.0e9

// computes the regularization term using lambda, alpha= exp(lambda)

double REG_H_LAMBDA(NumericVector lambda) {
  int i;
  double hk = 0.0;

  //from second element to the last
  for (i = 1; i < lambda.length(); i++) {
    double diff = exp(lambda[i]) - exp(lambda[i-1]);
    hk += diff * diff;
  }
  return hk;
}

// [[Rcpp::export]]
void disable_gsl_error_handler() {
 gsl_set_error_handler_off();
}

//computes the value of function g for j=index, input lambda values
double REG_DERIV_LAMBDA_G(NumericVector lambda, int index) {
  int prev = index - 1, next = index + 1;
  if (index == 0) prev = 0;
  if (index == lambda.length()-1) next = index;
  return 2*(2*exp(lambda[index]) - exp(lambda[prev]) - exp(lambda[next]));
}

/*Soft K-means clustering treats the cluster assignments as probability distributions
over the clusters. There is a connection between Euclidean distance and multivariate normal
models with a fixed covariance, soft K-means can be expressed as a multivariate
normal mixture model*/

// [[Rcpp::export]]
List soft_kmeans(NumericMatrix data,
                 int K,
                 bool verbose=false,
                 bool randomInit=true,
                 NumericMatrix centers=NumericMatrix(), //if not given, initialized randomly
                 double stiffness=50.0,
                 double threshold=1.0e-6,
                 int maxIt=1000,
                 bool rowNorm=true) {
    /*data=kmeans.binned.data, N x (M*S)
     centers=kmeanspp.centers, K x (M*S), can
     rowNorm=F */
    const int S = data.ncol(), N = data.nrow();
    RNGScope scope;

    int i, j, k, iter = 0;

    NumericVector W(K);
    NumericMatrix Z(K, N); //cluster membership probabilities
    List data_dimnames = data.attr("dimnames");
    Z.attr("dimnames") = List::create(CharacterVector(), data_dimnames[0]); //sample names
    NumericMatrix Mu(K, S);
    Mu.attr("dimnames") = List::create(CharacterVector(), data_dimnames[1]); //bin names

    double dMaxChange = BIG_DBL; //1.0e9

    if (verbose)
        Rprintf("  Soft kmeans\n");

    NumericMatrix Y(N, S);
    NumericVector Mu_row(S);

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

    Function sample("sample"); //access to R's sample function

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
      Mu.attr("dimnames") = List::create(CharacterVector(), data_dimnames[1]);

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

    return List::create(_["centers"] = Mu,
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
double neg_log_evidence_lambda_pi(NumericVector lambda, List lparams)
{
    int i, j;

    IntegerMatrix aanX = as<IntegerMatrix>(lparams["data"]); // N x L_x
    NumericVector adPi = as<NumericVector>(lparams["pi"]);   // size N, these are z_{km}
    double GAMMA_ITA = as<double>(lparams["eta"]);
    double GAMMA_NU = as<double>(lparams["nu"]);
    double GAMMA_ITA_H = as<double>(lparams["etah"]);
    double GAMMA_NU_H = as<double>(lparams["nuh"]);

    const int S = aanX.ncol(); // L_x
    const int N = aanX.nrow();
    
    /*dLogE collects the terms \sum_{n=1}^N \sum_{j=1}^S E[z_i] \log \gamma (x_ij+ alpha_j)
    and \sum_{n=1}^N E[z_i]* lng( \sum_j(x_ij+alpha-j) )*/
    double dLogE = 0.0;
    
    double dLogEAlpha = 0.0; // \log \gamma ( \sum_{j=1}^S  \alpha_{j} )
    double dSumAlpha = 0.0; // sum of alpha \sum_{j=1}^S \alpha_{j}
    double dSumLambda = 0.0; // sum of lambda?
    double dHk = 0.0;
    double dWeight = 0.0; // \sum_{n=1}^N E(z_i)

    NumericVector adSumAlphaN(N); // \sum_{j}^S \alpha_{j}+x_{ij}

    for (i = 0; i < N; i++) {
        adSumAlphaN[i] = 0.0;
        dWeight += adPi[i]; //sum over z_i
    }

    for (j = 0; j < S; j++) {
        const double dLambda = lambda[j];
        const double dAlpha = exp(dLambda);
        /* compute the logarithm of the Gamma function,
        dAlpha can not be zero or negative.
        Function computed using the real Lanczos method */
        dLogEAlpha += gsl_sf_lngamma(dAlpha); // \sum_{j=1}^S \log \gamma (\alpha_{j} )
        dSumLambda += dLambda;
        dSumAlpha += dAlpha;
        const double lngammaAlpha0 = gsl_sf_lngamma(dAlpha); //lngamma of \alpha_{j}
        for (i = 0; i < N; i++) {
            const double dN = aanX(i, j); // x_{ij}
            const double dAlphaN = dAlpha + dN; // \alpha_{j}+x_{ij}
            const double lngammaAlphaN = dN ? gsl_sf_lngamma(dAlphaN) : lngammaAlpha0; //if dN exists or is non-zero, compute lnGamma(dAlphaN), else compute lnGamma(DaLpha)
            adSumAlphaN[i] += dAlphaN; // \sum_{j}^S \alpha_{j}+x_{ij} , weight by pi
            dLogE -= adPi[i] * lngammaAlphaN; //weight by pi, -\sum_i E[z_i] *\sum_j lngamma(x_{ij}+alpha_j)
        }
    }
    dLogEAlpha -= gsl_sf_lngamma(dSumAlpha);// \sum_{j=1}^S \log \gamma (\alpha_{j} ) -\log \gamma ( \sum_{j=1}^S  \alpha_{j} )

    for(i = 0; i < N; i++)
        dLogE += adPi[i] * gsl_sf_lngamma(adSumAlphaN[i]); // \sum_{n=1}^N E[z_i]* lngamma( \sum_j(x_ij+alpha_j) )

    double reg_term;

    if ((GAMMA_ITA_H==0) && (GAMMA_NU_H==0)) {
      reg_term = 0.0;
    } else {
      double hk = REG_H_LAMBDA(lambda); // computes the value of the regularization term
      reg_term = GAMMA_NU_H*hk - (GAMMA_ITA_H-1)*log(hk); //This was nuh *hg-etah*log (hk), it is now corrected to (GAMMA_ITA_H-1)!!!
    }



    //should it be (GAMMA_ITA-1)*dSumLambda???
    return dLogE + dWeight*dLogEAlpha + // complete data likelihood term
      GAMMA_NU*dSumAlpha - (GAMMA_ITA - 1) * dSumLambda + reg_term; //prior term, ITA corrected to ITA-1 !!!
}

/* A function to return the gradient for the BFGS, derivative of the 
 * expected negative log posterior wrt lambda
 * lambda are the current values, lparams contains the current parameter values (E(z_i))
 * There might be again an error with GAMMA_ITA vs (GAMMA_ITA -1), now corrected
 * GAMMA_ITA_H for the regularization seems to as expected
 */

// [[Rcpp::export]]
NumericVector neg_log_derive_evidence_lambda_pi(NumericVector ptLambda,
        List lparams)
{
    IntegerMatrix aanX = as<IntegerMatrix>(lparams["data"]); // N x S
    NumericVector adPi = as<NumericVector>(lparams["pi"]); // N
    double GAMMA_ITA = as<double>(lparams["eta"]);
    double GAMMA_NU = as<double>(lparams["nu"]);
    double GAMMA_ITA_H = as<double>(lparams["etah"]);
    double GAMMA_NU_H = as<double>(lparams["nuh"]);
    List hkm = as<List>(lparams["hkm"]);
    IntegerVector hkm_index = as<IntegerVector>(lparams["hkm_index"]);

    int i, j;
    const int S = aanX.ncol(), N = aanX.nrow();

    NumericVector g(S);
    NumericVector adDeriv(S); //derivative for each j
    NumericVector adStore(N); // \sum_j^S x_ij + \sum_j^S \alpha_j
    NumericVector adAlpha(S);
    double dSumStore = 0.0; // \sum_n^N E[z_i]* psi( \sum_j^S x_ij + \sum_j^S \alpha_j )
    double dStore = 0.0; // sum of alpha over j
    double dWeight = 0; // sum of z_i over i/N

    for (i = 0; i < N; i++) {  
        adStore[i] = 0.0;
        dWeight += adPi[i];
    }

    for (j = 0; j < S; j++) {
        adAlpha[j] = exp(ptLambda[j]);
        dStore += adAlpha[j];
        adDeriv[j] = dWeight* gsl_sf_psi(adAlpha[j]);
        double alphaS0 = gsl_sf_psi(adAlpha[j]);
        for (i = 0; i < N; i++) {
            int dN = aanX(i, j);
            double dAlphaN = adAlpha[j] + dN;

            double psiAlphaN = dN ? gsl_sf_psi(dAlphaN) : alphaS0;
            adDeriv[j] -= adPi[i]*psiAlphaN;
            //            adDeriv[j] -= adPi[i]*gsl_sf_psi (dAlphaN);
            adStore[i] += dAlphaN; //  \sum_j^S x_ij + \sum_j^S \alpha_j
        }
    }

    for (i = 0; i < N; i++)
        dSumStore += adPi[i] * gsl_sf_psi(adStore[i]);
    dStore = dWeight * gsl_sf_psi(dStore);

    double hk = REG_H_LAMBDA(ptLambda);

    if (hkm_index[0] < hkm.size()) hkm[hkm_index[0]] = hk;
    else Rprintf("hkm_index exceeds hkm vector length!!! hkm_index: %i, hkm.length: %i\n",
                 hkm_index[0], hkm.size());

    hkm_index[0] += 1;

    for (j = 0; j < S; j++) {

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
  return g;
}
/* Computes the logarithm of -p(x_i^m|z_{ik}=1, \theta) for given i, k and m, 
 * but only those terms that depend on \alpha
 * I.e. computes the negative logarithm of the unnormalized DirichletMultinomial pdf value 
 * neg_log_evidence_i(data_row, lambda_matrix(k, _),LngammaLambda0_matrix(k, _));
 * dataRow=data_row 1 x S
 * Lambda=lambda_matrix (k, _) 1 x S
 * LnGammaLambda0=LngammaLambda0_matrix(k, _) 1 x S lngamma ( \alpha_{jk}^{(m)}  )
 *
*/

// [[Rcpp::export]]
double neg_log_evidence_i(IntegerVector dataRow, NumericVector Lambda,
                          NumericVector LnGammaLambda0)
{
    int j;
    const int S = dataRow.length();
    double dLogE = 0.0; // \sum_j^S lngamma ( x_{ij} + alpha_j  )
    double dLogEAlpha = 0.0; // \sum_j^S lngamma ( \alpha_{jk}^{(m)}  )
    double dSumAlpha = 0.0;  // \sum_j^S \alpha_{j}
    double dSumAlphaN = 0.0; // \sum_j^S (x_{ij} + \alpha_{j})

    for (j = 0; j < S; j++) {
        const double n = dataRow[j]; // x_{ij}
        const double dAlpha = exp(Lambda[j]); //alpha_j
        const double dAlphaN = n + dAlpha; // x_{ij} + alpha_j

        dLogEAlpha += LnGammaLambda0[j]; // \sum_j^S lngamma ( \alpha_{jk}^{(m)}  )
        dSumAlpha += dAlpha; // \sum_j^S \alpha_{j}
        dSumAlphaN += dAlphaN; // \sum_j^S (x_{ij} + \alpha_{j})
        dLogE -= n ? gsl_sf_lngamma(dAlphaN) : LnGammaLambda0[j] ;
    }

    dLogEAlpha -= gsl_sf_lngamma(dSumAlpha); // \sum_j^S lngamma ( \alpha_{jk}^{(m)}  ) - lngamma(  \sum_j^S (x_{ij} + \alpha_{j}) )
    dLogE += gsl_sf_lngamma(dSumAlphaN); // \sum_j^S (x_{ij} + \alpha_{j})

    return dLogE + dLogEAlpha;
}

/* Computes the E-step, i.e. the posterior probabilities of the cluster labels
 * Z=Ez K x N matrix (from the previous E-step, these are not used?)
 * data=binned.data list of M matrixes of size N x S
 * W=weights vector of length K, \pi_k^{old}
 * lambda=lambda list of M, each K times S matrix
 *
 * returns K x N matrix Z
 *
 * Why offset is subtracted?
 */

// [[Rcpp::export]]
NumericMatrix calc_z(NumericMatrix Z, List data,
                     NumericVector W, List Lambda)
{
    int i, j, k, m;
    IntegerMatrix temp = as<IntegerMatrix>(data[0]);
    const int N = temp.nrow();
    const int K = W.length();
    const int M = data.size();
    NumericMatrix evidence_matrix(M, K); // will contain the log p(\mathbf{x}_{i}^m| \theta) for sample i

    NumericVector S(M);
    for (i=0; i<M; i++) {
      temp = as<IntegerMatrix>(data[i]); //N x S matrix
      S[i] = temp.ncol(); //S
    }

    List LngammaLambda0(M); // lngamma ( \alpha_{jk}^{(m)}  )

    for (m = 0; m < M; m++) {
      NumericMatrix Lambda_matrix = as<NumericMatrix>(Lambda[m]); //K x S matrix
      NumericMatrix LngammaLambda0_matrix(K, S[m]); // lngamma ( \alpha_{jk}  )
      for(k = 0; k < K; k++){
          for(j = 0; j < S[m]; j++){
              const double dAlpha = exp(Lambda_matrix(k, j));
              LngammaLambda0_matrix(k, j) = gsl_sf_lngamma(dAlpha);
          }
      }
      LngammaLambda0[m] = clone(LngammaLambda0_matrix);
    }

    for (i = 0; i < N; i ++) {
        double dSum = 0.0;
        NumericVector offset(M); //save the smallest negLogEvidence for each cluster(largest absolute value)
        // Compute the evidence matrix, the DirichletMultinomial pdf for given i, m and k
        for (m = 0; m < M; m++) {
          offset[m] = BIG_DBL; //1.0e9

          IntegerMatrix data_matrix = as<IntegerMatrix>(data[m]); // N x S matrix
          NumericMatrix lambda_matrix = as<NumericMatrix>(Lambda[m]); //K x S matrix
          NumericMatrix LngammaLambda0_matrix = as<NumericMatrix>(LngammaLambda0[m]); // K x S lngamma ( \alpha_{jk}^{(m)}  )
          IntegerVector data_row = data_matrix(i, _); // one row of X of length S
          for (k = 0; k < K; k++) {
              //Computes the logarithm of -p(\mathbf{x}_i|z_{ik}=1,\bm{\theta}) but of only those
              // terms that depend on \alpha
              // I.e. computes the negative logarithm of the unnormalized DirichletMultinomial pdf value
              //logarithm due to numerical purposes
              double dNegLogEviI =neg_log_evidence_i(data_row, lambda_matrix(k, _),
                                     LngammaLambda0_matrix(k, _));
                  //data_row 1 x S
                  // lambda_matrix (k, _) 1 x S
                  // LngammaLambda0_matrix(k, _) 1 x S lngamma ( \alpha_{jk}^{(m)}  )


              if (dNegLogEviI < offset[m]) //1.0e9 Can this explode?
                  offset[m] = dNegLogEviI;
              evidence_matrix(m, k) = dNegLogEviI;
          } //over K
        } //over M
        
        for (k = 0; k < K; k++) {
          Z(k, i) = 0;
            for (m = 0; m < M; m++) {
              Z(k, i) += (-(evidence_matrix(m, k) - offset[m])); //why offset is substracted?? For numerical reasons?
          }
          Z(k, i) = W[k] * exp(Z(k, i)); //back from logaritm, multiply by the weights i.e. \pi_k
          dSum += Z(k, i); //normalization constant
        }
        for (k = 0; k < K; k++)
            Z(k, i) /= dSum;
    }
  return Z; //K x N matrix
}

/* Function whose convergence is checked in each EM iteration
 * Computes log p(\theta,Z| X) \propto log p(X,Z|theta) + log p(theta)
 * the terms in log p(theta) that do not depend on alpha do not affect the convergence, 
 * why they are computed? In addition, one constant term is missing?
 * (\eta-1) term needs to be corrected !!??
 * W weights 1xK, these are \pi_k values
 * lambda list of KxS matrices
 * binned.data NxS
*/

// [[Rcpp::export]]
double neg_log_likelihood(NumericVector W, List Lambda,
                          List data, double eta, double nu,
                          double etah, double nuh)
{
    IntegerMatrix temp = as<IntegerMatrix>(data[0]);
    const int N = temp.nrow();
    const int K = W.length();
    const int M = data.size();
    int i, j, k, m;

    NumericVector S(M);
    for (m = 0; m < M; m++) {
      temp = as<IntegerMatrix>(data[m]);
      S[m] = temp.ncol();
    }
    List LngammaLambda0(M);

    // return -dRet - dL5 - dL6 - dL7 - dL8 - regterm1 - regterm2;
    /* -log \sum_k^K(\pi_k*p(x|\theta_k))+M*S*K*lng(eta)- eta*K*M*S*log(nu)
    * + nu* \sum_j^S \alpha_kj^m -eta * \sum_j^S \lambda_kj^m
    -M * K * (etah * log(nuh) - gsl_sf_lngamma(etah)) - \sum_k^K ((etah - 1) * log(hkm) - nuh*hkm) */
    double dRet = 0.0; //dRet of log \sum_k^K(\pi_k*p(x|\theta_k)) without Gamma (J_i^m-1) term
    double dL5 = 0.0; // -M*S*K*lng(eta) from prior
    double dL6 = 0.0; // eta*K*M*S*log(nu) from prior
    double dL7 = 0.0; //-nu* \sum_j^S \alpha_kj^m from prior
    double dL8 = 0.0; // eta * \sum_j^S \lambda_kj^m from prior
    double regterm1 = 0.0; //M * K * (etah * log(nuh) - gsl_sf_lngamma(etah)); from prior
    double regterm2 = 0.0; // \sum_k^K ((etah - 1) * log(hkm) - nuh*hkm);  form prior

    NumericMatrix LogBAlpha(M, K); // // \sum_j^S lng (alpha_jk) -lngamma( \sum_j^S \alpha_jk)

    for (m = 0; m < M; m++) {
      NumericMatrix LngammaLambda0_matrix(K, S[m]); // lng (alpha_jk)
      NumericMatrix Lambda_matrix = as<NumericMatrix>(Lambda[m]);

      for (k = 0; k < K; k++){
          double dSumAlphaK = 0.0;
          LogBAlpha(m, k) = 0.0;
          for (j = 0; j < S[m]; j++){
              double dAlpha = exp(Lambda_matrix(k, j)); // \ alpha_jk
              double lngammaAlpha = gsl_sf_lngamma(dAlpha); // lng (alpha_jk)
              LngammaLambda0_matrix(k, j) = lngammaAlpha; // // lng (alpha_jk)

              dSumAlphaK += dAlpha; // \sum_j^S \alpha_jk
              LogBAlpha(m, k) += lngammaAlpha; // \sum_j^S lng (alpha_jk)
          } //over j
          LogBAlpha(m, k) -= gsl_sf_lngamma(dSumAlphaK); // \sum_j^S lng (alpha_jk) -lngamma( \sum_j^S \alpha_jk)
      } // over k
      LngammaLambda0[m] = clone(LngammaLambda0_matrix);
    } // over m
    for (i = 0; i < N; i++) {
        double dProb = 0.0;
        NumericMatrix LogStore(M, K);
        NumericVector offset(M);
        // Computes the log p(x_i^m|, z_ik=1, theta) for k and m, results in matrix LogStore
        for (m = 0; m < M; m++) {
          double dSum = 0.0; // \sum_j^S x_ij
          double dFactor = 0.0; // \sum_j^S lng( x_ij +1) - lngamma( \sum_j^S x_ij +1)
          offset(m) = -BIG_DBL;

          IntegerMatrix data_matrix = as<IntegerMatrix>(data[m]);
          NumericMatrix Lambda_matrix = as<NumericMatrix>(Lambda[m]);
          NumericMatrix LngammaLambda0_matrix = as<NumericMatrix>(LngammaLambda0[m]); // lng (alpha_jk)

          for (j = 0; j < S[m]; j++) {
              dSum += data_matrix(i, j); // dSum_m
              dFactor += gsl_sf_lngamma(data_matrix(i, j) + 1.0); // \sum_j^S lng( x_ij +1)
          }
          dFactor -= gsl_sf_lngamma(dSum + 1.0);  // \sum_j^S lng( x_ij +1) - lngamma( \sum_j^S x_ij +1)

          for (k = 0; k < K; k++) {
              double dSumAlphaKN = 0.0; // \sum_j^S \alpha_jk + x_ij
              double dLogBAlphaN = 0.0; // lng( \alpha_jk + x_ij ) - lng( \sum_j^S \alpha_jk + x_ij )
              for (j = 0; j < S[m]; j++) {
                  int countN = data_matrix(i, j); // x_ij
                  double dAlphaN = exp(Lambda_matrix(k, j)) + countN; // \alpha_jk + x_ij
                  dSumAlphaKN += dAlphaN; // \sum_j^S \alpha_jk + x_ij
                  dLogBAlphaN += countN ? gsl_sf_lngamma(dAlphaN) : LngammaLambda0_matrix(k, j); // lng( \alpha_jk + x_ij )
              }
              dLogBAlphaN -= gsl_sf_lngamma(dSumAlphaKN); // lng( \alpha_jk + x_ij ) - lng( \sum_j^S \alpha_jk + x_ij )

              //LogStore (m,k) seems to be p(\mathbf{x}_i|\theta) but without term \Gamma(J_i^{m}+1)
              // LogStore= lng( \alpha_jk + x_ij ) - lng( \sum_j^S \alpha_jk + x_ij )
              // -\sum_j^S lng (alpha_jk) + lngamma( \sum_j^S \alpha_jk) - \sum_j^S lng( x_ij +1)
              LogStore(m, k) = dLogBAlphaN - LogBAlpha(m, k) - dFactor; // the positive? log marginal likelihood of sample x_i^m for cluster k
              if (LogStore(m, k) > offset(m))
                  offset(m) = LogStore(m, k);
          } //over K
        } //over M
       
        //offset?
        for (k = 0; k < K; k++) {
            double piK = W[k]/N; //\pi weights not E[z]???
            dProb += piK*exp(sum(LogStore(_, k)) - sum(offset)); //sum of column k over m rows
        }
        dRet += log(dProb)+sum(offset); //logarithm of the normalization term, sum over all n
    } //over N

    dL5 = -sum(S) * K * gsl_sf_lngamma(eta); // -M*S*K*lng(eta)
    dL6 = eta * K * sum(S) * log(nu);  //eta*K*M*S*log(nu)

    if ((etah!=0) || (nuh!=0)) {
      regterm1 = M * K * (etah * log(nuh) - gsl_sf_lngamma(etah));
    }

    for (m = 0; m < M; m++) {
      NumericMatrix Lambda_matrix = as<NumericMatrix>(Lambda[m]); //K x S matrix

      for (i = 0; i < K; i++) {

        if ((etah!=0) || (nuh!=0)) {
          double hkm = sum(diff(exp(Lambda_matrix(i, _))) * diff(exp(Lambda_matrix(i, _))));
          regterm2 += (etah - 1) * log(hkm) - nuh*hkm; // \sum_k^K (etah - 1) * log(hkm) - nuh*hkm;
        }

        dL7 += sum(exp(Lambda_matrix(i, _))); // \sum_j^S \alpha_kj^m
        dL8 += sum(Lambda_matrix(i, _));  //  \sum_j^S \lambda_kj^m
      } // over K
    }  // over M
    dL7 *= -nu; //-nu* \sum_j^S \alpha_kj^m
    dL8 *= (eta-1); // eta * \sum_j^S \lambda_kj^m SHOULD BE (eta-1)
    return -dRet - dL5 - dL6 - dL7 - dL8 - regterm1 - regterm2;
    //

}

/* Used when optimizing the shift and flip. What does this try to optimize???
 * shift_dist: the amount of shift (-120) dist.candidates(1)
 * index of the sample
 * do we flip (F)
 * data =binned.data,
 * alpha=lapply(lambda, exp),
 * Z=Ez
 */

// [[Rcpp::export]]
double optimization_func(IntegerVector shift_dist,
                         IntegerVector index,
                         LogicalVector flip,
                         List data,
                         List alpha,
                         NumericMatrix Z) {

  double res = 0;
  int K = Z.nrow(); //cluster number
  Function shift_and_flip_signal("shift.and.flip.signal"); //this function is in util.R

  if (shift_dist.length() != index.length()) {
    throw std::length_error("Shift distance vector and index length do not match.");
  }

  if (flip.length() != index.length()) {
    throw std::length_error("Flip vector and index length not match.");
  }

//  Rprintf("shift_optimization: shift_dist:%d\n", shift_dist);

  //shift.and.flip.signal <- function(data=binned.data, indices=wrap(index), dist, flip)
  //shifts and flips the samples of index index by shift_dist, can be also used to flip the
  // samples, the shifted and flipped samples should be the same
  // This function is in util.R
  List shifteddata = shift_and_flip_signal(data,
                                           wrap(index), //onverting from C++ to R using Rcpp::wrap(obj)
                                           wrap(shift_dist),
                                           wrap(flip));

  //iterate over data types
  for (int m = 0; m < shifteddata.length(); m++) {
    NumericMatrix datamatrix = as<NumericMatrix>(shifteddata[m]);
    NumericMatrix alphamatrix = as<NumericMatrix>(alpha[m]);

    // iterate over rows/samples
    for (int i = 0; i < index.length(); i++) {
      int index_zerob = index[i] - 1;
      NumericVector Xi = datamatrix.row(index_zerob);

      // iterate over clusters
      for (int k = 0; k < K; k++) {
        NumericVector alphak = alphamatrix(k, _);

        double sumlgammaXialpha = 0.0; // \sum_j^S lng (x_ij + alpha_kj)
        double sumXialpha = 0.0;  // \sum_j^S (x_ij+\alpha_kj)
        double Xisum = 0.0; // \sum_j^S x_ij
        double sumlfactXi = 0.0; // \sum_j lngamma( x_ij +1 )

        // iterate over columns
        for(int s = 0; s < datamatrix.ncol(); s++) {
          double tmpsum = Xi[s] + alphak[s];
          sumlgammaXialpha += gsl_sf_lngamma(tmpsum); // \sum_j^S lng (x_ij + alpha_kj)
          sumXialpha += tmpsum; // \sum_j^S (x_ij+\alpha_kj)
          Xisum += Xi[s]; // \sum_j^S x_ij
          sumlfactXi += gsl_sf_lnfact(Xi[s]); // \sum_j lngamma( x_ij +1 )
        }

        double tmp = sumlgammaXialpha - // \sum_j^S lng (x_ij + alpha_kj)
                     gsl_sf_lngamma(sumXialpha) + // lngamma ( \sum_j^S (x_ij+\alpha_kj) )
                     gsl_sf_lnfact(Xisum) -      //lngamma ( \sum_j^S x_ij +1)
                     sumlfactXi;                 // \sum_j lngamma( x_ij +1 )
        res += Z(k, index_zerob)*tmp;
      }
    }
  }
//  Rprintf("shift_optimization: Result: %f\n", -res);
  return -res;
}



/* Computes the Hessian matrix, some terms missing? Error in the order of calculation
 * lambda[[m]][k, ], 1xS
* Pi=Ez[k, ], 1x1000
* binned.data[[m]], 1000xS
* nu 1
 */

// [[Rcpp::export]]
NumericMatrix hessian(NumericVector Lambda, NumericVector Pi,
                      IntegerMatrix data, double nu)
{
    const int S = data.ncol(), N = data.nrow();
    NumericMatrix Hessian(S, S);

    int i = 0, j = 0;
    NumericVector adAlpha(S);
    NumericVector adAJK(S);
    NumericVector adCJK(S);
    NumericVector adAJK0(S);
    NumericVector adCJK0(S);

    double dCK0 = 0.0, dAK0;
    double dCSum, dAlphaSum = 0.0, dW = 0.0, dCK = 0.0, dAK;

    for (j = 0; j < S; j++) {
        adAlpha[j] = exp(Lambda[j]);
        dAlphaSum += adAlpha[j];
        adAJK0[j] = adAJK[j] = adCJK0[j] = adCJK[j] = 0.0;
        const double dPsiAlpha = gsl_sf_psi(adAlpha[j]);
        const double dPsi1Alpha = gsl_sf_psi_1(adAlpha[j]);
        for (i = 0; i < N; i++) {
            const int n = data(i, j);
            //Rcout << adCJK0[j] << " ";
            adCJK0[j] += Pi[i] * n ? gsl_sf_psi(adAlpha[j] + n) : dPsiAlpha; //Pi=Ez[k,], Test that Pi[i]*n !=0 ?? or Pi[i] multiplied by the conditional
            //Rcout << adCJK0[j] << " " << Pi[i] << " "<< n << " " << Pi[i]*n << " " << gsl_sf_psi(adAlpha[j] + n) << " "<< dPsiAlpha << "\n";
            adAJK0[j] += Pi[i] * dPsiAlpha;
            adCJK[j] += Pi[i] * n ? gsl_sf_psi_1(adAlpha[j] + n): dPsi1Alpha;
            adAJK[j] += Pi[i] * dPsi1Alpha;
        }
    }

    for (i = 0; i < N; i++) { 
        dW += Pi[i];
        dCSum = 0.0;
        for (j = 0; j < S; j++)
            dCSum += adAlpha[j] + data(i, j);
        dCK  += Pi[i]*gsl_sf_psi_1(dCSum);
        dCK0 += Pi[i]*gsl_sf_psi(dCSum);
    }

    dAK = dW * gsl_sf_psi_1(dAlphaSum);
    dAK0 = dW * gsl_sf_psi(dAlphaSum);
    for (i = 0; i < S; i++) {
        for (j = 0; j < S; j++) {
            double dVal = 0.0;
            if (i == j) {
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
/*
static void group_output(struct data_t *data, double** aadZ)
{
    const int N = data->N, K = data->K;
    int i, k;
    for(k = 0; k < K; k++)
        for (i = 0; i < N; i++)
            data->group[k * N + i] = aadZ[k][i];
}
*/

// [[Rcpp::export]]
List mixture_output(IntegerMatrix data, NumericVector W,
                    NumericMatrix Lambda, NumericMatrix Err)
{
    const int N = data.nrow(), S = data.ncol(), K = W.length();
    int i, k;

    NumericVector mixture_wt(K);
    NumericMatrix fit_lower(K, S), fit_upper(K, S), fit_mpe(K, S);

    fit_lower.attr("dimnames") = Lambda.attr("dimnames");
    fit_upper.attr("dimnames") = Lambda.attr("dimnames");
    fit_mpe.attr("dimnames") = Lambda.attr("dimnames");

    for (k = 0; k < K; k++)
        mixture_wt[k] = W[k] / N;

    for (i = 0; i < S; i++) {
        for (k = 0; k < K; k++) {
            double dErr = Err(k, i), dL = 0.0, dU = 0.0;
            int bIll = FALSE;
            if (dErr >= 0.0) {
                dErr = sqrt(dErr);
                if (dErr < 100.0) {
                    dL =  exp(Lambda(k, i) - 2.0*dErr);
                    dU =  exp(Lambda(k, i) + 2.0*dErr);
                } else bIll = TRUE;
            } else bIll = TRUE;

            if (bIll)
                dL = dU = R_NaN;
            fit_lower(k, i) = dL;
            fit_mpe(k, i) = exp(Lambda(k, i));
            fit_upper(k, i) = dU;
        }
    }
    return List::create(_["Lower"]=fit_lower,
                        _["Upper"]=fit_upper,
                        _["Estimate"]=fit_mpe,
                        _["Mixture"]=mixture_wt);
}

/*
// [[Rcpp::export]]
List dirichlet_fit(IntegerMatrix data, int K, int seed = -1,
                        int maxIt=250, bool verbose=true,
                        double eta=0.1, double nu=0.1, double stiffness=50.0,
                        bool randomInit=true)
{
    const int N = data.nrow(), S = data.ncol();
    int i, j, k;
    RNGScope rng; //initialize RNG

    if(seed != -1) {
        Function setseed("set.seed");
        setseed(seed);
    }

    NumericMatrix Z(K, N), Lambda(K, S), Err(K, S);
    NumericVector W(K);

    // soft k means initialiser
    List kmeans_result = soft_kmeans(data, K, verbose,
                                     randomInit, stiffness);
    Lambda = as<NumericMatrix>(kmeans_result["centers"]);
    //W = kmeans_result["weights"];
    Z = as<NumericMatrix>(kmeans_result["labels"]);

    for (k = 0; k < K; k++) {
        for (i = 0; i < N; i++)
            W[k] += Z(k, i);
    }

    if (verbose)
        Rprintf("  Expectation Maximization setup\n");

    for (k = 0; k < K; k++) {
        for (j = 0; j < S; j++) {
            const double x = Lambda(k, j);
            Lambda(k, j) = (x > 0.0) ? log(x) : -10;
        }
        Lambda(k, _) = optimise_lambda_k(Lambda(k, _), data, Z(k, _), eta, nu);
    }

    // simple EM algorithm
    int iter = 0;
    double dNLL = 0.0, dNew = 0.0, dChange = BIG_DBL;

    if (verbose)
        Rprintf("  Expectation Maximization\n");

    while (dChange > 1.0e-6 && iter < maxIt) {
        if (verbose)
            Rprintf("  Calculating Ez...\n");

        Z = calc_z(Z, data, W, Lambda); // latent var expectation

        if (verbose)
            Rprintf("  Optimizing wrt Lambda values...\n");

        for (k = 0; k < K; k++) // mixture components, given pi
            Lambda(k, _) = optimise_lambda_k(Lambda(k, _), data, Z(k, _), eta, nu);

        for (k = 0; k < K; k++) { // current likelihood & weights
            W[k] = 0.0;
            for(i = 0; i < N; i++)
                W[k] += Z(k, i);
        }

        if (verbose)
            Rprintf("  Calculate negative loglikelihood...\n");

        dNew = neg_log_likelihood(W, Lambda, data, eta, nu);
        dChange = fabs(dNLL - dNew);
        dNLL = dNew;
        iter++;
        checkUserInterrupt();
        if (verbose && (iter % 1) == 0)
            Rprintf("    iteration %d change %.6f\n", iter, dChange);
    }

    // hessian
    if (verbose)
        Rprintf("  Hessian\n");

    gsl_permutation *p = gsl_permutation_alloc(S);
    double dLogDet = 0.0, dTemp;
    int signum, status;

    for (k = 0; k < K; k++) {
        if (k > 0)
            dLogDet += 2.0 * log(N) - log(W[k]);

        RcppGSL::matrix<double> Hessian(hessian(Lambda(k, _), Z(k, _), data, nu));
        RcppGSL::matrix<double> InverseHessian(S, S);

        status = gsl_linalg_LU_decomp(Hessian, p, &signum);
        gsl_linalg_LU_invert(Hessian, p, InverseHessian);
        for (j = 0; j < S; j++) {
            Err(k, j) = InverseHessian(j, j);
            dTemp = Hessian(j, j);
            dLogDet += log(fabs(dTemp));
        }
        Hessian.free();
        InverseHessian.free();
    }
    gsl_permutation_free(p);

    // results
    List result;
    double dP = K * S + K - 1;
    double laplace = dNLL + 0.5 * dLogDet - 0.5 * dP * log(2. * M_PI);
    double bic = dNLL + 0.5 * log(N) * dP;
    double aic = dNLL + dP;

    result["GoodnessOfFit"] = NumericVector::create(_["NLE"]=dNLL,
                                                    _["LogDet"]=dLogDet,
                                                    _["Laplace"]=laplace,
                                                    _["BIC"]=bic,
                                                    _["AIC"]=aic);
    //group and fit results must be transposed
    result["Group"] = Z;

    List mix_list = mixture_output(data, W, Lambda, Err);
    result["Mixture"] = List::create(_["Weight"] = mix_list["Mixture"]);

    result["Fit"] = List::create(_["Estimate"]=mix_list["Estimate"],
                                 _["Upper"]=mix_list["Upper"],
                                 _["Lower"]=mix_list["Lower"]);

    return result;
}
*/