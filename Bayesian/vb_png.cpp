# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>
#include <time.h>

#include "depend.hpp"
#include "plinkfun.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export()]]
List emMultiple (arma::mat X,
                 arma::mat Y,
                 Rcpp::Nullable<Rcpp::NumericMatrix> mu0     = R_NilValue,
                 Rcpp::Nullable<Rcpp::NumericVector> sigb0   = R_NilValue,
                 Rcpp::Nullable<Rcpp::NumericMatrix> Theta0  = R_NilValue,
                 Rcpp::Nullable<Rcpp::NumericVector> Lambda0 = R_NilValue,
                 Rcpp::Nullable<Rcpp::NumericMatrix> Alpha0  = R_NilValue,
                 bool   fixlambda     = true, //if Lambda0 is not NULL, fixlambda should be false
                 int    maxit         = 10^3,
                 double epsStopLogLik = 10^-5
) {
  // inputs
  const int N = Y.n_rows ;
  const int K = Y.n_cols ;
  const int p = X.n_cols ;

  int i=0, j=0, k=0;

  //centeralize X and Y
  mat xmean = mean(X);
  for (i = 0; i < N; i++) {
    for (j = 0; j < p; j++) {
      X(i,j) -= xmean(j);
    }
  }
  mat ymean = mean(Y);
  for (i = 0; i < N; i++) {
    for (k = 0; k < K; k++) {
      Y(i,k) -= ymean(k);
    }
  }


  // containers
  mat E (N, K, fill::zeros);

  mat mu(p, K, fill::zeros);
  if (mu0.isNotNull())  {
    mu = as<mat>(mu0);
  }

  rowvec mu_jNumer(K, fill::zeros);


  mat s (p, K, fill::ones) ;

  vec sigb(K, fill::ones);
  sigb = vectorise(var(Y, 0, 0)/2.0);
  if (sigb0.isNotNull()) {
    sigb = as<vec>(sigb0);
  }

  mat Theta(K, K, fill::zeros);
  Theta.diag() = vectorise(2.0/var(Y, 0, 0));
  if (Theta0.isNotNull()) {
    Theta = as<mat>(Theta0);
  }

  mat Sigma(K, K, fill::zeros);

  // Alpha
  mat Alpha(p, K, fill::zeros);
  Alpha.fill(0.2);
  if (Alpha0.isNotNull()) {
    Alpha = as<mat>(Alpha0);
  }

  vec alpha =  vectorise(mean(Alpha, 0));

  //Lambda
  vec Lambda(p, fill::ones);
  bool fixLambda = false;
  double lambda = 1;
  if (Lambda0.isNotNull()) {
    Lambda = as<vec>(Lambda0);
    lambda = mean(Lambda);
  }
  if (fixlambda & (lambda == 1.0))
    fixLambda = true;

  mat Beta(p, K, fill::zeros);
  // If mu_jk has initial values, then calculate.
  if (mu0.isNotNull()) {
    for (k=0; k < K; k++) {
      Beta.col(k) = Lambda % Alpha.col(k) % mu.col(k);
    }
  }
  // precomputation
  vec xx(p);
  for (j=0; j < p; j++) {
    xx(j) = sum(X.col(j) % X.col(j));
  }

  // initialization
  int it = 0;
  double L0 = -INFINITY;
  vec X_j(N, fill::zeros);
  vec ell(maxit, fill::zeros);
  mat ytilde = X * Beta;
  mat y_j(N, K, fill::zeros);
  mat y_j_new(N, K, fill::zeros);

  // algorithm
  while (it < maxit) {

    // if(it ==0){
      cout<< "mu:" << mu.submat( span(0,4), span(0,4)) << "\n";
      cout<< "Theta:" << Theta.submat( span(0,4), span(0,4)) << "\n";
      cout<< "sigb:" << sigb.rows( 0,4) << "\n";
      cout<< "Alpha:" << Alpha.submat( span(0,4), span(0,4)) << "\n";
      cout<< "a:" << alpha.rows( 0,4) << "\n";
      cout<< "s:" << s.submat( span(0,4), span(0,4)) << "\n";
    // }

    // E step
    double logitlambda = log(lambda/(1.0-lambda));
    vec logitalpha = log(alpha/(1.0-alpha));
    // update mu, s_jk and Alpha.
    for (j = 0; j < p; j++) {
      X_j = conv_to<vec>::from(X.col(j));
      //y_j = ytilde - Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
      for (k = 0; k < K; k++) {
        y_j = ytilde - Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
        if( j < 4 & k == 0 ){
          cout << y_j.submat(span(0,4), span(0,4)) << "\n";
        }
        // update s_jk
        s (j, k) = 1.0/(xx(j)*Theta(k,k) + 1.0/sigb(k));

        // update mu_jk
        vec residual_j = (Y - y_j) * Theta.col(k);
        double mu12 = sum(X.col(j) % residual_j);
        double mu3 = mu_jkThree(Theta.row(k), Alpha.row(j), mu.row(j), xx(j), k);
        mu(j, k) = (mu12 - mu3) * s(j, k);

        // update Alpha_jk
        double v_jk = logitalpha(k) + Lambda(j)/2.0 * (mu(j,k)*mu(j,k)/s(j,k) + log(s(j,k)/sigb(k)));
        Alpha(j, k) = 1.0 /(1.0+exp(-v_jk));
        ytilde = y_j + Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
      }

      // update Lambda_j
      // if (!fixLambda) {
      //   y_j_new = ytilde - Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
      //   for (k = 0; k < K; k++)	{
      //     vec residual_j = (Y - y_j_new) * Theta.col(k);
      //     double mu12 = sum(X.col(j) % residual_j);
      //     double mu3 = mu_jkThree(Theta.row(k), Alpha.row(j), mu.row(j), xx(j), k);
      //     mu_jNumer(k) = (mu12 - mu3);
      //     //ytilde = y_j + Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
      //   }
      //   double u_j2 = sum(Alpha.row(j) % mu.row(j) % mu_jNumer) +
      //     sum(Alpha.row(j) % (-mu.row(j)%mu.row(j)/s.row(j) + log(s.row(j)/sigb.t())))/2.0;
      //   double u_j3 = uThree(Theta, mu.row(j), Alpha.row(j), xx(j));
      //   double u_j  = logitlambda + u_j2 + u_j3;
      //   Lambda(j) = 1.0/(1.0+exp(-u_j));
      // }
      ytilde = y_j + Lambda(j) * kron(Alpha.row(j)%mu.row(j), X_j);
    }

    // M step
    for (k=0; k < K; k++) {
      E.col(k) = Y.col(k) - ytilde.col(k);
    }
    // update Theta^-1
    for (k=0; k < K; k++) {
      vec Siama_kk2 = Lambda%Alpha.col(k)%(mu.col(k)%mu.col(k)+s.col(k)) -
        Lambda%Lambda%Alpha.col(k)%Alpha.col(k)%mu.col(k)%mu.col(k);
      Sigma(k, k) = (sum(E.col(k) % E.col(k)) + sum(xx % Siama_kk2))/N;
      for (int t=k+1; t < K; t++) {
        vec Siama_kt2 = (Lambda - Lambda%Lambda) % Alpha.col(k) % mu.col(k) % Alpha.col(t) % mu.col(t);
        Sigma(k, t) = (sum(E.col(k) % E.col(t)) + sum(xx % Siama_kt2))/N;
        Sigma(t, k) = Sigma(k, t);
      }
    }
    Theta = inv_sympd(Sigma);

    // update sigma_{beta_k}
    for (k = 0; k < K; k++) {
      vec la = Lambda % Alpha.col(k);
      vec mu_k2 = mu.col(k) % mu.col(k);
      sigb(k) = sum(la % (mu_k2 + s.col(k))) / sum(la) ;
      alpha(k) = mean(Alpha.col(k));
    }
    if(!fixlambda) lambda = mean(Lambda);


    ell(it) = LogLikMultiple(xx, E, p, K, mu, s, Theta, sigb, Lambda, lambda, Alpha, alpha);
    // convergence
    if (ell(it) < L0){
      printf("Lowerbound decreasing, Error at iteration %d th iteration, diff=%g \n", it, ell(it) - L0);
      break;
    } else if ((ell(it) - L0) < epsStopLogLik) {
      printf("Converge at %d th iteration, LowerBound = %f \n", it, ell(it));
      break;
    }
    L0 = ell(it);

    it++;
  }



  if (it >= maxit) it = it - 1;
  vec LowerBound(it, fill::zeros);
  LowerBound = ell.subvec(0, it);

  // local FDR
  mat fdr(p, K, fill::zeros);
  fdr = 1 - Alpha.each_col() % Lambda;

  //for (k=0; k < K; k++) {
  Beta = (Lambda % Alpha.each_col()) % mu;
  //}
  // intercepts
  mat Beta0(ymean - xmean*Beta);
  // returns
  List ret ;
  ret["N"]      = N;
  ret["mu"]     = mu;
  ret["s"]      = s;
  ret["Beta0"]  = Beta0;
  ret["Beta"]   = Beta;
  ret["Alpha"]  = Alpha;
  ret["alpha"]  = alpha;
  ret["Lambda"] = Lambda;
  ret["lambda"] = lambda;
  ret["fdr"]    = fdr;
  ret["sigb"]   = sigb;
  ret["Sigma"]  = Sigma;
  ret["Theta"]  = Theta;
  ret["iter"]   = it+1;
  ret["eps"]    = ell(it) - ell(it-1);
  ret["LowerBound"] = LowerBound;

  return(ret) ;
}
