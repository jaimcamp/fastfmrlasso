//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
NumericVector rowSumsC(NumericMatrix x){
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);
  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
    }
    out[i] = total;
  }
  return out;
}

NumericVector colSumsC(NumericMatrix x){
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(ncol);
  for (int i = 0; i < ncol; i++) {
    double total = 0;
    for (int j = 0; j < nrow; j++) {
      total += x(j, i);
    }
    out[i] = total;
  }
  return out;
}

NumericVector cnloglikprob(
  NumericVector ncomp, NumericVector l1normphi,
  NumericVector prob, double lambda, double gamma){
    // Purpose: complete negative loglikelihood involving (pi1,...,pik) (rho,phi fixed).
    return -sum(ncomp*log(prob))+lambda*sum(pow(prob,gamma)*l1normphi);
    }
  
  

  
//cnloglikprob <- function(ncomp,l1normphi,prob,lambda,gamma=1)
//{
//  ## Purpose: complete negative loglikelihood involving (pi1,...,pik) (rho,phi fixed).
//  ## ----------------------------------------------------------------------
//  ## Arguments:         
//  ## ----------------------------------------------------------------------
//  ## Author: Nicolas Staedler
//  -sum(ncomp*log(prob))+lambda*sum((prob)^{gamma}*l1normphi)
//}

// [[Rcpp::export]]
List fmrlasso(
  NumericMatrix x, NumericVector y,
  int k, double lambda, double ssdini, NumericMatrix exini,
  double gamma=1,
  double term= 10e-6, int maxiter=1000,
  int actiter=10,
  bool warnings=true){
    int n = y.size();
    int p = x.ncol();
    NumericVector prob(k , 1.0/k);
    arma::mat beta(p,k, arma::fill::zeros);
    arma::vec ssd(k);
    ssd.fill(ssdini);
    NumericMatrix ex = exini;
    List act; //To store active set
    NumericMatrix xbeta(n,k);
    NumericMatrix dnregr(n,k);
    
    int i =0;
    double err1 = arma::math::inf(); //convergence of parameters
    double err2 = arma::math::inf(); //convergence of plik
    bool conv = false;
    double plik = arma::math::inf(); //sets plik to Inf
    double theta = arma::math::inf(); //sets the vector of estimated parameters to Inf 
    bool warn = false;
    bool allcoord = true;
    int actiteration = 0; // act.iter diff from the one in the arguments
    double del = 0.1;
    
    while ( ( (!conv)|(!allcoord) ) & (i<maxiter) & !warn ){
      //while  conv or allcord are false AND i is less than maxiter AND warn is false
      //M-STEP
      //update prob
      NumericVector ncomp = colSumsC(ex);
      beta(0,0)=-8;
      beta(2,1)=-13;
      beta(2,2)=23;
      arma::mat temp = sum(abs(beta.submat(1,0,beta.n_rows-1,beta.n_cols-1)),0);
      //NumericMatrix temp = beta( Range(1,beta.nrow()-1), Range(0,beta.ncol()-1));
      //int tempsize = temp.ncol()*temp.nrow();
      //for(int j=0; j<tempsize; j++){
      //  temp[j] = abs(temp[j]);
      //}
      arma::mat temp2 = temp % (1 / ssd);
      //NumericVector l1normphi = colSumsC(temp/ssd);
      //NumericVector probfeas = ncomp / n; //feasible point
      //NumericVector valueold = cnloglikprob(ncomp,l1normphi,prob,lambda,gamma);
      //NumericVector valueold2 = ncomp * l1normphi;
      List out = List::create(beta, temp,temp2,ssd);
      return out;
      warn = true;
    }
  }
