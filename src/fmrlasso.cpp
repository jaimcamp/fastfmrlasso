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

// [[Rcpp::export]]
List fmrlasso(
  NumericMatrix x, NumericVector y,
  int k, double lamda, double ssdini, NumericMatrix exini,
  double gamma=1,
  double term= 10e-6, int maxiter=1000,
  int actiter=10,
  bool warnings=true){
    int n = y.size();
    int p = x.ncol();
    NumericVector prob(k , 1.0/k);
    NumericMatrix beta(p,k);
    NumericVector ssd(k,ssdini);
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
      NumericMatrix temp = beta( Range(0,beta.nrow()-1), Range(1,beta.ncol()-1));
      int tempsize = temp.ncol()*temp.nrow();
      for(int j=0; j<tempsize; j++){
        temp[j] = abs(temp[j]);
      }
      NumericVector l1normphi = colSumsC(temp);
      NumericVector probfeas = ncomp / n; //feasible point
      List out = List::create(ex,beta,temp,l1normphi,probfeas);
      return out;
      warn = true;
    }
  }
