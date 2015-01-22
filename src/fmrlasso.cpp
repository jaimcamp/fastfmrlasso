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

List updatecoord(int a){
  return List::create(a+1);
  
}
//updatecoord <- function(phi,yy=yy,xx=xx,yx=yx,lambda=lambda,n=n,x=x){
//
//  #update rho
//  yxphi <- sum(yx*phi)
//  rho <- (yxphi+sqrt(yxphi^{2}+4*yy*n))/(2*yy)
//  #update phi[1] (not penalized)
//  phi[1] <- (rho*yx[1]-sum(xx[1,-1]*phi[-1]))/xx[1,1]
//
//  if (length(phi)>1){
//    for (j in 2:length(phi)){
//      phi[j] <- 0
//      s <- -rho*yx[j] + sum(xx[j,]*phi)
//      if (s > lambda){
//        phi[j] <- (lambda-s)/(xx[j,j])
//      }
//      if (s < -lambda){
//        phi[j] <- (-s-lambda)/(xx[j,j])
//      }
//    }
//  }
//  list(phi=phi,rho=rho)
//}


// [[Rcpp::export]]

double cnloglikprob(
  arma::vec ncomp, arma::vec l1normphi, arma::vec prob, double lambda, double gamma){
    // Purpose: complete negative loglikelihood involving (pi1,...,pik) (rho,phi fixed).
    return (-1*dot(ncomp,log(prob))) + (lambda * dot(pow(prob,gamma),l1normphi));
}
    
// [[Rcpp::export]]
    
List fmrlasso(
  arma::mat x, arma::vec y,
  int k, double lambda, double ssdini, arma::mat exini,
  double gamma=1,
  double term= 10e-6, int maxiter=1000,
  int actiter=10,
  bool warnings=true){
    int n = y.n_elem;
    int p = x.n_cols;
    arma::vec prob(k);
    prob.fill(1.0/k);
    arma::mat beta(p,k, arma::fill::zeros);
    arma::vec ssd(k);
    ssd.fill(ssdini);
    arma::mat ex = exini;
    List act; //To store active set
    NumericMatrix xbeta(n,k);
    NumericMatrix dnregr(n,k);
    List out;
    
    int i =0;
    //double err1 = arma::math::inf(); //convergence of parameters
    //double err2 = arma::math::inf(); //convergence of plik
    bool conv = false;
    //double plik = arma::math::inf(); //sets plik to Inf
    //double theta = arma::math::inf(); //sets the vector of estimated parameters to Inf 
    bool warn = false;
    bool allcoord = true;
    int actiteration = 0; // act.iter diff from the one in the arguments
    double del = 0.1;
    
    while ( ( (!conv)|(!allcoord) ) & (i<maxiter) & !warn ){
      //while  conv or allcord are false AND i is less than maxiter AND warn is false
      //M-STEP
      //update prob
      arma::vec ncomp = arma::conv_to<arma::vec>::from(sum(ex,0));
      beta(0,0)=-8; //To try
      beta(2,1)=-13;
      beta(2,2)=23;
      arma::mat temp = sum(abs(beta.submat(1,0,beta.n_rows-1,beta.n_cols-1)),0);
      arma::vec l1normphi = arma::conv_to<arma::vec>::from(temp)  % (1 / ssd);
      arma::vec probfeas = ncomp / n; //feasible point
      double valueold = cnloglikprob(ncomp,l1normphi,prob,lambda,gamma);
      double valuenew = cnloglikprob(ncomp,l1normphi,probfeas,lambda,gamma);
      double t = 1.0;
      arma::vec probnew(probfeas);
      out = List::create(ncomp,l1normphi,probfeas,lambda,gamma);
      printf("%lf\n", t);
      while ((valuenew-valueold) > 0){ //Modify the PI probability while the logLIK is growing
        t = t*del;
        probnew = (1-t)*prob+t*probfeas; // \pi^(m+1)
        valuenew = cnloglikprob(ncomp,l1normphi,probnew,lambda,gamma);
        printf("%lf\n", t);
      }    
      // Update phi,rho
      if ( (allcoord) & (i>0) ){
      allcoord = false;
      }
      
      if ( (actiteration==actiter) | conv ){
      actiteration = 0;
      allcoord = true;
      }
      
      printf("%s\n", allcoord ? "true" : "false");
      
      if (allcoord){
        for (int j=0; j < k; j++){
          printf("%s\n", "IF");
          arma::vec excol= ex.col(j);
          arma::mat xtilde(x);
          xtilde.each_col() %= sqrt(excol);
          arma::vec ytilde = y % sqrt(excol);
          arma::vec yy = sum(pow(ytilde,2));
          double yy2 = dot(ytilde,ytilde);
          arma::vec yx = xtilde.t() * ytilde;
          arma::mat xx = xtilde.t() * xtilde;
          beta.col(j)/ssd[j];
          List mstep = updatecoord(2);
                    
//          mstep <- updatecoord(phi=beta[,j]/ssd[j],yy=yy,xx=xx,yx=yx,lambda=lambda*(prob[j])^{gamma},n=sum(EX))
//          phi <- mstep$phi
//          rho <- mstep$rho
//          act[[j]] <-which(phi!=0)
//          beta[,j] <- phi/rho
//          ssd[j] <- 1/rho
          out = List::create(ex,excol,x,xtilde,ytilde,yy,yy2);
        }
      } else {
        printf("%s\n", "Else");
      }
          
 
      return out;
      warn = true;
    }
  }
