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

List updatecoord(arma::vec phi,double yy,arma::mat xx,arma::mat yx,
double lambdaupcoord, double n2,arma::mat x, std::string status){
  List out;
  //Updates  \rho  
  double yxphi = dot(phi,yx);
  double rho = (yxphi + sqrt( pow(yxphi,2)  + 4*yy*n2 ) ) / ( 2*yy  );
  int xxdim = xx.n_cols;
  //arma::vec tmp = arma::vec(phi.subvec(1,phi.n_elem-1));
//  printf("%s\n",status.c_str());
//  printf("%i\n",phi.n_elem);
//  printf("%i\n",xx.n_cols);
  arma::vec subxx = arma::vec(xx( arma::span(0,0), arma::span(1,xxdim-1)).t());
  phi(0) = ( rho * yx(0) - dot(subxx,phi.subvec(1,phi.n_elem-1))) / xx(0,0);
  if( phi.n_elem > 1){
    int philength = phi.n_elem;
    for(int j = 1; j < (philength ); j++){
      phi(j) = 0;
      double s = -1*rho*yx(j) + dot(xx.row(j),phi);
      if( s > lambdaupcoord){
        phi(j) = (lambdaupcoord - s )/( xx(j,j));
      } else if( s < -1*lambdaupcoord){
        phi(j) = (-lambdaupcoord - s )/( xx(j,j));
      }      
    }
  }
  out["phi"]=arma::vec(phi);
  out["rho"]=rho;
  return out;
}

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
    arma::mat act(p,k, arma::fill::zeros);; //To store active set
    arma::mat xbeta(n,k);
    arma::mat dnregr(n,k,arma::fill::zeros);
    List out;
    
    int i =0;
    double err1 = arma::math::inf(); //convergence of parameters
    double err2 = arma::math::inf(); //convergence of plik
    bool conv = false;
    double plik = arma::math::inf(); //sets plik to Inf
    arma::vec theta(k*p + 2*k);
    theta.fill(arma::math::inf()) ; //sets the vector of estimated parameters to Inf 
    bool warn = false;
    bool allcoord = true;
    int actiteration = 0; // act.iter diff from the one in the arguments
    double del = 0.1;
    
    arma::uvec jaime;
    
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
      while ((valuenew-valueold) > 0){ //Modify the PI probability while the logLIK is growing
        t = t*del;
        probnew = (1-t)*prob+t*probfeas; // \pi^(m+1)
        valuenew = cnloglikprob(ncomp,l1normphi,probnew,lambda,gamma);
      }
      prob = arma::vec(probnew); //Actually updates the probabilities \Pi
      
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
          arma::vec excol= ex.col(j);
          arma::mat xtilde(x);
          xtilde.each_col() %= sqrt(excol);
          arma::vec ytilde = y % sqrt(excol);
          double yy = dot(ytilde,ytilde);
          arma::vec yx = xtilde.t() * ytilde;
          arma::mat xx = xtilde.t() * xtilde; //Until now, everything the same as the R version
          arma::vec phi = beta.col(j)/ssd(j);
          double lambdaupcoord = lambda * pow(prob(j),gamma);
          List mstep = updatecoord(phi,yy,xx,yx,lambdaupcoord,sum(excol),x,"if");
          arma::vec tmp =  mstep["phi"];
          phi = tmp;
          double rho = mstep["rho"];
          jaime = phi!=0;
          act.col(j) = arma::conv_to<arma::vec>::from(phi != 0);
          beta.col(j) = arma::vec(phi/rho);
          ssd(j) = 1/rho;
        }
      } else {
        actiteration++ ;
        for(int j = 0; j<k; j++){
          arma::vec excol= ex.col(j);
          arma::vec t_act(act.col(j));
         // printf("%f\n",sum(t_act));
          ///printf("%i\n",j);
          //t_act.print();
          act.print();
          arma::mat xtilde(x.cols(arma::find(t_act>0)));
          xtilde.each_col() %= sqrt(excol);
          arma::vec ytilde = y % sqrt(excol);
          double yy = dot(ytilde,ytilde);
          arma::vec yx = xtilde.t() * ytilde;
          arma::mat xx = xtilde.t() * xtilde; //Until now, everything the same as the R version
          arma::uvec jtemp(1);
          jtemp(0) = j;
          arma::vec phi = beta(arma::find(t_act>0), jtemp) / ssd(j);
          double lambdaupcoord = lambda * pow(prob(j),gamma);
          List mstepact = updatecoord(phi,yy,xx,yx,lambdaupcoord,sum(excol),x,"else");
          arma::vec phiact =  mstepact["phi"];
          double rho = mstepact["rho"];
          act.col(j) = arma::conv_to<arma::vec>::from(phi != 0);
          beta(arma::find(t_act>0), jtemp)  = phiact/rho; //Same dimensions, the phi returned before is considering only the nonzero values
          ssd(j) = 1/rho;
        }
      }
      //E-Step:
      NumericVector value;
      xbeta = x * beta;
      for(int j = 0; j<k; j++){
        arma::vec mean(xbeta.col(j));
        //NumericVector xmean(xbeta.col(j));
        for(int h = 0; h<n; h++){
          value = dnorm(as<Rcpp::NumericVector>(wrap(y(h))), mean(h), ssd(j));
          dnregr(h,j) = value(0);
        }
      }
      arma::mat probdnregr(dnregr);
      probdnregr.each_row() %= arma::conv_to<arma::rowvec>::from(prob);
      arma::vec dmix = sum(probdnregr,1);
      ex = probdnregr;
      ex.each_col() /= dmix; //Everything ok 
      
      //Convergence criterion of \theta
      arma::vec thetaold = theta;
      theta = arma::join_cols(arma::join_cols(arma::vectorise(beta),ssd),prob);
      err1 = max(abs(theta-thetaold)/ (1+abs(theta)));
      
      //LogLig
      double loglik = sum(log(dmix));
      if(! arma::is_finite(loglik)){
        arma::get_stream_err2() << "Bad starting value, loglik = -inf" << arma::endl;
        warn = true;
        break;
      }
      
      //Convergence criterion of plik    
      double plikold = plik;
      plik = -loglik+lambda*sum(pow(prob,gamma) % arma::conv_to<arma::vec>::from(sum(abs(beta.submat(1,0,beta.n_rows-1,beta.n_cols-1)),0) ) / ssd);
      if( ( ( plik-plikold ) > 0 ) & ( i > 0 )) {
        //is plik reduced? remark: algorithm is constructed to reduce plik.
        if (warnings){
          arma::get_stream_err2() << "error: penalized negative loglik not reduced" << arma::endl;
        }
      }
      err2 = abs(plik-plikold) / (1+abs(plik));
      
      //converged?
      
      conv = ( (err1 < sqrt(term) ) & ( err2 < term ) );
      i++;     
      printf("%i\n",i);
    }
    out = List::create(conv,i,jaime,act,sum(jaime));
    return out;
  }
