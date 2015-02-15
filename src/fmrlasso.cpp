//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

List updatecoord_f(arma::vec phi,double yy,arma::mat xx,arma::mat yx,
double lambdaupcoord, double n2,arma::mat x, std::string status){
  List out;
  //Updates  \rho  
  double yxphi = dot(phi,yx);
  double rho = (yxphi + sqrt( pow(yxphi,2)  + 4*yy*n2 ) ) / ( 2*yy  );
  if( phi.n_elem <=1){
    phi(0) = ( rho * yx(0) ) / xx(0,0);
  } else {
    //xx.submat(0,1,xx.n_rows-1,xx.n_cols-1)
    //arma::vec tmp = arma::vec(phi.subvec(1,phi.n_elem-1));
    int xxdim = xx.n_cols;
    arma::vec subxx = arma::vec(xx( arma::span(0,0), arma::span(1,xxdim-1)).t());
    phi(0) = ( rho * yx(0) - dot(subxx,phi.subvec(1,phi.n_elem-1))) / xx(0,0);
  }
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

double cnloglikprob_f(
  arma::vec ncomp, arma::vec l1normphi, arma::vec prob, double lambda, double gamma){
    // Purpose: complete negative loglikelihood involving (pi1,...,pik) (rho,phi fixed).
    return (-1*dot(ncomp,log(prob))) + (lambda * dot(pow(prob,gamma),l1normphi));
}
    
// [[Rcpp::export]]
    
List fmrlasso_f(
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
    double loglik = 0;
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
    
   while ( ( (!conv)|(!allcoord) ) & (i<maxiter) & !warn ){
      //while  conv or allcord are false AND i is less than maxiter AND warn is false
      //M-STEP
      //update prob
      arma::vec ncomp = arma::conv_to<arma::vec>::from(sum(ex,0));
      arma::mat temp = sum(arma::abs(beta.submat(1,0,beta.n_rows-1,beta.n_cols-1)),0);
      arma::vec l1normphi = arma::conv_to<arma::vec>::from(temp)  % (1 / ssd);
      arma::vec probfeas = ncomp / n; //feasible point
      double valueold = cnloglikprob_f(ncomp,l1normphi,prob,lambda,gamma);
      double valuenew = cnloglikprob_f(ncomp,l1normphi,probfeas,lambda,gamma);
      double t = 1.0;
      arma::vec probnew(probfeas);
      while ((valuenew-valueold) > 0){ //Modify the PI probability while the logLIK is growing
        t = t*del;
        probnew = (1-t)*prob+t*probfeas; // \pi^(m+1)
        valuenew = cnloglikprob_f(ncomp,l1normphi,probnew,lambda,gamma);
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
      
      //printf("%s\n", allcoord ? "true" : "false");
      
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
          List mstep = updatecoord_f(phi,yy,xx,yx,lambdaupcoord,sum(excol),x,"if");
          arma::vec tmp =  mstep["phi"];
          phi = tmp;
          double rho = mstep["rho"];
          act.col(j) = arma::conv_to<arma::vec>::from(phi != 0);
          beta.col(j) = arma::vec(phi/rho);
          ssd(j) = 1/rho;
        }
      } else {
        actiteration++ ;
        for(int j = 0; j<k; j++){
          arma::vec excol= ex.col(j);
          arma::vec t_act(act.col(j));
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
          List mstepact = updatecoord_f(phi,yy,xx,yx,lambdaupcoord,sum(excol),x,"else");
          arma::vec phiact =  mstepact["phi"];
          double rho = mstepact["rho"];
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
      //xbeta.print("xbeta");
      //y.print("y");
      //ssd.print("ssd");
      arma::mat probdnregr(dnregr);
      probdnregr.each_row() %= arma::conv_to<arma::rowvec>::from(prob);
      arma::vec dmix = sum(probdnregr,1);
      ex = probdnregr;
      ex.each_col() /= dmix; //Everything ok 
      
      //Convergence criterion of \theta
      arma::vec thetaold = theta;
      theta = arma::join_cols(arma::join_cols(arma::vectorise(beta),ssd),prob);
      err1 = max(arma::abs(theta-thetaold)/ (1+arma::abs(theta)));
      
      //LogLig
      loglik = sum(log(dmix));
      if(! arma::is_finite(loglik)){
        arma::get_stream_err2() << "Bad starting value, loglik = -inf" << arma::endl;
        warn = true;
        break;
      }
      
      //Convergence criterion of plik    
      double plikold(plik);
      plik = -loglik+lambda*sum(pow(prob,gamma) % arma::conv_to<arma::vec>::from(sum(arma::abs(beta.submat(1,0,beta.n_rows-1,beta.n_cols-1)),0) ) / ssd);
      if( ( ( plik-plikold ) > 0 ) & ( i > 0 )) {
        //is plik reduced? remark: algorithm is constructed to reduce plik.
        if (warnings){
          arma::get_stream_err2() << "error: penalized negative loglik not reduced" << arma::endl;
        }
      }
      err2 = std::abs(plik-plikold) / (1+std::abs(plik));
      
      //converged?
      conv = ( (err1 < sqrt(term) ) & ( err2 < term ) );
      i++;     
      //printf("%i\n",i);
    }
    double n_zero = accu(beta == 0);
    double d = k*(p+1+1) - 1 - n_zero; //Degrees of freedom
    double bic = -2 * loglik + log(n)*d ; //BIC criterion
    
    arma::uvec cluster(n);
    arma::rowvec exrow;
    for (int j = 0; j<n; j++){
      exrow = ex.row(j);
      exrow.max(cluster(j));
    }
    cluster = cluster +1;
    out = List::create(Rcpp::Named("k")=k,Rcpp::Named("prob")=prob,
    Rcpp::Named("coef")=beta,Rcpp::Named("ssd")=ssd,
    Rcpp::Named("plik")=plik,Rcpp::Named("bic")=bic,
    Rcpp::Named("ex")=ex,Rcpp::Named("cluster")=cluster,
    Rcpp::Named("niter")=i,Rcpp::Named("warnings")=warn);
    return out;
  }



//List cvfmrlassopath(
//  arma::mat x, arma::vec y,
//  int k, arma::vec lambda, double ssdini, arma::mat exini,
//  double gamma=1,
//  double term= 10e-6, int maxiter=1000,
//  int actiter=10, int K = 10) {
//  //all.folds <- fmrlasso::cv.folds(length(y), K)
//  arma::mat errmat(lambda.n_elem, k, arma::fill::zeros);; 
//  for(int i = 1; i < K; i++){
//    omit <- all.folds[[i]]
//    fit <- fmrlassopath_f(x = x[-omit, ], y = y[-omit], k = k, 
//                          lambda = lambda, gamma = gamma, ssd.ini = ssd.ini, 
//                          ex.ini = as.matrix(ex.ini[-omit, ]), term = term, 
//                          maxiter = maxiter, actiter = actiter)
//    errmat[, i] <- fmrlasso::predloss(fit, x[omit, ], y[omit])$loss
//    cat("\n CV Fold", i, "\n\n")
//  }
//  cv <- apply(errmat, 1, mean)
//  cv.error <- sqrt(apply(errmat, 1, var)/K)
//  dimnames(errmat) <- list(lambda = lambda, fold = 1:K)
//  object <- list(lambda = lambda, cv = cv, cv.error = cv.error, 
//                 errmat = errmat)
//  invisible(object)
//}

// [[Rcpp::export]]

List fmrlassopath_f(
  arma::mat x, arma::vec y,
  int k, arma::vec lambda, double ssdini, arma::mat exini,
  double gamma=1,
  double term= 10e-6, int maxiter=1000,
  int actiter=10) {
    int n = y.n_elem;
    int p = x.n_cols -1;
    int lla = lambda.n_elem;
    arma::mat prob(k, lla, arma::fill::zeros);
    arma::cube coef(p+1, k, lla, arma::fill::zeros);
    arma::mat ssd(k, lla, arma::fill::zeros);
    arma::vec plik(lla,arma::fill::zeros);
    arma::vec bic(lla,arma::fill::zeros);
    arma::cube ex(n, k, lla, arma::fill::zeros);
    arma::mat cluster(n, lla, arma::fill::zeros);
    arma::vec niter(lla, arma::fill::zeros);
    List fit;
    for (int i=0; i<lla; i++) {
      fit = fmrlasso_f(x, y, k, lambda(i), ssdini, exini, gamma, term, 
                        maxiter, actiter);
      arma::vec tmp = fit["prob"];
      prob.col(i) = tmp;
      arma::mat tmp2 =fit["coef"];
      coef.slice(i) =  tmp2;
      arma::vec tmp3 = fit["ssd"];
      ssd.col(i) = tmp3;
      plik(i) = fit["plik"];
      bic(i) = fit("bic");
      arma::mat tmp4 = fit["ex"];
      ex.slice(i) = tmp4;
      arma::vec tmp5 = fit["cluster"];
      cluster.col(i) = tmp5;
      niter(i) = fit("niter");
    }
//    dimnames(coef) <- list(coef = 0:p, comp = 1:k, lambda = lambda)
//    dimnames(coef)[[1]][1] <- "intercept"
//    dimnames(prob) <- dimnames(ssd) <- list(comp = 1:k, lambda = lambda)
//    dimnames(ex) <- list(NULL, comp = 1:k, lambda = lambda)
//    dimnames(cluster) <- list(NULL, lambda = lambda)
    List res = List::create(Rcpp::Named("k")=k, Rcpp::Named("lambda") = lambda,  Rcpp::Named("prob") = prob,  Rcpp::Named("coef") = coef, 
                 Rcpp::Named("ssd") = ssd,  Rcpp::Named("plik") = plik,  Rcpp::Named("bic") = bic,  Rcpp::Named("ex") = ex,  Rcpp::Named("cluster") = cluster, 
                 Rcpp::Named("niter") = niter);
    return(res);
                 
}