######################################
##Coordinate Update for scaled LASSO##
######################################

updatecoord <- function(phi,yy=yy,xx=xx,yx=yx,lambda=lambda,n=n,x=x){

  #update rho
  yxphi <- sum(yx*phi)
  rho <- (yxphi+sqrt(yxphi^{2}+4*yy*n))/(2*yy)
  #update phi[1] (not penalized)
  phi[1] <- (rho*yx[1]-sum(xx[1,-1]*phi[-1]))/xx[1,1]

  if (length(phi)>1){
    for (j in 2:length(phi)){
      phi[j] <- 0
      s <- -rho*yx[j] + sum(xx[j,]*phi)
      if (s > lambda){
        phi[j] <- (lambda-s)/(xx[j,j])
      }
      if (s < -lambda){
        phi[j] <- (-s-lambda)/(xx[j,j])
      }
    }
  }
  list(phi=phi,rho=rho)
}


############################
##Initialisation of E-STEP##
############################

ini.ex <- function(k,n,weight=0.9){
# k: nr. of components
# n: nr. of observations
  
  EX <- matrix(NA,nrow=n,ncol=k)
  ex <- sample(1:k,n,replace=TRUE)
  for (i in 1:n){
    EX[i,ex[i]] <- weight
    EX[i,-ex[i]] <- 1-weight
    EX[i,] <- EX[i,]/sum(EX[i,])
  }
  EX
}


##########################
##LAMBDAMAX scaled LASSO##
##########################

lambdamax <- function(y,x){
# x: design matrix (with intercept)
# y: dependent variable
  x <- x[,-1]
  n <- length(y)
  rho <- sqrt(n/sum((y-mean(y))^{2}))
  mu <- mean(y)*rho
  grad <- -rho*crossprod(y,x)+mu*n*apply(x,2,mean)
  max(abs(grad))
}


#######################################################################
## FMRLASSO                                                          ##
## PENALTY=la*[pi1^{gamma}*|beta1/ssd1|+...+pik^{gamma}*|betak/ssdk|]##
#######################################################################

cnloglikprob <- function(ncomp,l1normphi,prob,lambda,gamma=1)
{
  ## Purpose: complete negative loglikelihood involving (pi1,...,pik) (rho,phi fixed).
  ## ----------------------------------------------------------------------
  ## Arguments:         
  ## ----------------------------------------------------------------------
  ## Author: Nicolas Staedler
  -sum(ncomp*log(prob))+lambda*sum((prob)^{gamma}*l1normphi)
}


fmrlasso <- function(x,y,k,lambda,gamma=1,ssd.ini,ex.ini,term=10^{-6},maxiter=1000,actiter=10,warnings=TRUE){

# ARGUMENTS:
# x: design matrix (with intercept)
# y: response vector
# k: nr of components
# lambda: 
# gamma: 
# ssd.ini: initialisation for variances (!=0)
# ex.ini: initialisation of aposterior probabilies (n*k matrix)
# T: termination criterion (see Philip Gill)
# maxiter: maximal nr of iterations
# warnings: should warnings (plik not reduced, loglik=-INF) be showed?
# actiter=10: update active set every 10th EM-iteration
  
  n <- length(y)
  p <- dim(x)[2]
  prob <- rep(1/k,k)
  beta <- matrix(0,nrow=p,ncol=k)
  ssd <- rep(ssd.ini,k)
  ex <- ex.ini
  act <- list()#to store the active set
  xbeta <- matrix(0,nrow=n,ncol=k)
  dnregr <- matrix(0,nrow=n,ncol=k)
  
  i <- 0
  err1 <- Inf #convergence of parameters
  err2 <- Inf #convergence of plik
  conv <- FALSE
  plik <- Inf #sets plik to Inf
  theta <- Inf #sets the vector of estimated parameters to Inf 
  warn <- FALSE
  allcoord <- TRUE
  act.iter <- 0
  del <- 0.1

  while (((!conv)|(!allcoord))&(i<maxiter)&!warn){

#M-STEP
    #update prob
    ncomp <- colSums(ex)
    l1normphi <- colSums(abs(beta[-1,,drop=FALSE]))/ssd
    probfeas <- ncomp/n #feasible point
    
    value.old <- cnloglikprob(ncomp,l1normphi,prob,lambda,gamma=gamma)
    value.new <- cnloglikprob(ncomp,l1normphi,probfeas,lambda,gamma=gamma)
    t <- 1
    probnew <- probfeas
    while ((value.new-value.old) > 0){
      t <- t*del
      probnew <- (1-t)*prob+t*probfeas
      value.new <- cnloglikprob(ncomp,l1normphi,probnew,lambda,gamma=gamma)
    }
    prob <- probnew

    #update phi,rho

    if ((allcoord)&(i>0)){
      allcoord <- FALSE
    }

    if ((act.iter==actiter)|conv){
      act.iter <- 0
      allcoord <- TRUE
    }
 
    if (allcoord){
      for (j in 1:k){
        EX <- ex[,j]
        xtilde <- sqrt(EX)*x
        ytilde <- sqrt(EX)*y
        yy <- sum(ytilde^{2})
        yx <- crossprod(ytilde,xtilde)
        xx <- crossprod(xtilde)
        mstep <- updatecoord(phi=beta[,j]/ssd[j],yy=yy,xx=xx,yx=yx,lambda=lambda*(prob[j])^{gamma},n=sum(EX))
        phi <- mstep$phi
        rho <- mstep$rho
        act[[j]] <-which(phi!=0)
        beta[,j] <- phi/rho
        ssd[j] <- 1/rho
      }
    }else{
      act.iter <- act.iter+1
      for (j in 1:k){
        EX <- ex[,j]
        t.act <- act[[j]]
        xtilde <- sqrt(EX)*x[,t.act,drop=FALSE]
        ytilde <- sqrt(EX)*y
        yy <- sum(ytilde^{2})
        yx <- crossprod(ytilde,xtilde)
        xx <- crossprod(xtilde)
        mstepact <- updatecoord(phi=beta[t.act,j]/ssd[j],yy=yy,xx=xx,yx=yx,lambda=lambda*(prob[j])^{gamma},n=sum(EX))
        phiact <- mstepact$phi
        rho <- mstepact$rho
        beta[t.act,j] <- phiact/rho
        ssd[j] <- 1/rho
      }
    }
    
   
#E-STEP
    xbeta <- x%*%beta
    for (j in 1:k){
      dnregr[,j] <- dnorm(y,mean=xbeta[,j],sd=ssd[j])
    }
    probdnregr <- t(prob*t(dnregr))
    dmix <- rowSums(probdnregr)
    ex <- probdnregr/dmix

#convergence criterion of theta
    thetaold <- theta
    theta <- c(as.vector(beta),as.vector(ssd),prob)
    err1 <- max(abs(theta-thetaold)/(1+abs(theta)))#Philip Gill termination criteria

#loglik    
    loglik <- sum(log(dmix))
    if (loglik==-Inf){warning("bad starting value (loglik=-Inf)");warn <- TRUE;break}

#convergence criterion of plik    
    plikold <- plik
    plik <- -loglik+lambda*sum((prob)^{gamma}*colSums(abs(beta[-1,,drop=FALSE]))/ssd)
    if(((plik-plikold)>0)&(i>0)) {  #is plik reduced? remark: algorithm is constructed to reduce plik.
     if (warnings) warning("error: penalized negative loglik not reduced")
    }
    err2 <- abs(plik-plikold)/(1+abs(plik))

    #converged?
    conv <- ((err1<sqrt(term))&(err2<term))
 
    i <- i+1
  
  }
  
  n.zero <- sum(beta==0)
  d <- k*(p+1+1)-1-n.zero #degree`s of freedom (sh DeSarbo)
  bic <- -2*loglik+log(n)*d #BIC criterion

  clust <- apply(-ex,1,which.min)

  #prepare result
  dimnames(beta) <- list(coef=0:(p-1),comp=1:k)
  dimnames(beta)[[1]][1] <- 'intercept'
  dimnames(ex) <- list(NULL,comp=1:k)
  
  res <- list(k=k,prob=prob,coef=beta,ssd=ssd,plik=plik,bic=bic,ex=ex,cluster=clust,niter=i,warnings=warn)
  
  res
  
}


##########
##TRACES##
##########

fmrlassopath <- function(x,y,k,lambda,gamma=1,ssd.ini,ex.ini,term=10^{-6},maxiter=1000,actiter=10){

# ARGUMENTS:
# x: design matrix (with intercept)
# y: response vector
# k: nr of components  
# lambda:
# gamma:
# ssd.ini: initialisation for variances (!=0)
# ex.ini: initialisation of aposterior probabilies (n*k matrix)
# T: termination criterion
# maxiter: maximal nr of iterations

  n <- length(y)
  p <- dim(x)[2]-1 
  l.la <- length(lambda)

  prob <- matrix(0,nrow=k,ncol=l.la)
  coef <- array(0,dim=c(p+1,k,l.la))
  ssd <- matrix(0,nrow=k,ncol=l.la)
  plik <- rep(0,l.la)
  bic <- rep(0,l.la)
  ex <- array(0,dim=c(n,k,l.la))
  cluster <- matrix(0,nrow=n,ncol=l.la)
  niter <- rep(0,l.la)

  for (i in 1:l.la){

    fit <- fmrlasso(x,y,k,lambda=lambda[i],gamma=gamma,ssd.ini=ssd.ini,ex.ini=ex.ini,term=term,maxiter=maxiter,actiter=actiter)
    prob[,i] <- fit$prob
    coef[,,i] <- fit$coef
    ssd[,i] <- fit$ssd
    plik[i] <- fit$plik
    bic[i] <- fit$bic
    ex[,,i] <- fit$ex
    cluster[,i] <- fit$cluster
    niter[i] <- fit$niter
    #cat('lambda:',lambda[i],'\n')
  }

  #prepare results
  dimnames(coef) <- list(coef=0:p,comp=1:k,lambda=lambda)
  dimnames(coef)[[1]][1] <- 'intercept'
  dimnames(prob) <- dimnames(ssd) <- list(comp=1:k,lambda=lambda)
  dimnames(ex) <- list(NULL,comp=1:k,lambda=lambda)
  dimnames(cluster) <- list(NULL,lambda=lambda)
  
  res <- list(k=k,lambda=lambda,prob=prob,coef=coef,ssd=ssd,plik=plik,bic=bic,ex=ex,cluster=cluster,niter=niter)
  class(res) <- 'fmrlassopath'
  res

}

###############
##PLOT TRACES##
###############
plotfmrlassopath <- function(object,col='black'){
  
  lambda <- object$lambda
  prob <- object$prob
  coef <- object$coef
  ssd <- object$ssd
  bic <- object$bic
  k <- dim(coef)[2]
  l.la <- length(lambda)
  minbic <- which.min(bic)
  
  par(mfrow=c(ceiling((k+3)/2),2),oma=c(0.5,0.5,4,0.5),mar=c(3,3,1,1),mgp=c(1.2,0.3,0),tck=-0.01)#mult.fig(k+3)
  
  for (i in 1:k){
    matplot(lambda,t(coef[,i,]),type="l",xlab=expression(lambda),ylab=paste("beta",i,sep=' '),col=col,pch='');abline(v=lambda[minbic],lty=4)
  }
  matplot(lambda,t(prob),type="l",xlab=expression(lambda),ylab="prior probability");abline(v=lambda[minbic],lty=4)
  text(lambda[l.la-1],t(prob)[l.la-1,],1:k)
  matplot(lambda,t(ssd),type="l",xlab=expression(lambda),ylab="ssd");abline(v=lambda[minbic],lty=4)
  text(lambda[l.la-1],t(ssd)[l.la-1,],1:k)
  plot(lambda,bic,type="l",xlab=expression(lambda),ylab="bic");abline(v=lambda[minbic],lty=4)
  mtext('estimated parameters',side=3,line=2,outer=TRUE)
}

###########
##LOGLOSS##
###########

logloss <- function(coef,ssd,prob,x,y){
  k <- length(prob)
  n <- length(y)
  dnregr <- matrix(0,nrow=n,ncol=k)
  xcoef <- x%*%coef
  for (j in 1:k){
    dnregr[,j] <- dnorm(y,mean=xcoef[,j],sd=ssd[j])
  }
  probdnregr <- t(prob*t(dnregr))
  dmix <- rowSums(probdnregr)
  -2*sum(log(dmix))
}

#################################################
##PREDICTION ERROR FOR TEST DATA: X.TEST,Y.TEST##
#################################################

predloss <- function(model,x,y){
  
  l <- length(model$lambda)
  coef <- model$coef
  ssd <- model$ssd
  prob <- model$prob
  loss <- rep(0,l)
  for (i in 1:l){
    loss[i] <- logloss(as.matrix(coef[,,i]),ssd[,i],prob[,i],x,y)
  }
  pred <- list(loss=loss,indexopt=which.min(loss))
  pred
}

##################
##cvfmrlassopath##
##################

cv.folds <- function (n, folds = 10) 
{
    split(sample(1:n), rep(1:folds, length = n))
}

cvfmrlassopath <- function(x,y,k,lambda,gamma=1,ssd.ini=0.5,ex.ini,term=10^{-6},maxiter=1000,actiter=10,K=10){
  all.folds <- cv.folds(length(y),K)
  errmat <- matrix(0,length(lambda),K)
  for (i in seq(K)){
    omit <- all.folds[[i]]
    fit <- fmrlassopath(x=x[-omit,],y=y[-omit],k=k,lambda=lambda,gamma=gamma,ssd.ini=ssd.ini,ex.ini=as.matrix(ex.ini[-omit,]),term=term,maxiter=maxiter,actiter=actiter)
    errmat[,i] <- predloss(fit,x[omit,],y[omit])$loss
    cat("\n CV Fold", i, "\n\n")
  }
  cv <- apply(errmat, 1, mean)
  cv.error <- sqrt(apply(errmat, 1, var)/K)
  dimnames(errmat) <- list(lambda=lambda,fold=1:K)
  
  object <- list(lambda = lambda, cv = cv, cv.error = cv.error,errmat=errmat)
  invisible(object)
}


