library(Rcpp)
Rcpp::sourceCpp("src/fmrlasso.cpp")
n<-5
k<-4
p<-8
x<-matrix(rexp(n*p),nrow = n,ncol = p)
y<-rexp(n)
exini<-fmrlasso::ini.ex(k = k,n = n)

(tmp<-fmrlasso(x = x,y = y,k = k,lambda = 1.08,ssdini = 0.5,
         gamma = 0.9,exini = exini))

#Comparison implementations cnloglikprob
library(microbenchmark)
#ncomp,l1normphi,probfeas,lambda,gamma
microbenchmark(times = 1000,
  cnloglikprob(tmp[[1]],tmp[[2]],tmp[[3]],1.28,0.9),
  fmrlasso::cnloglikprob(tmp[[1]],tmp[[2]],tmp[[3]],1.28,0.9)
)

#phi=phi,yy=yy,xx=xx,yx=yx,lambda=lambdaupcoord,n=sum(excol),x =x
microbenchmark(times = 10000,
updatecoord2(phi = tmp[[1]],yy = tmp[[2]],xx = tmp[[3]],yx = tmp[[4]],lambda = tmp[[5]],n = tmp[[6]],x = tmp[[7]]),
updatecoord(phi = tmp[[1]],yy = tmp[[2]],xx = tmp[[3]],yx = tmp[[4]],lambda = tmp[[5]],n = tmp[[6]],x = tmp[[7]])
)
