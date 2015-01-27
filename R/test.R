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
fmrlasso::fmrlasso(x = x,y = y,k = k,lambda = 1.08,  ssd.ini = 0.5,
                   gamma = 0.9,ex.ini = exini)

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

##Debuggin and testing:

save(tmp,file="R/test_cnlog_files.RData")
load(file="R/test_cnlog_files.RData")
cnloglikprob(tmp[[1]],tmp[[2]],tmp[[3]],1.28,0.9)
fmrlasso::cnloglikprob(tmp[[1]],tmp[[2]],tmp[[3]],1.28,0.9) #Alles gut!

##Testing time
microbenchmark(times = 10000,
               tmp<-fmrlasso(x = x,y = y,k = k,lambda = 1.08,ssdini = 0.5,gamma = 0.9,exini = exini),
               tmp2<-fmrlasso::fmrlasso(x = x,y = y,k = k,lambda = 1.08,  ssd.ini = 0.5,gamma = 0.9,ex.ini = exini)
)
#min        lq      mean    median        uq       max neval
#7.495873  7.755297  8.005111  7.899313  8.087654  53.49411 10000
#55.971208 57.939281 59.361686 58.792932 59.841854 115.67189 10000
 
#Something more complex 

n<-10
k<-4
p<-80
x<-matrix(rexp(n*p),nrow = n,ncol = p)
y<-rexp(n)
exini<-fmrlasso::ini.ex(k = k,n = n)
microbenchmark(times =1000,
               tmp<-fmrlasso(x = x,y = y,k = k,lambda = 1.08,ssdini = 0.5,gamma = 1,exini = exini),
               tmp2<-fmrlasso::fmrlasso(x = x,y = y,k = k,lambda = 1.08,  ssd.ini = 0.5,gamma = 1,ex.ini = exini)
)
#mean    median        uq       max neval
#67.62851  66.02392  67.34028  141.8617  1000
#800.94941 791.33235 813.40691 1123.8950  1000
> 
