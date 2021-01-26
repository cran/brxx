#'brxx_Cor: Bayesian Estimation of Reliability from Correlation
#'
#'This function estimates reliability from a correlation
#'
#'@param x First variable.
#'@param y Second variable.
#'@param alpha Prior true score variance (covariance between tests)
#'@param beta Prior error variance (product of standard deviations minus covariance)
#'@param iter Number of iterations for the Gibbs sampler.
#'@param burn Number of samples to burn in.
#'@param seed Seed for the Gibbs sampler
#'@param CI Credible interval quantile, as a decimal (ie, for 95 percent, 0.95).
#'@param mu0 Prior means for each column.
#'@param S0 Prior variance covariance matrix.
#'@param nu0 Prior degrees of freedom for inverse Wishart prior distribution.
#'@param items Number of test items.
#'
#'@import MCMCpack
#'@import MASS
#'@import TeachingDemos
#'
#'@return Returns median posterior estimates of the variance covariance matrix.
#'
#'@examples
#'set.seed(999)
#'your_data=mvrnorm(n=15,mu=c(0,0),Sigma=matrix(c(4,5,5,9),nrow=2,ncol=2))
#'x=your_data[,1]
#'y=your_data[,2]
#'Mu0=c(0,0)
#'Sigma0=matrix(c(1,0.6,0.6,4),nrow=2,ncol=2)
#'Nu0=3-1
#'brxx_Cor(x=x,y=y,iter=5000,burn=2500,seed=999,CI=0.95,
#'mu0=Mu0,S0=Sigma0,nu0=Nu0,items=10)
#'
#'@export

brxx_Cor=function(x,y,alpha,beta,iter,burn,seed,CI,S0,nu0,mu0,items){
 data=cbind(x,y)
 filler=matrix(nrow=ncol(data),ncol=ncol(data))
 for (a in 1:ncol(data)){
   for (b in 1:ncol(data)){
     filler[a,b]=ifelse(a==b,cov(data,use="pairwise.complete.obs")[a,b],0)}}
 filler1=matrix(nrow=ncol(data),ncol=ncol(data))
 for (a in 1:ncol(data)){
   for (b in 1:ncol(data)){
     filler1[a,b]=ifelse(missing(S0),filler[a,b],S0[a,b])}}
 S0=filler1
 L0=S0
 nu0=ifelse(missing(nu0),ncol(data)*(ncol(data)+1)/2-1,nu0)
 filler2=vector(length=ncol(data))
 for (a in 1:ncol(data)){
   filler2[a]=ifelse(missing(mu0),rep(0,ncol(data)),mu0)
 }
 mu0=filler2
 n=nrow(data)
 ybar=colMeans(data,na.rm=T)
 Sigma=cov(data,use="pairwise.complete.obs")
 seed=ifelse(missing(seed),999,seed)
 iter=ifelse(missing(iter),5000,iter)
 burn=ifelse(missing(burn),iter/2,burn)
 THETA=SIGMA=NULL
 print(noquote("Sampling, this may take a minute"))

 set.seed(seed)
 pct=rep(0,iter+1)
 for(s in 1:iter)
 {
   ###Update theta
   Ln=solve(solve(L0) + n*solve(Sigma))
   mun=Ln%*%(solve(L0)%*%mu0+n*solve(Sigma)%*%ybar)
   theta=mvrnorm(1,mun,Ln)
   ###Update sigma
   Sn=S0 + (t(data)-c(theta))%*%t( t(data)-c(theta))
   Sigma=solve(rwish(nu0+n, solve(Sn)))

   ###Save results
   THETA=rbind(THETA,theta)
   SIGMA=rbind(SIGMA,c(Sigma))
   pct[s+1]=(round(s/iter*10,1))*10
   if(pct[s+1]!=pct[s]){print(noquote(paste(pct[s+1],"%")))}

 }
 CI=ifelse(missing(CI),0.95,CI)
 CI=ifelse(CI>1,CI/100,CI)
 am=ifelse(missing(alpha),1,0)
 bm=ifelse(missing(beta),10,0)
 if(am+bm==0){alpha=alpha
 beta=beta}
 if(am+bm==1){alpha=(0.6*beta)/0.4
 beta=beta}
 if(am+bm==10){alpha=alpha
 beta=0.4*alpha/0.6}
 if(am+bm==11){alpha=3.51
 beta=1.75}
 print(noquote("Calculating reliability, this may take a minute"))
 pct=rep(0,nrow(SIGMA)-burn+1)
 rxx=NULL
 items=ifelse(missing(items),1,items)
 for (s in burn:nrow(SIGMA)){
   rxx[s-burn]=(alpha+(SIGMA[s,2])*items)/
      (alpha+beta+(sqrt(SIGMA[s,1]*SIGMA[s,4]))*items)
   num=(s-burn+1)
   denom=(nrow(SIGMA)-burn)
   pct[s-burn+2]=round((num/denom)*10,1)*10
   if(pct[s-burn+2]!=pct[s-burn+1]){print(noquote(paste(pct[s-burn+2],"%")))}

 }
 rxx1=round(c(quantile(probs=c(0.5),rxx),
              emp.hpd(rxx,conf=CI)[1],
              emp.hpd(rxx,conf=CI)[2]),4)
 Out=list()
 Out$THETA=THETA
 Out$SIGMA=SIGMA
 Out$rxx_s=rxx
 Out$rxx=rxx1
 names(Out$rxx)=c("Median","LL","UL")
 return(Out)

}
