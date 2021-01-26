#'bcov: Bayesian Estimation of the Variance Covariance Matrix
#'
#'This function estimates the variance covariance matrix for a
#'
#'@param data N by P data matrix.
#'@param iter Number of iterations for the Gibbs sampler.
#'@param burn Number of samples to burn in.
#'@param seed Seed for the Gibbs sampler
#'@param CI Credible interval quantile, as a decimal (ie, for 95 percent, 0.95).
#'@param mu0 Prior means for each column.
#'@param S0 Prior variance covariance matrix.
#'@param nu0 Prior degrees of freedom for inverse Wishart prior distribution.
#'
#'@import MCMCpack
#'@import MASS
#'@import TeachingDemos
#'
#'@return Returns median posterior estimates of the variance covariance matrix.
#'
#'@examples
#'\dontrun{
#'set.seed(999)
#'your_data=mvrnorm(n=15,mu=c(0,0),Sigma=matrix(c(4,3,3,9),nrow=2,ncol=2))
#'Mu0=c(0,0)
#'Sigma0=matrix(c(1,0.6,0.6,4),nrow=2,ncol=2)
#'Nu0=3-1
#'bcov(data=your_data,iter=5000,burn=2500,seed=999,CI=0.95,
#'     mu0=Mu0,S0=Sigma0,nu0=Nu0)}
#'
#'@export


bcov=function(data,iter,burn,seed,CI,S0,nu0,mu0){

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
VAR_M=NULL
VAR_SD=NULL
VAR_LL=NULL
VAR_UL=NULL
CI=ifelse(missing(CI),0.95,CI)
CI=ifelse(CI>1,CI/100,CI)
ll=(1-CI)/2
ul=1-ll
for (a in 1:ncol(SIGMA)){
VAR_M[a]=quantile(probs=c(0.5),SIGMA[burn:nrow(SIGMA),a])
VAR_SD[a]=sd(SIGMA[burn:nrow(SIGMA),a])
VAR_LL[a]=emp.hpd(SIGMA[burn:nrow(SIGMA),a],conf=CI)[1]
VAR_UL[a]=emp.hpd(SIGMA[burn:nrow(SIGMA),a],conf=CI)[2]

}


star_ll=ifelse(VAR_LL<0,1,0)
star_ul=ifelse(VAR_UL<0,1,0)
star=ifelse(star_ll+star_ul==1," ","*")
VAR_M1=paste(round(VAR_M,2),star)
COV=matrix(VAR_M1,nrow=ncol(data),ncol=ncol(data))
table=data.frame(COV)
colnames(table)=c(colnames(data))
rownames(table)=c(colnames(data))
Out=list()
Out$MU=THETA
Out$SIGMA=SIGMA
Out$VAR_M=matrix(VAR_M,nrow=ncol(data),ncol=ncol(data))
Out$VAR_SD=matrix(VAR_SD,nrow=ncol(data),ncol=ncol(data))
Out$VAR_LL=matrix(VAR_LL,nrow=ncol(data),ncol=ncol(data))
Out$VAR_UL=matrix(VAR_UL,nrow=ncol(data),ncol=ncol(data))
Out$table=table

return(Out)
}
