#'bcor: Bayesian Estimation of The Correlation Matrix
#'
#'This function estimates coefficient omega internal consistency reliability.
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
#'@return Returns median posterior estimates of the correlation matrix.
#'
#'@import MCMCpack
#'
#'@examples
#'set.seed(999)
#'your_data=mvrnorm(n=15,mu=c(0,0),Sigma=matrix(c(4,3,3,9),nrow=2,ncol=2))
#'Mu0=c(0,0)
#'Sigma0=matrix(c(1,0.6,0.6,4),nrow=2,ncol=2)
#'Nu0=1
#'bcor(data=your_data,iter=5000,burn=2500,seed=999,CI=0.95,
#'     mu0=Mu0,S0=Sigma0,nu0=Nu0)
#'
#'
#'@export


bcor=function(data,iter,burn,seed,CI,S0,nu0,mu0){

  filler=matrix(nrow=ncol(data),ncol=ncol(data))
  for (a in 1:ncol(data)){
    for (b in 1:ncol(data)){
      filler[a,b]=ifelse(a==b,cov(data,use="pairwise.complete.obs")[a,b],0)}}
  filler1=matrix(nrow=ncol(data),ncol=ncol(data))
  for (a in 1:ncol(data)){
    for (b in 1:ncol(data)){
      filler1[a,b]=ifelse(missing(S0),filler[a,b],S0[a,b])}}
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

  set.seed(seed)
  pct=rep(0,iter+1)
  print(noquote("Sampling, this may take a minute"))
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
    pct[s+1]=(round(s/iter*10))*10
    if(pct[s+1]!=pct[s]){print(noquote(paste(pct[s+1],"%")))}

  }
  VAR_M=NULL
  VAR_LL=NULL
  VAR_UL=NULL
  CI=ifelse(missing(CI),0.95,CI)
  CI=ifelse(CI>1,CI/100,CI)
  ll=(1-CI)/2
  ul=1-ll
  for (a in 1:ncol(SIGMA)){
    VAR_M[a]=quantile(probs=c(0.5),SIGMA[burn:nrow(SIGMA),a])
    VAR_LL[a]=quantile(probs=c(ll),SIGMA[burn:nrow(SIGMA),a])
    VAR_UL[a]=quantile(probs=c(ul),SIGMA[burn:nrow(SIGMA),a])

  }



  COR=NULL
  mat=matrix(nrow=ncol(data),ncol=ncol(data))
  cor=matrix(nrow=ncol(data),ncol=ncol(data),0)
  print(noquote("Calculating correlations, this may take a minute"))
  pct=rep(0,nrow(SIGMA)-burn+1)
  for (s in burn:nrow(SIGMA)){
    mat=matrix(SIGMA[s,],nrow=ncol(data),ncol=ncol(data))
    for (a in 1:ncol(data)){
      for (b in 1:ncol(data)){
        cor[a,b]=mat[a,b]/sqrt(mat[a,a]*mat[b,b])
        COR=rbind(COR,c(cor))
      }
    }
    num=(s-burn+1)
    denom=(nrow(SIGMA)-burn)
    pct[s-burn+2]=round((num/denom)*10)*10
    if(pct[s-burn+2]!=pct[s-burn+1]){print(noquote(paste(pct[s-burn+2],"%")))}

  }
  COR_M=NULL
  COR_LL=NULL
  COR_UL=NULL
  for (a in 1:ncol(COR)){
    COR_M[a]=quantile(probs=c(0.5),COR[1:nrow(COR),a])
    COR_LL[a]=quantile(probs=c(ll),COR[1:nrow(COR),a])
    COR_UL[a]=quantile(probs=c(ul),COR[1:nrow(COR),a])

  }

  star_ll=ifelse(COR_LL<0,1,0)
  star_ul=ifelse(COR_UL<0,1,0)
  star=ifelse(star_ll+star_ul==1," ","*")
  COR_M1=paste(round(COR_M,2),star)
  COR1=matrix(COR_M1,nrow=ncol(data),ncol=ncol(data))
  table=data.frame(COR1)
  colnames(table)=c(colnames(data))
  rownames(table)=c(colnames(data))
  diag(table)="1 "
  Out=list()
  Out$table=table
  Out$M=COR_M
  Out$LL=COR_LL
  Out$UL=COR_UL
  Out$THETA=THETA
  Out$SIGMA=SIGMA
  Out$COR=COR

  print(noquote(""))
  print(noquote(""))
  print(noquote("Correlation Matrix"))
  print(noquote("*Credible interval excludes 0"))
  Out$table
  return(Out)
}
