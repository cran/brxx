#'bomega: Bayesian Estimation of Coefficient Omega
#'
#'This function estimates coefficient omega internal consistency reliability.
#'
#'@param K The number of test items.
#'@param mod A measurement model estimated as a bsem object by blavaan.
#'@param alpha Prior true score variance.
#'@param beta Prior error variance.
#'@param CI Credible interval quantile, as a decimal (ie, for 95 percent, 0.95).
#'
#'@return Returns estimated median and quantile based credible limits for omega.
#'
#'@examples
#'
#'\donttest{
#'your_data=data.frame(mvrnorm(n=20,mu=c(0,0,0,0,0),
#'Sigma=matrix(c(4,2,2,2,2,
#'               2,4,2,2,2,
#'               2,2,4,2,2,
#'               2,2,2,4,2,
#'               2,2,2,2,4),
#'             nrow=5, ncol=5)))
#'colnames(your_data)=c("x1","x2","x3","x4","x5")
#'mod='tau=~x1+x2+x3+x4+x5'
#'fit=bsem(mod,data=your_data)
#'bomega(K=5,mod=fit,alpha=3.51,beta=1.75,CI=0.95)}
#'
#'@export
bomega=function(K,mod,alpha,beta,CI){

  s=unlist(summary(mod))
  out=NULL
  hold=NULL
  for (k in 1:K){

    hold=s[(length(s)-(3*K-k+2))]
    out=c(out,hold)
  }
  Lambda2=sum(as.numeric(out))^2
  out=NULL
  hold=NULL
  for (k in 1:K){

    hold=s[(length(s)-(3*K-k+2)+K)]
    out=c(out,hold)
  }
  Psi=sum(as.numeric(out))
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
  ll=(1-CI)/2
  ul=1-ll
  out=round(qbeta(c(ll,0.5,ul),alpha+Lambda2,beta+Psi),4)
  names(out)=c("LL","Median","UL")
  out
  return(out)
}
