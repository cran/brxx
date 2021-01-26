#'bomega_general: Bayesian Estimation of Coefficient Omega, General Form
#'
#'This function estimates coefficient omega internal consistency reliability.
#'
#'@param lambda vector of item loadings.
#'@param psi vector of item variances.
#'@param alpha Prior true score variance.
#'@param beta Prior error variance.
#'@param CI Credible interval quantile, as a decimal (ie, for 95 percent, 0.95).
#'
#'@import TeachingDemos
#'
#'@return Returns estimated median and quantile based credible limits for omega.
#'
#'@examples
#'lambda=c(0.7,0.5,0.6,0.7)
#'psi=c(0.2,0.4,0.3)
#'alpha=3.51
#'beta=1.75
#'bomega_general(lambda=lambda,psi=psi,alpha=alpha,beta=beta,CI=0.95)
#'
#'@export


bomega_general=function(lambda,psi,alpha,beta,CI){
  Lambda2=sum(lambda)^2
  Psi=sum(psi)
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
  Out=round(c(qbeta(c(0.5),alpha+Lambda2,beta+Psi),
              hpd(qbeta,shape1=alpha+Lambda2,shape2=beta+Psi,conf=CI)[1],
              hpd(qbeta,shape1=alpha+Lambda2,shape2=beta+Psi,conf=CI)[2]),4)
  names(Out)=c("Median","LL","UL")
  Out
  return(Out)
}
