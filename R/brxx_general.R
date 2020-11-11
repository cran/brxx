#'brxx_general: Bayesian Estimation of Reliability from Variance Estimates
#'
#'This function estimates reliability from given true and error variance estimates.
#'
#'@param a True score variance estimate.
#'@param b Error variance estimate.
#'@param alpha Prior true score variance.
#'@param beta Prior error variance.
#'@param CI Credible interval quantile, as a decimal (ie, for 95 percent, 0.95).
#'
#'@return Returns estimated median and quantile based credible limits for reliability.
#'
#'@examples
#'a=18.7
#'b=3.3
#'alpha=3.51
#'beta=1.75
#'brxx_general(a=a,b=b,alpha=alpha,beta=beta,CI=0.95)
#'
#'@export



brxx_general=function(a,b,alpha,beta,CI){
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
  out=round(qbeta(p=c(ll,0.5,ul),alpha+a,beta+b),4)
  names(out)=c("LL","Median","UL")
  out
}
