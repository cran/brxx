#'brxx_Cor_general: Bayesian Estimation of Reliability from Correlation, General Form
#'
#'This function estimates reliability from correlation given the correlation estimate.
#'
#'@param cor Correlation estimate.
#'@param alpha Prior true score variance.
#'@param beta Prior error variance.
#'@param CI Credible interval quantile, as a decimal (ie, for 95 percent, 0.95).
#'@param items Number of test items.
#'
#'@return Returns estimated median and quantile based credible limits for reliability.
#'
#'@examples
#'brxx_Cor_general(cor=0.85,alpha=3.51,beta=1.75,CI=0.95,items=10)
#'
#'@export



brxx_Cor_general=function(cor,alpha,beta,CI,items){
  CI=ifelse(missing(CI),0.95,CI)
  CI=ifelse(CI>1,CI/100,CI)
  am=ifelse(missing(alpha),1,0)
  bm=ifelse(missing(beta),10,0)
  items=ifelse(missing(items),1,items)
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
  Out=round(qbeta(c(ll,0.5,ul),alpha+cor*items,beta+items*(1-cor)),4)
  names(Out)=c("LL","Median","UL")
  Out
  return(Out)

}
