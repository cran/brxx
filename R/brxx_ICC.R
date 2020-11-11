#'brxx_ICC: Bayesian Estimation of Reliability from ICC
#'
#'This function estimates reliability from intraclass correlation coefficient
#'
#'@param mod A mixed effects model object estimated by blmer.
#'@param alpha Prior true score variance (subject variance)
#'@param beta Prior error variance (residual variance)
#'@param CI Credible interval quantile, as a decimal (ie, for 95 percent, 0.95).
#'
#'@return Returns estimated median and quantile based credible limits for ICC.
#'
#'@examples
#'\donttest{
#'your_data_wide=mvrnorm(20,c(0,0),matrix(c(1,0.8,0.8,1),nrow=2,ncol=2))
#'your_data_long=c(as.vector(your_data_wide[,1]),as.vector(your_data_wide[,2]))
#'time=c(rep(0,20),rep(1,20))
#'id=c(rep(1:20,2))
#'mod=blmer(your_data_long~time+(1|id))
#'brxx_ICC(mod=mod,alpha=3.51,beta=1.75,CI=0.95)}
#'
#'@export



brxx_ICC=function(mod,alpha,beta,CI){
  s=summary(mod)
  ws=as.numeric(s$varcor[1])
  er=s$sigma
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
  out=round(qbeta(c(ll,0.5,ul),alpha+ws,beta+er),4)
  names(out)=c("LL","Median","UL")
  out
}
