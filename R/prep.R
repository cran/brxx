#'prep: Prepare Data File for Bayesian Analysis
#'
#'This function prepares data for analysis using Stan factor analysis code.
#'
#'@param data N by P data matrix.
#'@param nfactors Number of factors to extract.
#'@param Prior Prior loading matrix.
#'
#'@return Returns a formatted data file for use with Stan MCMC sampler.
#'
#'@examples
#'set.seed(999)
#'your_data=data.frame(mvrnorm(n=20,mu=c(0,0,0,0,0),
#'                             Sigma=matrix(c(4,2,2,2,2,
#'                                            2,4,2,2,2,
#'                                            2,2,4,2,2,
#'                                            2,2,2,4,2,
#'                                            2,2,2,2,4),
#'                                          nrow=5, ncol=5)))
#'colnames(your_data)=c("x1","x2","x3","x4","x5")
#'your_data_miss=matrix(ncol=5,nrow=20)
#'for (i in 1:20){
#'for (p in 1:5){
#'    your_data_miss[i,p]=ifelse(runif(1,0,1)<0.2,NA,your_data[i,p])
#'  }
#'}
#'formatted_data=prep(your_data_miss,nfactors=3)
#'
#'
#'@export

prep=function(data,nfactors,Prior){
  dat2=data
  for(i in 1:nrow(dat2)){
    for(p in 1:ncol(dat2)){
      dat2[i,p]=ifelse(is.na(dat2[i,p]),0,dat2[i,p])
    }
  }
  N=nrow(dat2)
  P=ncol(dat2)
  Q=nfactors
  D_O=dat2
  Load_Prior=matrix(nrow=P,ncol=Q)
  for (p in 1:P){
    for (q in 1:Q){
      Load_Prior[p,q]=ifelse(missing(Prior),
                             (eigen(cor(dat2))$vectors*
                                sqrt(eigen(cor(dat2))$values))[p,q],
                             Prior[p,q])

    }
  }
  I=diag(P)
  R=matrix(nrow=N,ncol=P)
  for (i in 1:N){
    for (p in 1:P){
      R[i,p]=ifelse(is.na(data[i,p]),1,0)
    }
  }

  formatted_data=list(N=N,
                      P=P,
                      Q=Q,
                      D_O=D_O,
                      Load_Prior=Load_Prior,
                      I=I,
                      R=R)

  return(formatted_data)

}
