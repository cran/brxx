#'standardize: Standardization of Data Matrix
#'
#'This function standardizes an N by P data matrix, as is strongly recommended before
#'using any of the brxx reliability estimation functions
#'
#'@param data N by P data matrix.
#'
#'@return Returns an item level standardized data matrix.
#'
#'@examples
#'set.seed(999)
#'your_data=data.frame(mvrnorm(n=20,mu=c(0,0,0,0,0,0,0),
#'                             Sigma=matrix(c(4,2,2,2,2,2,2,
#'                                            2,4,2,2,2,2,2,
#'                                            2,2,4,2,2,2,2,
#'                                            2,2,2,4,2,2,2,
#'                                            2,2,2,2,4,2,2,
#'                                            2,2,2,2,2,4,2,
#'                                            2,2,2,2,2,2,4),
#'                                          nrow=7, ncol=7)))
#'your_data_miss=matrix(ncol=5,nrow=20)
#'for (i in 1:20){
#'for (p in 1:5){
#'    your_data_miss[i,p]=ifelse(runif(1,0,1)<0.2,NA,your_data[i,p])
#'  }
#'}
#'standardize(your_data_miss)
#'
#'@export


standardize=function(data){
  Out=data
  for(p in 1:ncol(data)){
    for(i in 1:nrow(data)){
      Out[i,p]=ifelse(is.na(Out[i,p]),NA,
                      (data[i,p]-mean(data[,p],na.rm=T))/(sd(data[,p],na.rm=T)))
    }
  }
  return(Out)
}

