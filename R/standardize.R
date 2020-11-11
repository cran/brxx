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
#'standardize(your_data)
#'
#'@export


standardize=function(data){out=matrix(nrow=nrow(data),ncol=ncol(data))
for(p in 1:ncol(data)){
  for(i in 1:nrow(data)){
    out[i,p]=(data[i,p]-mean(data[,p],na.rm=T))/(sd(data[,p],na.rm=T))
  }
}
out
}

