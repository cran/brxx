#'scree: Scree Plot with Pairwise Complete Cases
#'
#'This function provides a scree plot when data may be missing.
#'
#'@param data N by P data matrix.
#'
#'@return Returns eigenvalues and scree plot.
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
#'  for (p in 1:5){
#'    your_data_miss[i,p]=ifelse(runif(1,0,1)<0.2,NA,your_data[i,p])
#'  }
#'}
#'scree(your_data_miss)
#'
#'@export

scree=function(data){
  dat=data
  Value=round(eigen(cor(dat,use="pairwise.complete.obs"))$values,2)
  Pct=round(eigen(cor(dat,use="pairwise.complete.obs"))$values/ncol(dat)*100,2)
  Cum_Pct=cumsum(abs(Pct))
  table=cbind(round(Value,2),Pct,Cum_Pct)
  colnames(table)=c("Eigenvalue","% of variance","Cumulative %")
  rownames(table)=c(1:length(Value))
  Out=list()
  Out$Value=Value
  Out$Pct=Pct
  Out$Cum_Pct=Cum_Pct
  Out$Table=table
  plot(eigen(cor(dat,use="pairwise.complete.obs"))$values,
       main="Scree Plot",xlab="Factor",ylab="Eigenvalue")
  lines(eigen(cor(dat,use="pairwise.complete.obs"))$values)
  lines(rep(1,ncol(dat)),col="red",lty=2)


  return(Out)

}
