#'unpack: Unpack Stan output for factor analysis samples from Stan
#'
#'This function unpacks raw Stan samples output.
#'
#'@param Samples S by theta matrix of sample parameter estimates.
#'@param Format list formatted data file provided for Stan
#'
#'@return Returns four matrices:
#'@return 1). S by Q latent score matrix, x.
#'@return 2). S by Q*P loading matrix, lambda.
#'@return 3). S by P mean matrix, tau.
#'@return 4). S by P loading variance matrix, alpha.
#'
#'@examples
#'\dontrun{
#'your_data_s=standardize(your_data)
#'formatted_data=prep(your_data_s,nfactors=3)
#'out=sampling(model, data=formatted_data, iter=5000, seed=999)
#'res=as.matrix(out)
#'
#'
#'unpacked=unpack(Samples=res,Format=formatted_data)}
#'
#'
#'@export


unpack=function(Samples,Format){
  N=Format$N
  P=Format$P
  Q=Format$Q
  Out=list()
  Out$X_Matrix=Samples[,1:(Q*N)]
  Out$Loading_Matrix=Samples[,(Q*N+1):(Q*N+P*Q)]
  Out$Tau_Matrix=Samples[,(Q*N+P*Q+1):(Q*N+P*Q+P)]
  Out$Alpha_Matrix=Samples[,(Q*N+P*Q+P+1):(Q*N+P*Q+P+Q)]
  return(Out)

}

