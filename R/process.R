#'process: rotates and calulates reliability for Stan output
#'
#'This function processes Stan loading matrix data.
#'
#'@param Loading_Matrix S by P*Q matrix of loading samples.
#'@param Format list formatted data file provided for Stan
#'@param Rotate If Q>1, rotation (for options, see GPArotation package)
#'
#'@import GPArotation
#'
#'@return Returns rotated loadings, uniqueness, communality, and reliability.
#'
#'@examples
#'\dontrun{
#'your_data_s=standardize(your_data)
#'formatted_data=prep(your_data_s,nfactors=3)
#'out=sampling(model, data=formatted_data, iter=5000, seed=999)
#'res=as.matrix(out)
#'unpacked=unpack(Samples=res,Format=formatted_data)
#'
#'processed=process(Loading_Matrix=unpacked$Loading_Matrix,
#'                  Format=formatted_data,
#'                  Rotate="oblimin")}
#'
#'
#'@export



process=function(Loading_Matrix,Format,Rotate){

  P=Format$P
  Q=Format$Q
  Loading_Matrix=Loading_Matrix
  Rotate=ifelse(missing(Rotate),"oblimin",Rotate)
  Loadings=NULL
  Interfactor_Correlations=NULL
  G_Factor=NULL
  Communality=NULL
  Uniqueness=NULL
  Omega=NULL
  if(Q==1){
    for (s in 1:ncol(Loading_Matrix)){
      raw_mat=matrix(Loading_Matrix[s,],nrow=P,ncol=Q)
      LL=raw_mat%*%t(raw_mat)
      Com=diag(LL)
      Uni=1-Com
      omeg=sum(raw_mat)^2/(sum(raw_mat)^2+sum(1-raw_mat^2))
      Loadings=rbind(Loadings,raw_mat)
      Communality=rbind(Communality,Com)
      Uniqueness=rbind(Uniqueness,Uni)
      Omega=rbind(Omega,omeg)

    }
    Out=list()
    Out$Loadings=Loadings
    Out$Communality=Communality
    Out$Uniqueness=Uniqueness
    Out$Omega=Omega
    }
  else{
  for (s in 1:ncol(Loading_Matrix)){
    raw_mat=matrix(Loading_Matrix[s,],nrow=P,ncol=Q)
    Rotated=GPFoblq(raw_mat,method=Rotate,normalize=TRUE)
    R_Load=c(Rotated$loadings)
    Int_Cor=c(Rotated$Phi)
    LL=raw_mat%*%t(raw_mat)
    G_Fact=c(eigen(LL)$vectors[,1]*sqrt(eigen(LL)$values[1]))
    Com=diag(LL)
    Uni=1-Com
    omeg=sum(G_Fact)^2/(sum(G_Fact)^2+sum(1-G_Fact^2))
    Loadings=rbind(Loadings,R_Load)
    Interfactor_Correlations=rbind(Interfactor_Correlations,Int_Cor)
    G_Factor=rbind(G_Factor,G_Fact)
    Communality=rbind(Communality,Com)
    Uniqueness=rbind(Uniqueness,Uni)
    Omega=rbind(Omega,omeg)

  }
    Out=list()
    Out$Loadings=Loadings
    Out$Interfactor_Correlations=Interfactor_Correlations
    Out$G_Factor=G_Factor
    Out$Communality=Communality
    Out$Uniqueness=Uniqueness
    Out$Omega=Omega

    }
return(Out)
}


