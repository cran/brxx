#'summarize: Summarize Stan output as median, SD, and HPD quantiles
#'
#'This function converts raw MCMC sample data into matrix formatted summaries
#'
#'@param Samples S by theta matrix of sampled parameter estimates.
#'@param nrow Number of rows of target summary matrix
#'@param ncol Number of columns of target summary matrix
#'@param CI Creddible interval quantile, as a decimal (ie, for 95 percent, 0.95)
#'
#'@import TeachingDemos
#'
#'@return Returns median, SD, and HPD CI limits
#'
#'@examples
#'\dontrun{
#'your_data_s=standardize(your_data)
#'formatted_data=prep(your_data_s,nfactors=3)
#'out=sampling(model, data=formatted_data, iter=5000, seed=999)
#'res=as.matrix(out)
#'unpacked=unpack(Samples=res,Format=formatted_data)
#'processed=process(Loading_Matrix=unpacked$Loading_Matrix,
#'                  Format=formatted_data,
#'                  Rotate="oblimin")
#'
#'summarize(processed$Loadings,
#'          nrow=Formatted_data$P,
#'          ncol=Formatted_data$Q)$Table
#'summarize(processed$Communality,
#'          nrow=Formatted_data$P,
#'          ncol=1)$Table
#'summarize(processed$Uniqueness,
#'          nrow=Formatted_data$P,
#'          ncol=1)$Table
#'summarize(processed$G_Factor,
#'          nrow=Formatted_data$P,
#'          ncol=1)$Table
#'summarize(processed$Interfactor_Correlations,
#'          nrow=Formatted_data$Q,
#'          ncol=Formatted_data$Q)$Table
#'summarize(processed$Omega,
#'          nrow=1,
#'          ncol=1)$Table
#'summarize(unpacked$Tau_Matrix,
#'          nrow=Formatted_data$P,
#'          ncol=1)$Table}
#'
#'@export


summarize=function(Samples,nrow,ncol,CI){
  Obj=Samples
  CI=ifelse(missing(CI),0.95,CI)
  CI=ifelse(CI>1,CI/100,CI)

  Med=NULL
  SD=NULL
  LL=NULL
  UL=NULL
  for (c in 1:ncol(Obj)){
    Med[c]=quantile(Obj[,c],0.5)
    SD[c]=sd(Obj[,c])
    LL[c]=emp.hpd(Obj[,c],conf=CI)[1]
    UL[c]=emp.hpd(Obj[,c],conf=CI)[2]
  }

  Rows=nrow
  Columns=ncol

  table_Median=matrix(round(Med,2),nrow=Rows,ncol=Columns)
  table_SD=matrix(round(SD,4),nrow=Rows,ncol=Columns)
  table_LL=matrix(round(LL,2),nrow=Rows,ncol=Columns)
  table_UL=matrix(round(UL,2),nrow=Rows,ncol=Columns)
  sig=ifelse(ifelse(LL<0,1,0)+ifelse(UL<0,1,0)==1," ","*")
  table_sig=matrix(sig,nrow=Rows,ncol=Columns)
  table=data.frame(matrix(noquote(paste(table_Median,
                                        " (",table_SD,")",
                                        table_sig,sep="")),
                          nrow=Rows,ncol=Columns))

  row.names(table)=colnames(Obj)
  colnames(table)=c(1:Columns)
  table

  Out=list()
  Out$Medians=table_Median
  Out$SDs=table_SD
  Out$LLs=table_LL
  Out$ULs=table_UL
  Out$Table=table
  return(Out)
}
