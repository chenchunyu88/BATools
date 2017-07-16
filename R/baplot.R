#' Plot the results for the Bayesian models
#' @title plot result for baFit
#' @param BAout \code{ba} object contains all the output of the Bayesian model
#' @param type
#'  \code{string} defines the type of the plot, which can be "pre","trace"
#'  \tabular{ll}{
#'    \tab \code{"pre"}  plot the predicted value against the true value\cr
#'    \tab \code{"trace"}  Traceplot \cr
#'  }
#' @param op \code{op} object created by `create.options` that is used for the analysis
#' @param iterStart the starting iteration for traceplot
#' @param col color for the plot in cross-validation
#' @examples \dontrun{
#'  rm(list=ls())
#'  library(BATools)
#'  data("Pig")
#'  #Standardize genotype matrix
#   geno=std_geno(PigM,method="s",freq=PigAlleleFreq)
#   init=set_init("driploss",data=PigPheno,geno=geno,"id",df=5,pi_snp=0.001,h2=0.5,c=1000,model="SSVS",centered=TRUE)
#'  #or set your own starting values using 
#'  #init=list(df=5,scale=0.01,pi=1) 
#'  run_para=list(niter=2000,burnIn=1000,skip=10)
#'  print_mcmc=list(piter=500)
#'  update_para=list(df=FALSE,scale=TRUE,pi=F)
#'  op<-create.options(model="SSVS",method="MCMC",seed=1,priors=NULL,init=init,
#'                   update_para=update_para,run_para=run_para,save.at="SSVS",print_mcmc=print_mcmc)
#'  #### Cross-validation using BATools
#'  set.seed(1234)
#'  PigPheno=createCV(data = PigPheno,k=5,"driploss")
#'  cvSSVS<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1)
#'  par(mfrow=c(1,1))
#'  plot(cvSSVS)
#' }
#' @export
baplot <- function(BAout=NULL,type=c("pre","trace"),op=NULL,iterStart=NULL,col=c("black","red"),...){

  if(is.null(BAout)) stop("BAout must be provided for the function")
	type=match.arg(type)
	
	
	if(type=="trace"){
	  if(is.null(op)) stop("options need to known for traceplot")
	  if(is.null(iterStart)) sample_idx=c((op$run_para$burnIn %/% op$run_para$skip+1):(op$run_para$niter/op$run_para$skip))
	  else sample_idx=c((iterStart %/% op$run_para$skip+1):(op$run_para$niter/op$run_para$skip))
	}

		

	
		if(type=="pre")
		{
			if(!is.null(BAout$train)) {
			  plot(BAout$y,BAout$yhat,xlab="true values of the phenotype",ylab="predicted values",col=col[(!BAout$train%%2)+1],pch = 20,...)
			  legend("topleft",legend=c("Train","Val"),col = col,pch=20)
			}
			  else plot(BAout$y,BAout$yhat,xlab="true values of the phenotype",ylab="predicted values",col=col[1],pch = 20,...)
		  abline(a=0,b=1)
			#if(!is.null(yNa))
			#{
			#	whichNa=which(is.na(yNa))
			#	points(x=y[whichNa],y=BAout$yhat[whichNa],col=2,cex=.8,pch=19)
			#}
		}else if(type=="trace")
		{
			plot(as.mcmc(BAout$hypers[sample_idx,]))
		}
}



