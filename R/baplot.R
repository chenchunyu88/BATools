#' Plot the results for the Bayesian models
#'  @title plot result for bafit
#'  @export
#'  @param dataobj \code{baData} object 
#'  @param BAout \code{BAout} object contains all the output of the Bayesian model
#'  @param y \code{numeric} vector of phenotypes
#'  @param Z \code{matrix} of genotypes
#'  @param yNa \code{numeric} vector of phenotypes used in cross validation
#'  @param type
#'  \code{string} defines the type of the plot, which can be "pre","man","trace"
#'  \tabular{ll}{
#'    \tab \code{"pre"}  plot the predicted value against the true value\cr
#'    \tab \code{"man"}  Manhattan plot  \cr
#'    \tab \code{"trace"}  Traceplot \cr
#'  }
#'  @param op \code{op} object created by create.options that is used for the analysis
#'  @details You can either provide \code{baData} or y and Z for this function. If they're provided at the same time, \code{dataobj} will be used for the analysis
#' @examples 
#' baplot(dataobj=pig,BAout=ba,type="pre")
#' baplot(dataobj=pig,BAout=ba,type="trace",op=op)
baplot <- function(dataobj=NULL,BAout=NULL,y=NULL,Z=NULL,yNa=NULL,type=NULL,op=NULL){
	if(is.null(BAout)) stop("BAout must be provided for the function")
	plot_types=c("pre","man","trace")
	if(type=="trace"){
	if(is.null(op)) sample_idx=1:dim(BAout$hypers)[1]
	else sample_idx=c((op$run_para$burnIn/op$run_para$skip+1):(op$run_para$niter/op$run_para$skip))
	}
	if(!is.null(type))
	{
		if(!(prod(names(type) %in% plot_types)))
		{
				tmp<-paste(plot_types, collapse=", ")
				stop("The type of the plot must be ",tmp)
		}
		
	}else{
		stop("please input the type of the plot (pre, man, or trace) you want")
	}
	if(!is.null(dataobj))
	{
		if(type=="pre")
		{
			plot(dataobj$pheno[BAout$idx,BAout$trait,],BAout$yhat,xlab="true values of the phenotype",ylab="predicted values")
			abline(a=0,b=1)
		}else if(type=="trace")
		{
			plot(as.mcmc(BAout$hypers[sample_idx,]))
		}else if(type=="man")
		{
			
		}
	}else if((!is.null(y)) && (!is.null(Z))){
		if(type=="pre")
		{
			plot(y,BAout$yhat,xlab="true values of the phenotype",ylab="predicted values")
			abline(a=0,b=1)
			if(!is.null(yNa))
			{
				whichNa=which(is.na(yNa))
				points(x=y[whichNa],y=BAout$yhat[whichNa],col=2,cex=.8,pch=19)
			}
		}else if(type=="trace")
		{
			plot(as.mcmc(BAout$hypers[sample_idx,]))
		}else if(type=="man")
		{
			
		}
	}else{
		stop("baData or y and Z must be provided for the function")
	}

}



