#' Plot the results for the Bayesian models
#'  @title plot result for bafit
#' @param BAout \code{ba} object contains all the output of the Bayesian model
#' @param type
#'  \code{string} defines the type of the plot, which can be "pre","man","trace"
#'  \tabular{ll}{
#'    \tab \code{"pre"}  plot the predicted value against the true value\cr
#'    \tab \code{"trace"}  Traceplot \cr
#'  }
#' @param op \code{op} object created by create.options that is used for the analysis
#' @param iterStart the starting iteration for traceplot
#' @details You can either provide \code{baData} or y and Z for this function. If they're provided at the same time, \code{dataobj} will be used for the analysis
#' @examples 
#' baplot(dataobj=pig,BAout=ba,type="pre")
#' baplot(dataobj=pig,BAout=ba,type="trace",op=op)
#' @export
plot.ba <- function(BAout=NULL,type=c("pre","trace"),op=NULL,iterStart=NULL,col=c("black","red")){

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
			  plot(BAout$y,BAout$yhat,xlab="true values of the phenotype",ylab="predicted values",col=col[(!BAout$train%%2)+1],pch = 20)
			  legend("topleft",legend=c("Train","Val"),col = col,pch=20)
			}
			  else plot(BAout$y,BAout$yhat,xlab="true values of the phenotype",ylab="predicted values",col=col[1],pch = 20)
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



