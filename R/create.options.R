#' Create object for options of bayesian model
#' @export
#' @title create \code{\link{options}} object
#' @param  model string indicate the model for the analysis, \code{model} can be "GBLUP","rrBLUP,"BayesA", "BayesB","SSVS","ssBayesA", "ssBayesB", or "ssSSVS","anteBayesA" and "anteBayesB"
#' @param method string indicate the method for the analysis, \code{model} can be "MCMC", "MAP" or "REML"
#' @param ssGBLUPvar string indicate the variance component treatment in ssGBLUP \code{ssGBLUPvar} can be c("homVAR", "hetVAR"), default is "homVAR"
#' @param priors \code{list} contains priors for the Bayesian model, elements in \code{priors} can be "nu_e","tau2_e",\cr"shape_scale","rate_scale","cdef","alphapi","betapi","mu_m_t","sigma2_m_t",\cr
#'  "df_var_t","scale_var_t"
#'  
#' \code{nu_e}    \code{numeric}, prior degrees of freedom for residual variance, default value -1\cr
#' \code{tau2_e}   \code{numeric},  prior scale for residual variance,default value 0  \cr
#' \code{shape_scale}  \code{numeric},  prior shape for scale parameter, default value 0.1 \cr
#' \code{rate_scale}    \code{numeric},  prior rate for scale parameter, default value 0.1 \cr 
#' \code{cdef}   \code{numeric},  cdef is scale parameter for proposal density on df, default value 0.5\cr
#' \code{alphapi}   \code{numeric}, prior \eqn{\alpha} for \code{pi}, default value 1 \cr   
#' \code{betapi}    \code{numeric}, prior \eqn{\beta} for \code{pi}, default value 9\cr 
#' \code{mu_m_t}    \code{numeric}, prior mean for \eqn{\mu_t}, default value 0\cr 
#' \code{sigma2_m_t}   \code{numeric}, prior variance for \eqn{\mu_t}, default value 0.01 \cr  
#' \code{df_var_t}  \code{numeric}, prior degrees of freedom for \eqn{\sigma^2_t}, default value -1\cr    
#' \code{scale_var_t}     \code{numeric}, prior scale for \eqn{\sigma^2_t}, default value 0, the scale and degrees of freedom for \eqn{\sigma^2_t} are non-informative priors based on the original paper\cr 
  
#' @param init \code{list} contains initial values for the Bayesian model, elements in \code{init} can be "df","scale","pi","mut","vart"\cr
  
#' \code{df}    \code{numeric}, the starting value of degrees of freedom parameter for SNP effect variance, default value 5\cr
#' \code{scale}   \code{numeric},  the starting value of scale parameter for SNP effect variance, default value 0.02  \cr
#' \code{pi}  \code{numeric}, the starting value of \eqn{\pi}, which is the percentage of SNP that has variantion to the phenotype, default value 0.1 for BayesB, 1 for other models \cr
#' \code{mut}    \code{numeric},  the starting value of \eqn{\mu_t}, default value 0 \cr 
#' \code{vart}   \code{numeric},  the starting value of \eqn{\sigma^2_t}, default value 0.5\cr 
#' @param D \code{string} indicate use relative variances ("V") or relative precisions ("P") in MAP BayesA, default is "P" 
#' @param update_para \code{list} of \code{logical} indicate whether a parameter is sampled in the Bayesian model, elements in \code{update_para} can be "df","scale","pi","mut","vart"
  
#' \code{df}    \code{logical}, \code{TRUE} if degrees of freedom parameter for SNP effect variance needs to be sampled\cr
#' \code{scale}   \code{logical},   \code{TRUE} if scale parameter for SNP effect variance  needs to be sampled\cr
#' \code{pi}  \code{logical}, \code{TRUE} if \eqn{\pi} needs to be sampled\cr
#' \code{mut}    \code{logical}, \code{TRUE} if  \eqn{\mu_t} needs to be sampled\cr 
#' \code{vart}   \code{logical}, \code{TRUE} if \eqn{\sigma^2_t} needs to be sampled\cr 
#' \code{vare}   \code{logical}, \code{TRUE} if \eqn{\sigma^2_e} needs to be sampled\cr 
#' @param run_para
#' \code{list} elements in \code{run_para} can be "niter","burnIn","skip"
#' \code{niter}     \code{numeric}, the number of iterates for MCMC sampling \cr
#' \code{burnIn}    \code{numeric},   burnIn period for MCMC sampling \cr
#' \code{skip}   \code{numeric} , skip for MCMC sampling\cr
#' @param save.at \code{string} define the directory to save the results
#' @param print_mcmc \code{list} define the monitor options for the Bayesian model, elements in \code{print_mcmc} can be "piter","time_est","print_to"
#' \tabular{ll}{
#'    \tab \code{piter}    numeric, print status every piter. If piter=0, the program don't print any status\cr
#'    \tab \code{time_est}   logical, \code{TRUE} mean display time left estimates \cr
#'    \tab \code{print_to}  string, "disk" means print to disk in file "log.txt"; "screen" means print to screen \cr
#'  }
#' @param ncore  \code{numeric}, the number of cpu cores for the analysis
#' @param seed \code{numeric}, seed for random number generator
#' @param convcrit \code{numeric}, convergence criteria for EM, default is 1E-4
#' @return a object of \code{options} for BATools, it's basicaly a list
#' @examples \dontrun{
#'   op=create.options(model="SSVS") #This will use the default settings
#'
#'  #For set your values, check out one of the demos, for example
#'  demo(SSVS)
#'  }
set.options <- function(model=c("GBLUP","rrBLUP","BayesA","BayesB","SSVS","ssGBLUP","ssBayesA","ssBayesB","ssSSVS",
                                   "anteBayesA","anteBayesB","anteSSVS","IWBayesA"),method=c("MCMC","MAP","REML"),ssGBLUPvar=c("homVAR","hetVAR"),
                                    priors=NULL,init=NULL,D="P",update_para=NULL,
                                    run_para=NULL,save.at=NULL,print_mcmc=NULL,ncore=1,seed=1,convcrit=1E-4){

    model<-match.arg(model)
    method<-match.arg(method)
    ssGBLUPvar<-match.arg(ssGBLUPvar)
   
  	
  	#define methods
  	if(model=="GBLUP") method="REML"
	  #methods=c("MCMC","MAP","REML")
  	if(is.null(method)) method="MCMC"

  
	#define priors
	prior_names=c("nu_e","tau2_e","nu_s","tau2_s","shape_scale","rate_scale","cdef","cscalea","alphapi","betapi","df0","Sig0","def_sigmag0","scale_sigmag0")
	prior_names_ante=c("nu_e","tau2_e","nu_s","tau2_s","shape_scale","rate_scale","cdef","cscalea","alphapi","betapi","mu_m_t","sigma2_m_t","df_var_t","scale_var_t")

	if(is.null(priors$nu_e)) priors$nu_e=-2
	if(is.null(priors$tau2_e)) priors$tau2_e=0
	
	#flat prior for scale in REML	
	if(method %in% c("MAP","REML"))
	{
		if(is.null(priors$nu_s)) priors$nu_s=-2
		if(is.null(priors$tau2_s)) priors$tau2_s=0	
	}

	
	if(model %in% c("BayesA","ssBayesA"))
	{
		if(is.null(priors$shape_scale)) priors$shape_scale=0.5
		if(is.null(priors$rate_scale)) priors$rate_scale=0	
	}	
	if(model %in% c("BayesB","ssBayesB"))
	{
		if(is.null(priors$shape_scale)) priors$shape_scale=0.1
		if(is.null(priors$rate_scale)) priors$rate_scale=0.1	
	}
	
	if(is.null(priors$cdef)) priors$cdef=0.5
	if(is.null(priors$cscalea)) priors$cscalea=0.5
	
	if(model%in%c("SSVS","BayesB","ssSSVS","ssBayesB"))
	{
		if(is.null(priors$alphapi)) priors$alphapi=1
		if(is.null(priors$betapi)) priors$betapi = 9
	}
	
	if(model=="IWBayesA")
	{
		if(is.null(priors$df0)) priors$df0 = 4
		if(is.null(priors$Sig0)) priors$Sig0 = matrix(c(0.5,0,0,0.5),2,2)	
	}	
	if(substr(model,1,2)=="ss"){
	  if(is.null(priors$def_sigmag0)) priors$def_sigmag0=-1
	  if(is.null(priors$scale_sigmag0)) priors$scale_sigmag0=0
	}
	
	if(substr(model,1,4)=="ante")
	{
		if(!(prod(names(priors) %in% prior_names_ante)))
		{
			tmp<-paste(prior_names_ante, collapse=", ")
			stop("Prior names must be ",tmp)
		}	
	  if(is.null(priors$nu_e)) priors$nu_e=-1
	  if(is.null(priors$tau2_e)) priors$tau2_e=0
	  if(is.null(priors$shape_scale)) priors$shape_scale=0.1
	  if(is.null(priors$rate_scale)) priors$rate_scale=0.1
		if(is.null(priors$mu_m_t)) priors$mu_m_t=0
		if(is.null(priors$sigma2_m_t)) priors$sigma2_m_t=0.01
		if(is.null(priors$df_var_t)) priors$df_var_t=-1
		if(is.null(priors$scale_var_t)) priors$scale_var_t=0	
		if(is.null(priors$alphapi)) priors$alphapi=1
		if(is.null(priors$betapi)) priors$betapi = 9
	}else{
		if(!(prod(names(priors) %in% prior_names)))
		{
			tmp<-paste(prior_names, collapse=", ")
			stop("Prior names must be ",tmp)
		}
	}	

        #define initial values
	initial_names=c("df","scale","pi","varu","d","Sig","vare","g","beta","c","post_prob","varGenetic")
	initial_names_ante=c("df","scale","pi","varu","mut","vart","vare","g","beta","c","post_prob","varGenetic")  

	if(is.null(init$df)) init$df=5
	if(is.null(init$scale)) init$scale=0.02
	if(is.null(init$varGenetic)) init$varGenetic=init$scale
	if(is.null(init$pi)) 
	{
		if(model=="BayesB") init$pi=0.1
		else init$pi=1
	}
	if(is.null(init$pi)) init$varu=2.1


	if(substr(model,1,4)=="ante")
	{
		if(is.null(init$mut)) init$mut=0
		if(is.null(init$vart)) init$vart=0.5
		if(!(prod(names(init) %in% initial_names_ante)))
		{
			tmp<-paste(initial_names_ante, collapse=", ")
			stop("Initial names must be ",tmp)
		}		
	}else
	{
		if(!(prod(names(init) %in% initial_names)))
		{
			tmp<-paste(initial_names, collapse=", ")
			stop("Initial names must be ",tmp)
		}
	}
	
	if(model=="BayesB" & init$pi==1) 
	{
		stop("starting value of pi for BayesB must be less than 1\n")
	}
	
	#define update parameters
	update_names=c("df","scale","pi","vare")
	update_names_ante=c("df","scale","pi","mut","vart","vare")
	
	
	if(is.null(update_para$df)) {update_para$df=FALSE}
	if(is.null(update_para$scale)) {update_para$scale=TRUE}
	if(is.null(update_para$vare)) {update_para$vare=TRUE}
	if(is.null(update_para$pi)) 
	{
		if(model=="BayesB") update_para$pi=TRUE else update_para$pi=FALSE
	}
	
	

	if(substr(model,1,4)=="ante")
	{
		if(is.null(update_para$mut)) {update_para$mut=TRUE}
		if(is.null(update_para$vart)) {update_para$vart=TRUE}
		if(!(prod(names(update_para) %in% update_names_ante)))
		{
			tmp<-paste(update_names_ante, collapse=", ")
			stop("Update parameter names must be ",tmp)
		}
	}else
	{
		if(!(prod(names(update_para) %in% update_names)))
		{
			tmp<-paste(update_names, collapse=", ")
			stop("Update parameter names must be ",tmp)
		}
	}
		
	#define running parameters
	run_names_MCMC=c("niter","burnIn","skip")
	if(!is.null(run_para))
	{
		if(method=="MCMC")
		{
			if(!(prod(names(run_para) %in% run_names_MCMC)))
			{
				tmp<-paste(run_names_MCMC, collapse=", ")
				stop("Running parameter names for MCMC must be ",tmp)
			}
			if(is.null(run_para$niter)) {run_para$niter=50000}
			if(is.null(run_para$burnIn)) {run_para$burnIn=10000}
			if(is.null(run_para$skip)) {run_para$skip=10}
			if(run_para$niter<2000)
			{
				stop("niter must be larger than 2000 to get good estimates")
			}
			if(run_para$niter<=run_para$burnIn)
			{
				stop("niter must be larger than burnIn for MCMC")
			}
		}else{
			if(!(names(run_para)=="maxiter"))
			{
				stop("Running parameter names for MAP must be maxiter")
			}
			if(is.null(run_para$maxiter)) {run_para$maxiter=50}
		}	
	}else{
		if(method=="MCMC")
		{
			if(is.null(run_para$niter)) {run_para$niter=50000}
			if(is.null(run_para$burnIn)) {run_para$burnIn=30000}
			if(is.null(run_para$skip)) {run_para$skip=10}
		}else{
			if(is.null(run_para$maxiter)) {run_para$maxiter=50}
		}
	}
	
	#define save.at
	if(is.null(save.at)) save.at=model
	
	
	#define print options
	print_names=c("piter","time_est","print_to")
	if(!is.null(print_mcmc))
	{
		if(!(prod(names(print_mcmc) %in% print_names)))
		{
			tmp<-paste(run_names_MCMC, collapse=", ")
			stop("Printing parameter names must be one of ",tmp)
		}
		if(is.null(print_mcmc$piter)) print_mcmc$piter=500
		if(is.null(print_mcmc$time_est)) print_mcmc$time_est=T
		if(is.null(print_mcmc$print_to)) print_mcmc$print_to="screen"
		if(!is.numeric(print_mcmc$piter)) stop("piter must be numeric")
		if(!is.logical(print_mcmc$time_est)) stop("time_est must be logical")
		print_tos=c("screen","disk")
		if(!(prod(print_mcmc$print_to %in% print_tos))) stop("print_to must be screen or disk")	
	}else{
		print_mcmc$piter=10000
		print_mcmc$time_est=T
		print_mcmc$print_to="screen"	
	}
	
	#define ncore
	if(!is.null(ncore))
	{
		if(!is.numeric(ncore)) stop("ncore must be numeric")
		
	}else{
		ncore=1
	}
  	# return object of class 'options'
	obj=list(model=model,method=method,priors=priors,init=init,D=D,ssGBLUPvar=ssGBLUPvar,update_para=update_para,
	                          run_para=run_para,save.at=save.at,print_mcmc=print_mcmc,ncore=ncore,seed=seed,convcrit=convcrit)
  	class(obj) <- "options"
  	return(obj)
}

