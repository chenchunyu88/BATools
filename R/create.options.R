# create object for options of bayesian model
#'  @export
#'  @title create \code{\link{op}} object
create.options <- function(model=NULL,method=NULL,ante=NULL,poly=NULL,priors=NULL,init=NULL,D="V",update_para=NULL,
                          run_para=NULL,save.at=NULL,cv=NULL,print_mcmc=NULL,ncore=1,seed=1,convcrit=1E-4){
	#Testing code
	#model=NULL;method=NULL;ante=NULL;priors=NULL;init=NULL;update_para=NULL;
	#run_para=NULL;save.at=NULL;cv=NULL;print_mcmc=NULL;ncore=1;poly=NULL	
  	#define models
	models=c("GBLUP","rrBLUP","BayesA","BayesB","BayesC","IWBayesA")
  	if(!is.null(model))
  	{
		if(!(prod(model %in% models)))
		{
			tmp<-paste(models, collapse=", ")
			stop("Model must be one of ",tmp)
		}
  	}else model="GBLUP"

  	#define methods
	methods=c("MCMC","EM")
  	if(!is.null(method))
  	{
		if(!(prod(method %in% methods)))
		{
			tmp<-paste(methods, collapse=" or ")
			stop("Method must be ",tmp)
		}
  	} else method="MCMC"

  	if(!is.null(ante))
	{
		if(class(ante)!="logical") stop("Ante must be TRUE or FALSE")
  	}else ante=FALSE

	if(!is.null(poly))
	{
		if(class(poly)!="logical") stop("poly must be TRUE or FALSE")
		if(poly & is.null(Ainv) & method=="MCMC"){
			cat("Ainv matrix is required for MCMC to sample polygenenic effects, program to run without sampling it\n")
			poly=FALSE
		} 
		if(poly & is.null(A) & method=="EM"){
			cat("A matrix is required to include polygenenic effects in EM, program to run without it\n")
			poly=FALSE
		}
	}else poly=FALSE

	#define priors
	prior_names=c("nu_e","tau2_e","nu_s","tau2_s","shape_scale","rate_scale","cdef","cscalea","alphapi","betapi","df0","Sig0")
	prior_names_ante=c("nu_e","tau2_e","nu_s","tau2_s","shape_scale","rate_scale","cdef","cscalea","alphapi","betapi","mu_m_t","sigma2_m_t","df_var_t","scale_var_t")

	if(is.null(priors$nu_e)) priors$nu_e=-2
	if(is.null(priors$tau2_e)) priors$tau2_e=0
	
	#flat prior for scale in REML	
	if(method=="EM")
	{
		if(is.null(priors$nu_s)) priors$nu_s=-2
		if(is.null(priors$tau2_s)) priors$tau2_s=0	
	}

	
	if(model=="BayesA")
	{
		if(is.null(priors$shape_scale)) priors$shape_scale=0.5
		if(is.null(priors$rate_scale)) priors$rate_scale=0	
	}	
	if(model=="BayesB")
	{
		if(is.null(priors$shape_scale)) priors$shape_scale=0.1
		if(is.null(priors$rate_scale)) priors$rate_scale=0.1	
	}
	
	if(is.null(priors$cdef)) priors$cdef=0.5
	if(is.null(priors$cscalea)) priors$cscalea=0.5
	
	if(model%in%c("BayesC","BayesB"))
	{
		if(is.null(priors$alphapi)) priors$alphapi=1
		if(is.null(priors$betapi)) priors$betapi = 9
	}
	
	if(model=="IWBayesA")
	{
		if(is.null(priors$df0)) priors$df0 = 4
		if(is.null(priors$Sig0)) priors$Sig0 = matrix(c(0.5,0,0,0.5),2,2)	
	}	
	
	if(ante)
	{
		if(!(prod(names(priors) %in% prior_names_ante)))
		{
			tmp<-paste(prior_names_ante, collapse=", ")
			stop("Prior names must be ",tmp)
		}		
		if(is.null(priors$mu_m_t)) priors$mu_m_t=0
		if(is.null(priors$sigma2_m_t)) priors$sigma2_m_t=0.01
		if(is.null(priors$df_var_t)) priors$df_var_t=-1
		if(is.null(priors$scale_var_t)) priors$scale_var_t=0				
	}else{
		if(!(prod(names(priors) %in% prior_names)))
		{
			tmp<-paste(prior_names, collapse=", ")
			stop("Prior names must be ",tmp)
		}
	}	

        #define initial values
	initial_names=c("df","scale","pi","varu","d","Sig","vare","g","b","c")
	initial_names_ante=c("df","scale","pi","varu","mut","vart")	   

	if(is.null(init$df)) init$df=5
	if(is.null(init$scale)) init$scale=0.02
	if(is.null(init$pi)) 
	{
		if(model=="BayesB") init$pi=0.1
		else init$pi=1
	}
	if(is.null(init$pi)) init$varu=2.1


	if(ante)
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
	update_names=c("df","scale","pi")
	update_names_ante=c("df","scale","pi","mut","vart")
	
	
	if(is.null(update_para$df)) {update_para$df=TRUE}
	if(is.null(update_para$scale)) {update_para$scale=TRUE}
	if(is.null(update_para$pi)) 
	{
		if(model=="BayesB") update_para$pi=TRUE else update_para$pi=FALSE
	}
	
	

	if(ante)
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
				stop("Running parameter names for EM must be maxiter")
			}
			if(is.null(run_para$maxiter)) {run_para$maxiter=500}
		}	
	}else{
		if(method=="MCMC")
		{
			if(is.null(run_para$niter)) {run_para$niter=50000}
			if(is.null(run_para$burnIn)) {run_para$burnIn=10000}
			if(is.null(run_para$skip)) {run_para$skip=10}
		}else{
			if(is.null(run_para$maxiter)) {run_para$maxiter=500}
		}
	}
	
	#define save.at
	if(is.null(save.at)) save.at=""
	
	#define cv
	if(!is.null(cv)) 
	{
		if(!is.logical(cv)) stop("cv must be logical")
	}else
	{
		#add auto cv code here
	}
	
	#define print options
	print_names=c("piter","time_est","print_to")
	if(!is.null(print_mcmc))
	{
		if(!(prod(names(print_mcmc) %in% print_names)))
		{
			tmp<-paste(run_names_MCMC, collapse=", ")
			stop("Printing parameter names must be one of ",tmp)
		}
		if(is.null(print_mcmc$piter)) print_mcmc$piter=10000
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
	obj=list(model=model,method=method,ante=ante,poly=poly,priors=priors,init=init,D=D,update_para=update_para,
	                          run_para=run_para,save.at=save.at,cv=cv,print_mcmc=print_mcmc,ncore=ncore,seed=seed,convcrit=convcrit)
  	class(obj) <- "options"
  	return(obj)
}



