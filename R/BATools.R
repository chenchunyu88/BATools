#' bafit function can fit various Bayesian models including rrBLUP, BayesA/B/C, antedependence models and etc. Input data can be either a `baData` object or specify the fixed and random effects seperately.
#' @title Fitting various Bayesian models
#' @param dataobj A list of baData including phenotypes, genotypes and etc.
#' @param op A list of options to run Bayesian models
#' @param y A numeric vector of phenotypes
#' @param A matrix of genotypes
#' @param trait A string indicating the trait for analysis
#' @return The result of the analysis as a ba object, a list of estimate of fixed and random effects as well as variance components.
#' @examples
#' ###########Loading and preparing data##########
#' library(BATools)
#' data("MSUPRP_sample")
#' summary(MSUPRP_sample)
#' pheno<-data.frame(MSUPRP_sample$pheno[,,]) 
#' geno<-MSUPRP_sample$geno[,1:500]
#' ped<-MSUPRP_sample$pedigree
#' map=MSUPRP_sample$map
#' sex<-ped$sex
#' sex<-as.factor(sex)
#' x<-model.matrix( ~ sex -1,contrasts.arg=list(sex=contrasts(sex, contrasts=F)))
#' colnames(x)<-c("female","male")
#' rownames(x)<-ped$ID
#' pig=create.baData(pheno=pheno,geno=geno,map=map,pedigree=ped,fixed=x,makeAinv=F)
#' ######################BayesA####################
#' ##############Setting up options################
#' init=list(df=5,scale=0.01,pi=1)
#' run_para=list(niter=2000,burnIn=1000,skip=10)
#' print_mcmc=list(piter=200)
#' update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
#' op<-create.options(model="BayesA",method="MCMC",ante=FALSE,priors=NULL,init=init,
#'                    update_para=update_para,run_para=run_para,save.at="BayesA",cv=NULL,print_mcmc=print_mcmc)
#' #################run BayesA MCMC################
#' ba<-bafit(dataobj=pig,op=op,trait="driploss")
#' ba
#'  @export
#'  @useDynLib BATools
bafit <- function(dataobj=NULL,op=NULL,y=NULL,Z=NULL,X=NULL,trait=NULL){
	if(is.null(dataobj)) 
	{
		if(is.null(y) || is.null(Z)) stop("must have baData or y and Z to start the function")
		if(!is.null(Z) && class(Z)!="matrix") {
			warning("Z is not a matrix, converting it to a matrix, this may cause issue in the estimates")
			Z=as.matrix(Z)
		}
		if(!is.null(X) && class(X)!="matrix") {
			warning("Z is not a matrix, converting it to a matrix, this may cause issue in the estimates")
			X=as.matrix(X)
		}
	}else{
	 	if(class(dataobj)!="baData") stop("dataobj must be baData")
	}
	
	if(is.null(op)) cat("no options inputed, use the default options.\n")
	
	if(op$ante)
	{
		if(is.null(dataobj$map))
		{
			cat("no map is provided for antedepedence model, program will assume the markers are in the same chromosome.\n")
			ante_p=NULL
		}
		else{
			if(max(dataobj$map$chr)>1)
			{
				map=dataobj$map[which(rownames(dataobj$map)%in%colnames(dataobj$geno)),]
				ante_p=rep(0,max(map$chr))
				ante_p[1]=sum(map$chr==1)
				for(i in 2:length(ante_p))
				{
				 	ante_p[i]=sum(dataobj$map$chr==i)+ante_p[i-1]
				}
				#print(ante_p)
			}else{
				ante_p=NULL
			}

		}
		
	}
	
	if(op$model=="BayesA")
	{
		if(op$method=="MCMC")
		{
			if(op$ante) {cat("Model underdevelopment\n")}#res<-anteBayesAm(dataobj,op,ante_p,y,Z,trait)
			else res<-BayesM(dataobj,op,y,Z,X,trait=trait)
		}else{
			if(op$ante) {cat("Model underdevelopment\n")}#res<-anteBayesAe(dataobj,op,ante_p,y,Z,trait)
			else res<-BayesE(dataobj,op,y,Z,X,trait=trait)	
		}
		
	}
	
	if(op$model=="BayesB"){
		if(op$method=="MCMC")
		{
			if(op$ante) {cat("Model underdevelopment\n")}#res<-anteBayesBm(dataobj,op,ante_p,y,Z,trait)
			else res<-BayesM(dataobj,op,y,Z,X,trait=trait)
		}else{
			if(op$ante) {cat("Model underdevelopment\n")}#res<-anteBayesBe(dataobj,op,ante_p,y,Z,trait)
			else {cat("Model underdevelopment\n")} #res<-BayesBe(dataobj,op,y,Z,X,trait)	
		}
	}
	
	if(op$model=="rrBLUP"){
		if(op$method=="MCMC")
		{
			res<-BayesM(dataobj,op,y,Z,X=X,trait=trait)
		}else{
			res<-BayesE(dataobj,op,y,Z,X=X,trait=trait)	
			
		}
	}
	
	if(op$model=="SSVS")
		{
			if(op$method=="MCMC")
			{
				if(op$ante) {cat("Model underdevelopment\n")}#res<-anteBayesAm(dataobj,op,ante_p,y,Z,X,trait)
				else res<-BayesM(dataobj,op,y,Z,X,trait=trait)
			}else{
				if(op$ante) {cat("Model underdevelopment\n")}#res<-anteBayesAe(dataobj,op,ante_p,y,Z,X,trait)
				else res<-BayesE(dataobj,op,y,Z,X,trait=trait)	
			}

		}
	
	if(op$model=="IWBayesA")
	{
		if(op$method=="MCMC")
		{
			{cat("Model underdevelopment\n")}#res<-IWBayesAm(dataobj,op,y,Z,X,trait)
		}else{
			
			#res<-BayesAe(dataobj,op,trait)	
		}
		
	}
	
	return(res)
}





BATools <- function(y=NULL,Z=NULL,X=NULL,Ainv=NULL,W=NULL,model=NULL,param=NULL,nIter=2000,
      burnIn=1000,skip=10,seed=1,saveAt=NULL,R2=0.5,verbose=NULL,ncore=1,fixed=list(df0=T,s0=F,pi=T)){
  if(is.null(y)) stop("must have phenotype y to start the function")
  if(is.null(Z)) stop("must have genotype Z to start the function")
  if(any(is.na(Z)))
  { 
    stop("Z cannot have NA values")
  }
  
  if(!is.null(Z) && class(Z)!="matrix") {
    warning("Z is not a matrix, converting it to a matrix, this may cause issues")
    Z=as.matrix(Z)
  }
  
  nR=1
  #see if X needs to be processed here 
  if(is.null(X)){
    X = matrix(1,length(y),1)
    rownames(X)=names(y)
  }else{
    if(class(X)=="matrix"){
       if(dim(X)[2]>1) nR=nR+1
    }else stop("Fixed effect X must be a matrix")
  }
  
  #see if A needs to be processed here 
  if(!is.null(Ainv)){
    if(any(is.na(Ainv)))
    { 
      stop("Ainv cannot have NA values")
    }
    
   
    if(class(Ainv)=="matrix") {
      nR=nR+1
    }else stop("Ainv must be a matrix")
  }
  
  #see if W needs to be processed here 
  if(!is.null(W)){
    if(any(is.na(W)))
    { 
      stop("W cannot have NA values")
    }
    if(class(W)=="matrix") {
      nR=nR+1
    }else stop("W must be a matrix")
  }
    
    
  if(is.null(model) && !is.null(Z)) {
    cat("model is NULL for marker effects, use the default options to fit a RR model.\n")
    model=list(type="RR")
  } 
  
  if(!is.null(model) && !is.null(Z)){
    model_names=c("RR","BayesA","BayesB","SSVS","anteBayesA","anteBayesB","ssBayesA","ssBayesB","ssSSVS","IWBayesA","CDBayesA")
    if(!(model %in% model_names))
    {
      tmp<-paste(model_names, collapse=", ")
      stop("model must be ",tmp)
    }
    param=create.param(model,param,y,Z,Ainv=Ainv,W=W,R2=R2,nR=nR)
  }
  
  
  #default for saveAt
  if(is.null(saveAt)) saveAt=""
  
  #default for verbose
  if(is.null(verbose)) verbose=list(v=T,iter=100,to="screen")
  else{
    verbose_names=c("v","iter","to")
    if(!(prod(names(verbose) %in% verbose_names)))
    {
      tmp<-paste(verbose_names, collapse=", ")
      stop("Verbos parameter names must be ",tmp)
    }
    if(is.null(verbose$iter)) verbose$iter=100
    if(is.null(verbose$v)) verbose$v=T
    if(is.null(verbose$to)) verbose$to="screen"
    if(!is.numeric(verbose$iter)) stop("iter in verbose must be numeric")
    if(!is.logical(verbose$v)) stop("v in verbose must be logical")
    print_tos=c("screen","disk")
    if(!(prod(verbose$to %in% print_tos))) stop("to in verbose must be screen or disk")	
    
  }
  
  if(model %in% c("RR","BayesA","BayesB","SSVS")) {
    BayesM2(y=y,Z=y,X=y,Ainv=Ainv,W=W,model=model,param=param,nIter=nIter,
           burnIn=burnIn,skip=skip,seed=seed,saveAt=saveAt,R2=R2,
           verbose=verbose,ncore=ncore,fixed)
  }

  
  return(res)
}


AIREML<-function(y=NULL,Z=NULL,X=NULL,A=NULL,W=NULL,model=NULL,
                 maxIter=50,saveAt=NULL,R2=0.5){
 model=c("GBLUP","EMBayesA","EMSSVS")
}


##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "2.15.0"))
    stop("This package requires R 2.15.0 or later")
  assign(".BATools.home", file.path(library, pkg),
         pos=match("package:BATools", search()))
  BATools.version = "1.0.2 (2016-11-30), build 11"
  assign(".BATools.version", BATools.version, pos=match("package:BATools", search()))
  if(interactive())
  {
    packageStartupMessage(paste("Package 'BATools', ", BATools.version, ". ",sep=""),appendLF=TRUE)
    packageStartupMessage("Type 'help(BATools)' for summary information",appendLF=TRUE)
  }
  ###Add citation infomation
  invisible()
}

calc.MSx<-function(X){
  x2=apply(X,2,function(a) sum(a^2))
  sumMeanSq= sum((apply(X,2,mean))^2)
  MSx=sum(x2)/dim(X)[1]-sumMeanSq
}


create.param<-function(model,param,y,Z,Ainv=Ainv,W=W,R2=R2,nR=nR){
  #if(is.null(param)) {
  #  param=list(df0=4.2,df0e=4.2)
  #  cat(" Set degree of freedom of varG and varE to default value = 4.2 \n")
  #}
  param_names=c("s0","df0","R2","s0e","df0e","R2e","s0a","df0a","R2a","s0u","df0u","R2u")
  
  #For REML equivalent, set all df=-2 and scale=0, that's flat prior 
  
  if(!(prod(names(param) %in% param_names)))
  {
    tmp<-paste(param_names, collapse=", ")
    stop("model parameter names can only be ",tmp," for RR")
  }
  
  vary=var(y,na.rm=TRUE)
  
  if(is.null(param$df0)) {
    param$df0=4.2
    cat(" Set degree of freedom of varG to default value = 4.2 \n")
  }
  
  
  
  if(is.null(param$s0)){
    if(is.null(param$R2)) param$R2=R2/nR
    MSz=calc.MSx(Z)
    param$s0=vary*param$R2*(param$df0+2)/MSz
    cat(paste(" Set scale of varG to default value = ",param$s0," \n",sep=""))
  }
  

  if(is.null(param$df0e)) {
    param$df0e=4.2
    cat(" Set degree of freedom of varE to default value = 4.2 \n")
  }
  if(is.null(param$s0e)){
    param$s0e=vary*(1-R2)*(param$df0e)
    cat(paste(" Set scale of varE to default value = ",param$s0e," \n",sep=""))
  }

  if(is.null(param$df0a) && (!is.null(Ainv))) {
    param$df0a=4.2
    cat(" Set degree of freedom of varA to default value = 4.2 \n")
  }
  if(is.null(param$s0a) && (!is.null(W))){
    param$s0a=param$s0
    cat(paste(" Set scale of varA to default value = ",param$s0a," \n",sep=""))
  }
  
  if(is.null(param$df0a) && (!is.null(W))) {
    param$df0u=4.2
    cat(" Set degree of freedom of varU to default value = 4.2 \n")
  }
  if(is.null(param$s0u) && (!is.null(W))){
    if(is.null(param$R2)) param$R2u=R2/nR
    MSw=calc.MSx(W)
    param$s0u=vary*param$R2u*(param$df0+2)/MSw
    cat(paste(" Set scale of varU to default value = ",param$s0u," \n",sep=""))
  }
  
  param
}

check.BayesA<-function(model,W=W,nR=nR){
  if(is.null(param)) param=list(df0=4.2,df0e=4.2)
  model_names=c("varG","SNPeff","s0","df0","s0e","df0e","gwa","win")
  #For GBLUP, df0=-2,s0=0,df0e=-2,s0e=0
  if(!(prod(names(model) %in% model_names)))
  {
    tmp<-paste(model_names, collapse=", ")
    stop("model parameter names can only be ",tmp," for GBLUP")
  }
  if(is.null(model$df0)) model$df0=4.2
  if(is.null(model$df0e)) model$df0e=4.2
  
  
}


check.GBLUP<-function(model){
  model_names=c("type","s0","df0","s0e","df0e","gwa","win")
  #For GBLUP, df0=-2,s0=0,df0e=-2,s0e=0
  if(!(prod(names(model) %in% model_names)))
  {
    tmp<-paste(model_names, collapse=", ")
    stop("model parameter names can only be ",tmp," for GBLUP")
  }
    
    
}

