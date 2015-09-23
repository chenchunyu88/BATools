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
bafit <- function(dataobj=NULL,op=NULL,y=NULL,Z=NULL,X=NULL,trait=NULL){
	if(is.null(dataobj)) 
	{
		if(is.null(y) || is.null(Z)) stop("must have baData or y and Z to start the function")
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
	
	if(op$model=="BayesC")
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
  BATools.version = "0.1.0 (2015-09-17), build 10"
  assign(".BATools.version", BATools.version, pos=match("package:BATools", search()))
  if(interactive())
  {
    packageStartupMessage(paste("Package 'BATools', ", BATools.version, ". ",sep=""),appendLF=TRUE)
    packageStartupMessage("Type 'help(BATools)' for summary information",appendLF=TRUE)
  }
  invisible()
}



