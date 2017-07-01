#' bafitTest function can fit various Bayesian models including rrBLUP, BayesA/B/C, antedependence models and etc. Input data can be either a `baData` object or specify the fixed and random effects seperately.
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
#' @export
#' @useDynLib BATools

baTest<-function(formula, data, geno, genoid,randomFormula=NULL,map=NULL,PedA=NULL,options=NULL,train=NULL,GWA=FALSE){
  id<-model.frame(genoid,data=data,na.action = na.pass)
  id <- eval(id, parent.frame())
  id <-as.character(t(id))
  if(length(id)!=length(unique(id))) stop("BATools does not support repeated measures yet, please keep one record per individual or contact us")
  
  if(is.null(options)) cat("no options inputed, use the default options.\n")
  
  if(substr(options$model,1,2)=="ss" | substr(options$model,5,6)=="ss"){
    if(is.null(PedA)) stop("Pedigree based additive relationship matrix is required for single-step approach")
    
  }else{
    M<-geno[id,]
    if(nrow(M)<length(id)){
      warning("The number of phenotyped individuals are larger than genotyped, 
               only genotyped individual will be used. Please consider to use single-step approach")
      data<- data %>% filter (id %in% rownames(M))
    }   
  }
  
  mf <- model.frame(formula, data = data, na.action = na.pass)
  mf <- eval(mf, parent.frame())
  y <- model.response(mf)
  names(y)=id
  X <- model.matrix(formula,data=data)
  rownames(X)=id
  if (!is.null(randomFormula)) {
    reff <- model.frame(randomFormula, data = data, na.action = na.pass)
    reff <- eval(reff, parent.frame())
    Z<-list()
    if (ncol(reff) == 1) {
      Z[[1]]=model.matrix(as.formula(paste0("y~0+",colnames(reff)[1])),data=data)
    }else {
      for(i in 1:ncol(reff)){
        Z[[i]]=model.matrix(as.formula(paste0("y~0+",colnames(reff)[i])),data=data)
      }
    }
  }else {
    Z=NULL
  }
  if(!is.null(train)){
    vtrain=model.frame(train, data = data, na.action = na.pass)
    vtrain <-as.vector(t(vtrain))
  }else vtrain=NULL

  if(options$model=="rrBLUP"){
    if(options$method=="MCMC") res<-BayesM(op,y,M,X,vtrain,GWA)
    if(options$method %in% c("REML","MAP")) {
      if(dim(M)[1]> dim(M)[2]) res<-BayesE(op,y,M,X,vtrain,GWA) else res<-aBayesE(op,y,M,X,vtrain,GWA)
    } 
  }
  
  if(options$model=="GBLUP"){
    res<-aBayesE(op,y,M,X,vtrain,GWA)
  }
  
  if(options$model=="BayesA"){
    if(options$method=="MCMC") res<-BayesM(op,y,M,X,vtrain,GWA)
    if(options$method =="MAP") {
      if(dim(M)[1]> dim(M)[2]) res<-BayesE(op,y,M,X,vtrain,GWA) else res<-aBayesE(op,y,M,X,vtrain,GWA)
    } 
  }
  
  if(options$model=="BayesB"){
    if(options$method=="MCMC") res<-BayesM(op,y,M,X,vtrain,GWA)
    if(options$method=="MAP") stop("BayesB MAP approach is not available yet")
  }
  
  if(options$model=="SSVS"){
    if(options$method=="MCMC") res<-BayesM(op,y,M,X,vtrain,GWA)
    if(options$method=="MAP") {
      if(dim(M)[1]> dim(M)[2]) res<-BayesE(op,y,M,X,vtrain,GWA) else res<-aBayesE(op,y,M,X,vtrain,GWA)
    } 
  }
  
  if(options$model=="ssBayesA"){
    
  }
  
  if(options$model=="ssBayesB"){
    
  }
  
  if(options$model=="ssSSVS"){
    
  }
  
  if(options$model=="anteBayesA"){
    
  }
  
  if(options$model=="anteBayesB"){
    
  }
  
  if(options$model=="anteSSVS"){
    
  }
  
  if(options$model=="antessBayesA"){
    
  }
  
  if(options$model=="antessBayesB"){
    
  }
  
  if(options$model=="antessSSVS"){
    
  }
  if(op$model=="IWBayesA"){
  
  }  
  res
  #list(y=y,X=X,M=M,Z=Z,id=id)
}

#res<-baTest(driploss~sex,data=PigPheno)#,geno = PigM,~id,~litter+slgdt_cd,options=list(model="BayesA"))
#' @export
getMatrix<-function(formula, data, geno, genoid,randomFormula=NULL,map=NULL,PedA=NULL,options=NULL,train=NULL){
  id<-model.frame(genoid,data=data,na.action = na.pass)
  id <- eval(id, parent.frame())
  id <-as.character(t(id))
  if(length(id)!=length(unique(id))) stop("BATools does not support repeated measures yet, please keep one record per individual or contact us")
  
  if(is.null(options)) cat("no options inputed, use the default options.\n")
  
  if(substr(options$model,1,2)=="ss" | substr(options$model,5,6)=="ss"){
    if(is.null(PedA)) stop("Pedigree based additive relationship matrix is required for single-step approach")
    
  }else{
    M<-geno[id,]
    if(nrow(M)<length(id)){
      warning("The number of phenotyped individuals are larger than genotyped, 
              only genotyped individual will be used. Please consider to use single-step approach")
      data<- data %>% filter (id %in% rownames(M))
    }   
  }
  
  mf <- model.frame(formula, data = data, na.action = na.pass)
  mf <- eval(mf, parent.frame())
  y <- model.response(mf)
  names(y)=id
  X <- model.matrix(formula,data=data)
  rownames(X)=id
  if (!is.null(randomFormula)) {
    reff <- model.frame(randomFormula, data = data, na.action = na.pass)
    reff <- eval(reff, parent.frame())
    Z<-list()
    if (ncol(reff) == 1) {
      Z[[1]]=model.matrix(as.formula(paste0("y~0+",colnames(reff)[1])),data=data)
    }else {
      for(i in 1:ncol(reff)){
        Z[[i]]=model.matrix(as.formula(paste0("y~0+",colnames(reff)[i])),data=data)
      }
    }
  }else {
    Z=NULL
  }
  if(!is.null(train)){
  vtrain=model.frame(train, data = data, na.action = na.pass)
  vtrain <-as.vector(t(vtrain))
  }else vtrain=NULL
  list(y=y,X=X,M=M,Z=Z,id=id,vtrain=vtrain)
  }

