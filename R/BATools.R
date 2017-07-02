#' baFit function can fit various Bayesian models including rrBLUP, BayesA/B/C, antedependence models and etc. Input data can be either a `baData` object or specify the fixed and random effects seperately.
#' @title Fitting various Bayesian models
#' @param dataobj A list of baData including phenotypes, genotypes and etc.
#' @param op A list of options to run Bayesian models
#' @param y A numeric vector of phenotypes
#' @param A matrix of genotypes
#' @param trait A string indicating the trait for analysis
#' @return The result of the analysis as a `ba` object, a list of estimate of fixed and random effects as well as variance components.
#' @examples \dontrun{
#' ###########Loading and preparing data##########
#' #This code demonstrate SSVS model
#' rm(list=ls())
#' library(BATools)
#' data("Pig")
#' #Standardize genotype matrix
#' geno=std_geno(PigM,method="s",freq=PigAlleleFreq)
#' init=set_init("driploss",data=PigPheno,geno=geno,"id",df=5,pi_snp=0.001,h2=0.5,c=1000,model="SSVS",centered=TRUE)
#' #or set your own starting values using 
#' #init=list(df=5,scale=0.01,pi=1) 
#' run_para=list(niter=20000,burnIn=10000,skip=10)
#' print_mcmc=list(piter=1000)
#' update_para=list(df=FALSE,scale=TRUE,pi=F)
#' op<-create.options(model="SSVS",method="MCMC",priors=NULL,init=init,
#'                    update_para=update_para,run_para=run_para,save.at="SSVS",print_mcmc=print_mcmc)
#' 
#' SSVS<-baTest(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op)
#' SSVS
#' #### Cross-validation using BATools
#' set.seed(1234)
#' PigPheno=createCV(data = PigPheno,k=5,"driploss")
#' cvSSVS<-baTest(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1)
#' plot(cvSSVS)
#' cvSSVS
#' 
#' #For GBLUP/rrBLUP:
#' demo(GBLUP)
#' demo(RR)
#' 
#' #For BayesA:
#' demo(BA)
#' demo(mapBA)
#' 
#' #For SSVS:
#' demo(SSVS)
#' demo(mapSSVS)
#' 
#' #For BayesB:
#' demo(BB)
#' 
#' }
#' @export
#' @useDynLib BATools
baFit<-function(formula, data, geno, genoid,randomFormula=NULL,map=NULL,
                 PedA=NULL,options=NULL,train=NULL,GWA=c("No","SNP","Win")){
  GWA<-match.arg(GWA)
  if(GWA!="No") {
    if(is.null(map)) stop("provide map for GWA")
    else{
      if(sum(colnames(map)%in%c("chr","pos","idw"))<3) stop("map should be a data.frame with colnames of chr, pos and idw")
    }
  }
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
    if(options$method=="MCMC") res<-BayesM(op,y,M,X,vtrain,GWA,map)
    if(options$method %in% c("REML","MAP")) {
      if(dim(M)[1]> dim(M)[2]) stop("nObservation>nSNP is not available for GBLUP")#res<-BayesE(op,y,M,X,vtrain,GWA,map) 
      else res<-aBayesE(op,y,M,X,vtrain,GWA,map)
    } 
  }
  
  if(options$model=="GBLUP"){
    res<-aBayesE(op,y,M,X,vtrain,GWA,map)
  }
  
  if(options$model=="BayesA"){
    if(options$method=="MCMC") res<-BayesM(op,y,M,X,vtrain,GWA,map)
    if(options$method =="MAP") {
      if(dim(M)[1]> dim(M)[2]) stop("nObservation>nSNP is not available for MAP")#res<-BayesE(op,y,M,X,vtrain,GWA,map) 
      else res<-aBayesE(op,y,M,X,vtrain,GWA,map)
    } 
  }
  
  if(options$model=="BayesB"){
    if(options$method=="MCMC") res<-BayesM(op,y,M,X,vtrain,GWA,map)
    if(options$method=="MAP") stop("BayesB MAP approach is not available yet")
  }
  
  if(options$model=="SSVS"){
    if(options$method=="MCMC") res<-BayesM(op,y,M,X,vtrain,GWA,map)
    if(options$method=="MAP") {
      if(dim(M)[1]> dim(M)[2]) stop("nObservation>nSNP is not available for MAP")#res<-BayesE(op,y,M,X,vtrain,GWA,map) 
      else res<-aBayesE(op,y,M,X,vtrain,GWA,map)
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

