#' baFit function can fit various Bayesian models including rrBLUP, BayesA/B/C, antedependence models and etc. Input data can be either a `baData` object or specify the fixed and random effects seperately.
#' @title Fitting various Bayesian models
#' @param formula A list of baData including phenotypes, genotypes and etc.
#' @param data A list of options to run Bayesian models
#' @param geno A numeric vector of phenotypes
#' @param genoid matrix of genotypes
#' @param randomFormula A string indicating the trait for analysis
#' @param map
#' @param PedAinv
#' @param options
#' @param train
#' @param GWA
#' @return The result of the analysis as a `ba` object, a list of estimate of fixed and random effects as well as variance components.
#' @examples \dontrun{
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
                 PedAinv=NULL,options=NULL,train=NULL,GWA=c("No","SNP","Win")){
  GWA<-match.arg(GWA)
  if(GWA!="No") {
    if(is.null(map)) stop("provide map for GWA")
    else{
      if(sum(colnames(map)%in%c("chr","pos","idw"))<3) stop("map should be a data.frame with colnames of chr, pos and idw")
    }
  }
  
  if(substr(options$model,1,4)=="ante"){
  	if(is.null(map)) stop("provide map for antedepedence model") else{
  		if(max(map$chr)>1){
  			Tmpmap=map[which(rownames(map)%in%colnames(geno)),]
			ante_p=rep(0,max(map$chr))
			ante_p[1]=sum(map$chr==1)
			for(i in 2:length(ante_p))
			{
			 	ante_p[i]=sum(map$chr==i)+ante_p[i-1]
			}
  		}else{
  			ante_p=NULL
  		}
  	}
  }
  
  
  id<-model.frame(genoid,data=data,na.action = na.pass)
  id <- eval(id, parent.frame())
  id <-as.character(t(id))
  if(length(id)!=length(unique(id))) stop("BATools does not support repeated measures yet, please keep one record per individual or contact us")
  
  if(is.null(options)) cat("no options inputed, use the default options.\n")
  
  if(substr(options$model,1,2)=="ss" || substr(options$model,5,6)=="ss"){
    if(is.null(PedAinv)) stop("Inverse of pedigree based additive relationship matrix is required for single-step approach")
    if(nrow(data)!=nrow(PedAinv)){
      A=Matrix::solve(PedAinv,sparse=TRUE,tol=1e-16)
      A=as.matrix(A)
      colnames(A)=rownames(A)=colnames(PedAinv)
      A=A[id,id]
      Ainv=solve(A)
    }else{
      Ainv=as.matrix(PedAinv)
    }
    idgeno<-id[which(id %in% rownames(geno))]
    M<-geno[idgeno,]
  }else{
    id<-id[which(id %in% rownames(geno))]
    M<-geno[id,]
    if(nrow(M)<nrow(data)){
      warning("The number of phenotyped individuals are larger than genotyped, 
              only genotyped individual will be used. Please consider to use single-step approach")
      data<- data %>% filter (id %in% rownames(M))
    }   
  }
  
  mf <- model.frame(formula, data = data, na.action = na.pass)
  mf <- eval(mf, parent.frame())
  y <- model.response(mf)
  
  if(nrow(M)>=length(y) && substr(options$model,1,2)=="ss") stop("Number of genotyped individual is not less than observation, cannot to single-step")
  
  names(y)=id
  X <- model.matrix(formula,data=data)
  rownames(X)=id
  if (!is.null(randomFormula)) {
    stop("random formula is not implemented yet, contact us for more information")
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
    names(vtrain)=names(y)
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
  if(options$model=="ssGBLUP"){
    res<-ssBayesE(op,y,M,X,vtrain,GWA,map,Ainv)
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
    res<-ssBayesM(op,y,M,X,vtrain,GWA,map,Ainv)
  }
  
  if(options$model=="ssBayesB"){
    res<-ssBayesM(op,y,M,X,vtrain,GWA,map,Ainv)
  }
  
  if(options$model=="ssSSVS"){
    res<-ssBayesM(op,y,M,X,vtrain,GWA,map,Ainv)
  }
  
  if(options$model=="anteBayesA"){
     res<-anteBayesAm(op,y,M,X,vtrain,GWA,map,ante_p)
  }
  
  if(options$model=="anteBayesB"){
    res<-anteBayesBm(op,y,M,X,vtrain,GWA,map,ante_p)
  }
  
  if(options$model=="anteSSVS"){
	  stop("anteSSVS not yet available")
  }
  
  if(options$model=="antessBayesA"){
    stop("antessBayesA not yet available")
  }
  
  if(options$model=="antessBayesB"){
    stop("antessBayesB not yet available")
  }
  
  if(options$model=="antessSSVS"){
    stop("antessSSVS not yet available")
  }
  if(op$model=="IWBayesA"){
    
  }  
  res
  #list(y=y,X=X,M=M,Z=Z,id=id)
  }

#res<-baTest(driploss~sex,data=PigPheno)#,geno = PigM,~id,~litter+slgdt_cd,options=list(model="BayesA"))
#' @export
getMatrix<-function(formula, data, geno, genoid,randomFormula=NULL,map=NULL,PedAinv=NULL,options=NULL,train=NULL){
  id<-model.frame(genoid,data=data,na.action = na.pass)
  id <- eval(id, parent.frame())
  id <-as.character(t(id))
  if(length(id)!=length(unique(id))) stop("BATools does not support repeated measures yet, please keep one record per individual or contact us")
  
  if(is.null(options)) cat("no options inputed, use the default options.\n")
  
 
  
  if(substr(options$model,1,2)=="ss" | substr(options$model,5,6)=="ss"){
    if(is.null(PedAinv)) stop("Pedigree based additive relationship matrix is required for single-step approach")
    if(nrow(data)!=nrow(PedAinv)){
      A=Matrix::solve(PedAinv,sparse=TRUE,tol=1e-16)
      A=as.matrix(A)
      colnames(A)=rownames(A)=colnames(PedAinv)
      A=A[id,id]
      Ainv=solve(A)
    }else{
      Ainv=as.matrix(PedAinv)
    }
    idgeno<-id[which(id %in% rownames(geno))]
    M<-geno[idgeno,]
    if(GWA!="No" && !is.null(vtrain)) stop("Cannot do both GWA and cross-validation at the same time")
    
  }else{
    id<-id[which(id %in% rownames(geno))]
    M<-geno[id,]
    if(nrow(M)<nrow(data)){
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
    names(vtrain)=names(y)
  }else vtrain=NULL
  list(y=y,X=X,M=M,Z=Z,id=id,vtrain=vtrain,Ainv=Ainv)
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
  BATools.version = "0.2.3 (2017-07-02), build 13"
  assign(".BATools.version", BATools.version, pos=match("package:BATools", search()))
  if(interactive())
  {
    packageStartupMessage(paste("Package 'BATools', ", BATools.version, ". ",sep=""),appendLF=TRUE)
    packageStartupMessage("Type 'help(BATools)' for summary information and citations",appendLF=TRUE)
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

