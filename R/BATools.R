#' BATools
#' 
#' BATools: An R package that extends 'Bayesian Alphabet' for Whole Genome Prediction and Genome-wide association study
#'
#' @details 
#' \describe{
#'  Package: BATools
#'  
#'  Version: 1.1.3
#'  
#'  Date: 2017-07-04
#'  
#'  License: GPL-3
#' }
#'  
#' @author
#'  Chunyu Chen and Robert J. Tempelman
#'  
#'  Maintainer: Chunyu Chen <chench57@msu.edu> 
#' @seealso \code{\link{baFit}} for fitting the model  
#' @references 
#' W. Yang, R. J. Tempelman. 2012. A Bayesian Antedependence Model for Whole Genome Prediction. \emph{Genetics}, vol. 190 (4) pp. 1491-1501.
#' 
#' C. Chen, R. J. Tempelman. 2015. An Integrated Approach to Empirical Bayesian Whole Genome Prediction Modeling. \emph{JABES}, Vol. 20 Issue 4, p491-511.
#' 
#' W. Yang, C. Chen, R. J. Tempelman. 2015. Improving the computational efficiency of fully Bayes inference and assessing the effect of misspecification of hyperparameters in whole-genome prediction models. \emph{Genetics Selection Evolution}, Vol. 47	Issue 1.
#' 
#' C. Chen, J. P. Steibel and R. J. Tempelman. 2017. Genome Wide Association Analyses Based on Broadly Different Specifications for Prior Distributions, Genomic Windows, and Estimation Methods.\emph{Genetics}, https://doi.org/10.1534/genetics.117.202259.
#' 
#' @docType package
#' @name BATools
NULL





#' baFit function can fit various Bayesian models including rrBLUP/GBLUP, BayesA/B, SSVS, single-step SSVS/BayesA/B, antedependence models and etc.
#' @title Fitting various Bayesian models
#' @param formula a formula for the fixed effects
#' @param data a `data.frame` containing the phenotypes, should have at least two columns including phenotype and id of the individual
#' @param geno a matrix of genotypes with rownames correpsonds to the id in data
#' @param genoid a forumla indicating the column of `data` contains the IDs related to genotype
#' @param randomFormula a formula for random effects 
#' @param map genomic map with column of chr (chromosome number) and pos (postion); for window based approach idw (window id) is also required. Refer to PigMap in the pig data as an example 
#' @param PedAinv the inverse of the additive relationship matrix based on pedigree; this is for single-step approach only
#' @param options the `options` object created by `create.options` function to run Bayesian models
#' @param train a formula indicating the column of `data` as reference for trainning and validation
#' @param GWA a character of in one of c("No","SNP","GWA") indicating what type of GWA for the model; default is "No"; `map` is required for both "SNP" and "Win"
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
  BATools.version = "1.1.3 (2017-07-04), build 13"
  assign(".BATools.version", BATools.version, pos=match("package:BATools", search()))
  if(interactive())
  {
    packageStartupMessage(paste("Package 'BATools', ", BATools.version, ". ",sep=""),appendLF=TRUE)
    packageStartupMessage("Type 'help(BATools)' for summary information and citations",appendLF=TRUE)
  }
  ###Add citation infomation
  invisible()
}

