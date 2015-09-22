# read genomic prediction data
#'  @export
#'  @title create \code{\link{ba}} object
create.baData <- function(synbreedobj=NULL,pheno=NULL,geno=NULL,map=NULL,pedigree=NULL,family=NULL,covar=NULL,
                          reorderMap=TRUE,map.unit="cM",repeated=NULL,modCovar=NULL,A=NULL,
			  Ainv=NULL,fixed=NULL,random=NULL,makeA=FALSE,makeAinv=FALSE){
  
  # create gpData first with some checks on data
  if(is.null(synbreedobj))
  {
	obj<-create.gpData(pheno=pheno,geno=geno,map=map,pedigree=pedigree,family=family,
	covar=covar,reorderMap=reorderMap,map.unit="cM",repeated=repeated,modCovar=modCovar)
  }else{
	if(class(synbreedobj)!="gpData") stop("synbreed data object must be a gpData")
	else obj=synbreedobj
  }	
  #add A to the object	
  if(!is.null(A)) {
	if(!is.matrix(A)) stop("A must be a matrix or data.frame, not a ", class(A))
    	if(nrow(A)!=nrow(pheno)) stop("A matrix must have the same number of row with phenotype, which should be ", nrow(obj$pheno))
	else{obj$A=A}  
  }
  else obj$A=NULL
	
  #add A inverse to the object	
  if(!is.null(Ainv)) {
	if(!is.matrix(Ainv)) stop("A inverse must be a matrix or data.frame, not a ", class(Ainv))
    	if(nrow(Ainv)!=nrow(obj$pheno)) stop("A inverse matrix must have the same number of row with phenotype, which should be ", nrow(obj$pheno))	
	else{obj$Ainv=Ainv}  
  }
  else obj$Ainv=NULL

  #add other fixed effect to the object	
  if(!is.null(fixed)) {
	if(class(fixed)=="factor") 
	{
		mat.fixed=model.matrix(~0+fixed)
		obj$fixed=mat.fixed
	}
	if(!is.matrix(fixed)) stop("Fixed effects must be a matrix, not a ", class(fixed))
    	#if(nrow(fixed)!=nrow(obj$pheno)) stop("Fixed effects matrix must have the same number of row with phenotype, which should be ", dim(obj$pheno)[1])
	    obj$fixed=fixed #else{obj$fixed=fixed}  
  }
  else obj$fixed=matrix(rep(1,nrow(obj$pheno)),nrow=nrow(obj$pheno),ncol=1)

  #add other random effect to the object	
  if(!is.null(random)) {
	if(class(random)=="factor") 
	{
		mat.random=model.matrix(~0+random)
		obj$random=mat.random
	}
	else 
	{
		if(!is.matrix(random)) stop("Random effects must be a matrix or data.frame, not a ", class(random))
    		if(nrow(random)!=nrow(pheno)) stop("Random effects matrix must have the same number of row with phenotype, which should be ", dim(obj$pheno)[1])
		obj$random= random	
	}


  }
  else obj$random=NULL

  	if(makeAinv)
	{
		obj$A=createA(obj)
		obj$Ainv=createAinv(obj$A)
	}
	
	if(makeA)
	{
		obj$A=createA(obj)	
	}

  # return object of class 'gpData'
  class(obj) <- "baData"
  return(obj)
}

arraytodataframe<-function(a1)
{
	if(class(a1)!="array") {
		stop("Input must be an array, not a ", class(random))
	}else
	{
		df1=data.frame(matrix(c(a1),nrow=dim(a1)[1],ncol=dim(a1)[2]))
		rownames(df1)=rownames(a1)
		colnames(df1)=colnames(a1)
	}
	return(df1)
}

#makeA
createA<-function(obj,index=NULL)
{
	cat("Creating A matrix...\n")
	Amat<-kin(obj,ret="add")
	if(is.null(index)) index_c<-which(rownames(Amat) %in% rownames(obj$geno))
	else index_c=index
	A<-Amat[index_c,index_c]
	return(A)	
}

createAinv<-function(A)
{
	cat("Creating A inverse matrix...\n")
	Ainv=solve(A)
	return(Ainv)
}


