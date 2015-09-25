# This function combines all raw data sources and/or pre-created synbreed data object in a single, unified data object of class baData. This is a list with elements for phenotypic, genotypic, marker map, pedigree and further covariate data. All elements are optional.
#' @title create \code{\link{ba}} object
#' @param synbreedobj A pre-created synbreed data object
#' @param pheno \code{data.frame} with individuals organized in rows and traits organized in columns. For unrepeated measures unique rownames should identify individuals. For repeated measures, the first column identifies individuals and a second column indicates repetitions (see also argument \code{repeated}).
#' @param geno \code{matrix} with individuals organized in rows and markers organized in columns. Genotypes could be coded arbitrarily. Missing values should be coded as \code{NA}. Colums or rows with only missing values not allowed. Unique \code{rownames} identify individuals and unique \code{colnames} markers. If no \code{rownames} are available, they are taken from element \code{pheno} (if available and if dimension matches).  If no \code{colnames} are used, the \code{rownames} of \code{map} are used if dimension matches.
#' @param map \code{data.frame} with one row for each marker and two columns (named \code{chr} and \code{pos}). First columns gives the chromosome (\code{numeric} or \code{character} but not \code{factor}) and second column the position  on the chromosome in centimorgan or the physical distance relative to the reference sequence in basepairs. Unique \code{rownames} indicate the marker names which should match with marker names in \code{geno}. Note that order and number of markers must not be identical with the order in \code{geno}. If this is the case, gaps in the map are filled with \code{NA} to ensure the same number and order as in element \code{geno} of the resulting \code{gpData} object.
#' @param pedigree Object of class \code{pedigree}.
#' @param family \code{data.frame} assigning individuals to families with names of individuals in \code{rownames} This information could be used for replacing of missing values with function \code{codeGeno} in synbreed.
#' @param covar \code{data.frame} with further covariates for all individuals that either appear in \code{pheno}, \code{geno} or \code{pedigree$ID},  e.g. sex or age. \code{rownames} must be specified to identify individuals. Typically this element is not specified by the user.
#' @param reorderMap \code{logical}. Should markers in \code{geno} and \code{map} be reordered by chromosome number and position within chromosome according to \code{map} (default = \code{TRUE})?
#' @param map.unit Character. Unit of position in \code{map}, i.e. 'cM' for genetic distance or 'bp' for physical distance (default = 'cM').
#' @param repeated This column is used to identify the replications of the phenotypic values. The unique values become the names of the third dimension of the pheno object in the \code{gpData}. This argument is only required for repeated measurements.
#' @param modCovar \code{vector} with \code{colnames} which identify columns with covariables in \code{pheno}. This argument is only required for repeated measurements.
#' @param A \code{matrix} of additive genetic relationship matrix 
#' @param Ainv \code{matrix} of inverse of A
#' @param fixed \code{matrix} or \code{factor} of fixed effects
#' @param random \code{matrix} or \code{factor} of random effects
#' @param makeA  \code{logical} indicating if create A matrix based on the pedigree information
#' @param makeAinv \code{logical} indicating if create A inverse matrix based on the pedigree information, A will be created before creating Ainv is created
#' @details The object \code{baData} is designed to provide a unified data object for whole genome prediction. It inherits the \code{gpData} object from synbreed package and is extended with more optional data item for whole genome prediction. So \code{baData} can be used in any functions in synbreed package including data imputation and etc. 
#' @examples 
#' #####create baData from gpData######
#' data("MSUPRP_sample") #load Michigan State University pig data
#' pheno<-data.frame(MSUPRP_sample$pheno[,,]) 
#' geno<-MSUPRP_sample$geno[,1:500]
#' ped<-MSUPRP_sample$pedigree
#' map=MSUPRP_sample$map
#' sex<-ped$sex
#' sex<-as.factor(sex)
#' x<-model.matrix( ~ sex -1,contrasts.arg=list(sex=contrasts(sex, contrasts=F)))
#' colnames(x)<-c("female","male")
#' rownames(x)<-ped$ID
#' #create baData object using synbreed object
#' pig.syn<-create.baData(synbreedobj=MSUPRP_sample,fixed=x,makeAinv=T) 
#' names(pig.syn)
#' #create baData object using raw data
#' pig=create.baData(pheno=pheno,geno=geno,map=map,pedigree=ped,fixed=x,makeAinv=F)
#'  @export
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
  else {
	  obj$fixed=matrix(rep(1,nrow(obj$pheno)),nrow=nrow(obj$pheno),ncol=1)
	  names(obj$fixed)=rownames(pheno)
  }

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


