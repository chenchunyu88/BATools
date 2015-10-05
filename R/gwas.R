
#calculate SNP effect variance
varsnp<-function(zf,iG,iGCG){
      t1<-zf%*%iG%*%zf
      t2<-zf%*%iGCG%*%zf
      return(c(t1,t2))
}

#' A function to get p-values from \code{`ba`} object 
#' @title get p-values 
#'  @export
get_pvalues<-function(data,type="random"){
	g<-data$ghat	
	if(type=="random"){
		varg<-diag(data$Cgg)
		zscore=g/sqrt(abs(varg))
		pvalues=2*(1-pnorm(abs(zscore)))
	}
	if(type=="fixed"){
		varg<-data$hyper_est[2]-diag(data$Cgg)
		zscore=g/sqrt(abs(varg))
		pvalues=2*(1-pnorm(abs(zscore)))
	}
	pvalues
}

#' A function to create manhattan plots from GWA and QTL mapping analyses.  
#' @title Manhattan plots of p-values 

#' @param pvalues A vector of pvalues with SNP, marker or position names
#' @param map a genetic map: a data.frame with SNP or marker names in rows, chromosome name (chr column) and position (pos column)
#' @param threshold the threshold of significance 
#' @param col a vector of strings with color names for painting SNP according to their chromosome
#' @param ... additional graphical parameters used by plot 
#' @return NULL. A graphic is be generated
#'  @export
manhattan_plot <- function(pvalues, map, threshold = 0.01, col = c("black", "red"), chrom = NULL, ...) {
    if (!("chr" %in% colnames(map))) 
        stop("chr should be a column of map")
    if (sum(rownames(map) == names(pvalues)) < length(pvalues)) 
        stop("make sure that the values are ordered, rownames of map and names of p-value vector should agree")
    if ((max(pvalues) > 1) | (min(pvalues) < 0)) 
        stop("pvalues should be in the range [0,1]")
    if (sum(pvalues == 0) > 0) {
        pvalues[pvalues == 0] <- min(pvalues[pvalues > 0])
        warning("some p-values were exactly zero and they were replaced with the next smallest value")
    }
    if (!is.null(chrom)) {
        idx <- as.character(map$chr) %in% as.character(chrom)
        map <- map[idx, ]
        pvalues <- pvalues[idx]
    }
    if (length(col) < 2) 
        stop("two colors are needed in the color vector")
    mps <- as.numeric(as.factor(map$chr))
    plot(-log(pvalues, 10), pch = 20, col = col[(mps%%2) + 1], ylab = "-log10-pvalue", xlab = "chr", ..., axes = F)
    axis(2)
    lns <- (by(map, map$chr, nrow))
    axis(1, at = c(0, cumsum(lns)[-length(lns)]) + as.vector(lns/2), labels = names(lns))
    box()
    abline(h = -log(threshold, 10), lwd = 2, col = "green")
}


#' A function to create posterior probability plots from GWA and QTL mapping analyses.  
#' @title Manhattan plots of posterior probability 

#' @param postprob A vector of posterior probability with SNP, marker or position names
#' @param map a genetic map: a data.frame with SNP or marker names in rows, chromosome name (chr column) and position (pos column)
#' @param col a vector of strings with color names for painting SNP according to their chromosome
#' @param ... additional graphical parameters used by plot 
#' @return NULL. A graphic is be generated
#'  @export
postprob_plot <- function(postprob, map, col = c("black", "red"), chrom = NULL, ...) {
  if (!("chr" %in% colnames(map))) 
    stop("chr should be a column of map")
  if (sum(rownames(map) == names(postprob)) < length(postprob)) 
    stop("make sure that the values are ordered, rownames of map and names of p-value vector should agree")
  if (!is.null(chrom)) {
    idx <- as.character(map$chr) %in% as.character(chrom)
    map <- map[idx, ]
    postprob <- postprob[idx]
  }
  if (length(col) < 2) 
    stop("two colors are needed in the color vector")
  mps <- as.numeric(as.factor(map$chr))
  plot(postprob, pch = 20, col = col[(mps%%2) + 1], ylab = "Posterior Probability", xlab = "chr", ..., axes = F)
  axis(2)
  lns <- (by(map, map$chr, nrow))
  axis(1, at = c(0, cumsum(lns)[-length(lns)]) + as.vector(lns/2), labels = names(lns))
  box()
}