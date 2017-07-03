#' A dataset contains an example of MSUPRP F2 Pig data 
#'
#' A dataset containing the MSUPRP example dataset with 176 pig with 20597 SNPs have both genotype and phenotype infomation.
#' The example data is a subset of MSUPRP (Gualdrón Duarte et al. 2014). The phenotype used in the example is drip loss in pork.
#'
#' @format 
#' \itemize{
#'   \item PigPheno \code{data.frame} of the phenotypes
#'   \item PigM  \code{matrix} of SNP marker genotypes coded as 0,1,2
#'   \item PigMap \code{data.frame} containing chromosome, postition and window for SNPs in PigM
#'   \item PigPed \code{data.frame} of pedigree of the pig population
#'   \item PigAlleleFreq allele frequency of F0
#'   \item PigAinv \code{matrix} of the inverse of the additive relationship matrix based on Pig pedigree
#'   \item PigA \code{matrix} of the additive relationship matrix based on Pig pedigree
#' }
#' @docType data
#' @keywords datasets
#' @name Pig
#' @usage data(Pig)
#' @source \url{https://github.com/steibelj/gwaR/blob/master/data/MSUPRP_sample.RData}
#' @source How the data is prepared for BATools: \url{https://github.com/chenchunyu88/batoolsdata/blob/master/MSUPRP.ipynb}
#' @references Jose L Gualdron Duarte, Rodolfo JC Cantet, Ronald O Bates, Catherine W Ernst, Nancy E Raney and Juan P Steibel. 2014. Rapid screening for phenotype-genotype associations by linear transformations of genomic evaluations. BMC Bioinformatics, 15:246.
NULL