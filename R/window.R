#' create fix size window  
#' @param map genomic map has at least two column with column name \code{chr} for chromosome id and \code{pos} for postition in the chromosome in base pair (bp)
#' @param len window length
#' @param unit the unit for window length can be either Mb or kb or count
#' @return a map with new column \code{idw}
#' @examples \dontrun{
#'  rm(list=ls())
#'  library(BATools)
#'  data("Pig")
#'  map<-PigMap
#'  map<-set.win(map = map,len=1,unit = "Mb")
#'  map<-set.win(map = map,len = 5,unit="count")
#' }
#' @export
set.win<-function(map,len=NULL,unit=c("Mb","kb","count"))
{
  unit<-match.arg(unit)
  if(unit=="Mb"){
    map$idw=get_win_bp(len,map,ldiff=1000000)
  }else if(unit=="kb")
  {
    map$idw=get_win_bp(len,map,ldiff=1000)
  }else if(unit=="count"){
    map$idw=get_win_snp(len,map)
  }
  map
}


get_win_snp<-function(len,map){
  win=c()
  chr_nwin0=0
  chr_nwin1=0
  for(i in 1:max(map$chr)){
    chr_len=sum(map$chr==i)
    chr_nwin1<-chr_len%/%len
    for(j in 1:chr_nwin1)
      win=c(win,rep(j+chr_nwin0,len))
    if(chr_len%%len!=0){
      chr_nwin1=chr_nwin1+1
      j=j+1
      win=c(win,rep(j+chr_nwin0,chr_len%%len))
    }
    chr_nwin0=chr_nwin0+chr_nwin1
  }
  win
}


get_win_bp<-function(len,map,start=0,ldiff=1000000){
  #start!=0 not tested
  pos=map$pos/ldiff
  chr=map$chr
  win=c()
  chr_nwin0=0
  for(i in 1:max(chr)){
    subpos=pos[which(chr==i)]
    upper_limit=max(subpos)%/%len*len+len
    if(max(subpos)%%len==0) upper_limit=upper_limit-len
    if(start==0) seqwin=seq(len,upper_limit,len) else{
      seqwin=seq(len+len+start,upper_limit+start,len)
      countSNP=sum(subpos<=(len+start))
      if(countSNP!=0) {
        win=c(win, rep(1+chr_nwin0,countSNP))
      }
    }
    jw=0
    for(j in seqwin){
      countSNP=sum(subpos>(j-len) & subpos<=j)
      if(countSNP!=0) {
        jw=jw+1
        win=c(win, rep(jw+chr_nwin0, countSNP))
      }
    }
    chr_nwin0=chr_nwin0+max(win)
  }
  win
}
