#This code demonstrate SSVS model
rm(list=ls())
devtools::install("../BATools",build_vignettes=T)
library(BATools)
data("Pig")
#Standardize genotype matrix
geno=std_geno(PigM,method="s",freq=PigAlleleFreq)
init=set_init("driploss",data=PigPheno,geno=geno,"id",df=5,pi_snp=0.001,h2=0.5,c=1000,model="SSVS",centered=TRUE)
#or set your own starting values using 
#init=list(df=5,scale=0.01,pi=1) 
run_para=list(niter=20000,burnIn=10000,skip=10)
print_mcmc=list(piter=1000)
update_para=list(df=FALSE,scale=TRUE,pi=F)
op<-create.options(model="SSVS",method="MCMC",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="SSVS",print_mcmc=print_mcmc)

SSVS<-baTest(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op)
SSVS
#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(data = PigPheno,k=5,"driploss")
cvSSVS<-baTest(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1)
plot(cvSSVS)
cvSSVS