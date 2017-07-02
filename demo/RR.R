#This code demonstrate rrBLUP model
rm(list=ls())
library(BATools)
data("Pig")
#Standardize genotype matrix
geno=std_geno(PigM,method="s",freq=PigAlleleFreq)

init=set_init("driploss",data=PigPheno,geno=geno,"id",df=5,pi_snp=1,h2=0.5,c=NULL,model="rrBLUP",centered=TRUE)
run_para=list(niter=2000,burnIn=1000,skip=10)
print_mcmc=list(piter=100)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-create.options(model="rrBLUP",method="MCMC",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="rrBLUP",print_mcmc=print_mcmc)

rrBLUP<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op)

#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(data = PigPheno,k=5,"driploss")
cvrrBLUP<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1)
plot(cvrrBLUP)
cvrrBLUP

