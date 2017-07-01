#This code demonstrate GBLUP model
rm(list=ls())
library(BATools)
data("Pig")
#Standardize genotype matrix
geno=std_geno(PigM,method="s",freq=PigAlleleFreq)

init=set_init("driploss",data=PigPheno,geno=geno,"id",df=5,pi_snp=1,h2=0.5,c=NULL,model="GBLUP",centered=TRUE)
#or set your own starting values using 
#init=list(df=5,scale=0.01,pi=1) 
run_para=list(maxiter=100)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-create.options(model="GBLUP",method="REML",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="GBLUP",print_mcmc=NULL)

###Tested it's the same with other REML packages using the default settings
gblup<-baTest(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op)
gblup
#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(data = PigPheno,k=5,"driploss")
cvgblup<-baTest(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1)
cvgblup
plot(cvgblup)