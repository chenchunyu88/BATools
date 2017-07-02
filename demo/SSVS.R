#This code demonstrate SSVS model
rm(list=ls())
library(BATools)
data("Pig")
#Standardize genotype matrix
geno=std_geno(PigM,method="s",freq=PigAlleleFreq)
init=set_init("driploss",data=PigPheno,geno=geno,"id",df=5,pi_snp=0.001,h2=0.5,c=1000,model="SSVS",centered=TRUE)
#or set your own starting values using 
#init=list(df=5,scale=0.01,pi=1) 
run_para=list(niter=2000,burnIn=1000,skip=10)
print_mcmc=list(piter=500)
update_para=list(df=FALSE,scale=TRUE,pi=F)
op<-create.options(model="SSVS",method="MCMC",seed=1,priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="SSVS",print_mcmc=print_mcmc)

SSVS<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op,map=PigMap,GWA="Win")
SSVS

par(mfrow=c(1,2))
man_plot_prob(SSVS)
man_plot_prob(SSVS,type="Win")


#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(data = PigPheno,k=5,"driploss")
cvSSVS<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1)
plot(cvSSVS)
cvSSVS