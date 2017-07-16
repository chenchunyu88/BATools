#This code demonstrate BayesA model
rm(list=ls())
library(BATools)
data("Pig")
#Standardize genotype matrix
geno=std_geno(PigM,method="s",freq=PigAlleleFreq)

#Mask some genotype as missing to test single-step approach
set.seed(1001)
n=dim(geno)[1]
indexng<-sort(sample(1:n,n%/%5))
genoNew=geno[-indexng,]

init=set.init(~driploss,data=PigPheno,geno=genoNew,~id,df=5,pi_snp=1,h2=0.5,c=NULL,model="ssBayesA",centered=TRUE)
#or set your own starting values using 
#init=list(df=5,scale=0.01,pi=1) 
run_para=list(niter=2000,burnIn=1000,skip=10)
print_mcmc=list(piter=500)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-set.options(model="ssBayesA",method="MCMC",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="ssBayesA",print_mcmc=print_mcmc)

ssBA<-baFit(driploss~sex,data=PigPheno,geno=genoNew ,genoid = ~id,options = op,map=PigMap,GWA="Win",ped = PigPed)
ssBA
par(mfrow=c(1,2))
man_plot_prob(ssBA)
man_plot_prob(ssBA,type="Win")
#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(~driploss,data = PigPheno,k=5)
cvssBA<-baFit(driploss~sex,data=PigPheno,geno=genoNew ,genoid = ~id,options = op, train=~cv1,ped=PigPed)
cvssBA
par(mfrow=c(1,1))
baplot(cvssBA)
