#This code demonstrate BayesA model
rm(list=ls())
library(BATools)
data("Pig")
#Standardize genotype matrix
geno=std_geno(PigM,method="s",freq=PigAlleleFreq)


init=set.init(~driploss,data=PigPheno,geno=geno,~id,df=5,pi_snp=1,h2=0.5,c=NULL,model="BayesA",centered=TRUE)
#or set your own starting values using 
#init=list(df=5,scale=0.01,pi=1) 
run_para=list(niter=2000,burnIn=1000,skip=10)
print_mcmc=list(piter=500)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-set.options(model="BayesA",method="MCMC",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="BayesA",print_mcmc=print_mcmc)

BA<-baFit(driploss~1,data=PigPheno,geno=geno ,genoid = ~id,options = op,map=PigMap,GWA="Win")
BA
par(mfrow=c(1,2))
man_plot_prob(BA)
man_plot_prob(BA,type="Win")
#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(~driploss,data = PigPheno,k=5)
cvBA<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1)
cvBA
par(mfrow=c(1,1))
baplot(cvBA)