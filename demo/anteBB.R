#This code demonstrate BayesB model
rm(list=ls())
library(BATools)
data("Pig")
#Standardize genotype matrix
geno=std_geno(PigM,method="s",freq=PigAlleleFreq)

init=set.init(~driploss,data=PigPheno,geno=geno,~id,df=5,pi_snp=0.001,h2=0.5,c=NULL,model="anteBayesB",centered=TRUE)
run_para=list(niter=2000,burnIn=1000,skip=10)
print_mcmc=list(piter=500)
update_para=list(df=FALSE,scale=TRUE,pi=F)
op<-set.options(model="anteBayesB",method="MCMC",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="anteBayesB",print_mcmc=print_mcmc)

anteBB<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op,map=PigMap,GWA="Win")
anteBB
par(mfrow=c(1,2))
man_plot_prob(anteBB)
man_plot_prob(anteBB,type="Win")
#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(~driploss,data = PigPheno,k=5)
head(PigPheno)
cvanteBB<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1,map=PigMap)
cvanteBB
par(mfrow=c(1,1))
baplot(cvanteBB)
