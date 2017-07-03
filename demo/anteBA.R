#This code demonstrate BayesA model
rm(list=ls())
library(BATools)
data("Pig")
#Standardize genotype matrix
geno=std_geno(PigM,method="s",freq=PigAlleleFreq)


init=set_init("driploss",data=PigPheno,geno=geno,"id",df=5,pi_snp=1,h2=0.5,c=NULL,model="anteBayesA",centered=TRUE)
#or set your own starting values using 
#init=list(df=5,scale=0.01,pi=1) 
run_para=list(niter=2000,burnIn=1000,skip=10)
print_mcmc=list(piter=500)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-create.options(model="anteBayesA",method="MCMC",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="anteBayesA",print_mcmc=print_mcmc)

anteBA<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op,map=PigMap,GWA="Win")
anteBA
par(mfrow=c(2,2))
man_plot_prob(anteBA)
man_plot_prob(anteBA,type="Win")
plot(anteBA)
#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(data = PigPheno,k=5,"driploss")
cvanteBA<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1,map = PigMap)
cvanteBA
par(mfrow=c(1,1))
plot(cvanteBA)
