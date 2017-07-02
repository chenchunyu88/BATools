#This code demonstrate BayesB model
rm(list=ls())
library(BATools)
data("Pig")
#Standardize genotype matrix
geno=std_geno(PigM,method="s",freq=PigAlleleFreq)

init=set_init("driploss",data=PigPheno,geno=geno,"id",df=5,pi_snp=0.001,h2=0.5,c=NULL,model="BayesB",centered=TRUE)
run_para=list(niter=2000,burnIn=1000,skip=10)
print_mcmc=list(piter=500)
update_para=list(df=FALSE,scale=TRUE,pi=F)
op<-create.options(model="BayesB",method="MCMC",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="BayesB",print_mcmc=print_mcmc)

BB<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op,map=PigMap,GWA="Win")
BB
par(mfrow=c(1,2))
man_plot_prob(BB)
man_plot_prob(BB,type="Win")
#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(data = PigPheno,k=5,"driploss")
head(PigPheno)
cvBB<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1)
cvBB
plot(cvBB)
