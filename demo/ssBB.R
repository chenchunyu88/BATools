#This code demonstrate BayesB model
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


init=set_init("driploss",data=PigPheno,geno=genoNew,"id",df=5,pi_snp=0.001,h2=0.5,c=NULL,model="ssBayesB",centered=TRUE)
run_para=list(niter=2000,burnIn=1000,skip=10)
print_mcmc=list(piter=500)
update_para=list(df=FALSE,scale=TRUE,pi=F)
op<-create.options(model="ssBayesB",method="MCMC",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="ssBayesB",print_mcmc=print_mcmc)

ssBB<-baFit(driploss~sex,data=PigPheno,geno=genoNew ,genoid = ~id,options = op,map=PigMap,GWA="Win",PedAinv = PigAinv)
ssBB
par(mfrow=c(1,2))
man_plot_prob(ssBB)
man_plot_prob(ssBB,type="Win")
#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(data = PigPheno,k=5,"driploss")
head(PigPheno)
cvssBB<-baFit(driploss~sex,data=PigPheno,geno=genoNew ,genoid = ~id,options = op, train=~cv1,PedAinv = PigAinv)
par(mfrow=c(1,1))
cvssBB
plot(cvssBB)
