#Vignette demo
rm(list=ls())
library(BATools)
library(regress)

data("MSUPRP_sample")
summary(MSUPRP_sample)



pheno<-data.frame(MSUPRP_sample$pheno[,,]) 
geno<-MSUPRP_sample$geno[,1:500]
ped<-MSUPRP_sample$pedigree
map=MSUPRP_sample$map

sex<-ped$sex
sex<-as.factor(sex)
x<-model.matrix( ~ sex -1,contrasts.arg=list(sex=contrasts(sex, contrasts=F)))
colnames(x)<-c("female","male")
#x<-model.matrix(~1,data=sex)
rownames(x)<-ped$ID
pig=create.baData(pheno=pheno,geno=geno,map=map,pedigree=ped,fixed=x,makeAinv=F)


##################run rrBLUP REML#####################
init=list(df=NULL,scale=NULL,vare=NULL)
run_para=list(maxiter=100)
update_para=list(df=FALSE,scale=TRUE)
op<-create.options(model="rrBLUP",method="EM",ante=FALSE,priors=NULL,init=init,
    update_para=update_para,run_para=run_para,save.at="rrBLUP",cv=NULL,print_mcmc=NULL,convcrit=1E-4)

rr<-bafit(dataobj=pig,op=op,trait="driploss")




##################run rrBLUP MCMC#####################
init=list(df=5,scale=NULL,vare=NULL)
run_para=list(niter=200000,burnIn=100000,skip=10)
print_mcmc=list(piter=10000)
update_para=list(df=FALSE,scale=TRUE)
op<-create.options(model="rrBLUP",method="MCMC",ante=FALSE,priors=NULL,init=init,
    update_para=update_para,run_para=run_para,save.at="rrBLUP",cv=NULL,print_mcmc=print_mcmc,convcrit=1E-4)
	

rrm<-bafit(dataobj=pig,op=op,trait="driploss")
rrm

plot(rrm$ghat,rr$ghat)
abline(a=0,b=1)

######################BayesA####################
###############Setting up options###############
init=list(df=5,scale=0.01,pi=1)
run_para=list(niter=20000,burnIn=10000,skip=10)
print_mcmc=list(piter=2000)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-create.options(model="BayesA",method="MCMC",ante=FALSE,priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="BayesA",cv=NULL,print_mcmc=print_mcmc)

##################run BayesA MCMC###################


ba<-bafit(dataobj=pig,op=op,trait="driploss")
ba

baplot(dataobj=pig,BAout=ba,type="pre")
baplot(dataobj=pig,BAout=ba,type="trace",op=op)

##################run BayesA EM#####################
init=list(df=5,scale=rr$hyper_est[2],vare=rr$hyper_est[1],g=rr$ghat,b=rr$betahat,pi=1)
run_para=list(maxiter=100)
#print_mcmc=list(piter=2000)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-create.options(model="BayesA",method="EM",ante=FALSE,priors=NULL,init=init,D="V",
    update_para=update_para,run_para=run_para,save.at="BayesA",cv=NULL,print_mcmc=NULL,convcrit=1E-4)

ba_em<-bafit(dataobj=pig,op=op,trait="driploss")
ba_em

baplot(dataobj=pig,BAout=ba_em,type="pre")

plot(ba$ghat,ba_em$ghat,xlab="MCMC",ylab="EM",main="BayesA MCMC v.s. EM")
abline(a=0,b=1)



##################run BayesC MCMC###################
init=list(df=5,scale=0.2,pi=0.1,c=1000)
run_para=list(niter=2000,burnIn=1000,skip=1)
print_mcmc=list(piter=200)
update_para=list(df=F,scale=T,pi=T)
priors=list(shape_scale=5,rate_scale=0.1,alphapi=1,betapi=9)
op<-create.options(model="BayesC",method="MCMC",
ante=FALSE,priors=NULL,init=init,update_para=update_para,
run_para=run_para,save.at="BayesCC",cv=NULL,print_mcmc=print_mcmc)

bc<-bafit(dataobj=pig,op=op,trait="driploss")

#bm<-BayesM(dataobj=pig,op=op,trait="driploss")


##################run BayesC EM#####################
init=list(df=5,scale=rr$hyper_est[2],vare=rr$hyper_est[1],g=rr$ghat,b=rr$betahat,pi=0.1,c=1000)
run_para=list(maxiter=100)
update_para=list(df=FALSE,scale=TRUE,pi=T)
op<-create.options(model="BayesC",method="EM",ante=FALSE,priors=NULL,init=init,
    update_para=update_para,run_para=run_para,save.at="BayesC",cv=NULL,print_mcmc=NULL,convcrit=1E-4)

#bc_em<-BayesCe(dataobj=pig,op=op,trait="driploss")
bc_em<-bafit(dataobj=pig,op=op,trait="driploss")

##################run BayesB MCMC###################
init=list(df=5,scale=0.01,pi=0.11)
run_para=list(niter=2000,burnIn=100,skip=10)
print_mcmc=list(piter=100,time_est=T,print_to="screen")
update_para=list(df=F,scale=TRUE,pi=TRUE)
op<-create.options(model="BayesB",method="MCMC",ante=FALSE,priors=NULL,init=init,
update_para=update_para,run_para=run_para,save.at="BayesB",cv=NULL,print_mcmc=print_mcmc)

bb<-bafit(dataobj=pig,op=op,trait="driploss")


