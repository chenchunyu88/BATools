#Vignette demo
rm(list=ls())
devtools::install("../BATools")
library(BATools)
library(regress)

data("MSUPRP_sample")
summary(MSUPRP_sample)



pheno<-data.frame(MSUPRP_sample$pheno[,,]) 
geno<-MSUPRP_sample$geno[,seq(1,dim(MSUPRP_sample$geno)[2],20)]
ped<-MSUPRP_sample$pedigree
map=MSUPRP_sample$map

sex<-ped$sex
sex<-as.factor(sex)
x<-model.matrix( ~ sex -1,contrasts.arg=list(sex=contrasts(sex, contrasts=F)))
colnames(x)<-c("female","male")
#x<-model.matrix(~1,data=sex)
rownames(x)<-ped$ID
pig_tmp=create.baData(pheno=pheno,geno=geno,map=map,pedigree=ped,fixed=x,makeAinv=F)

pig<-std_geno(pig_tmp,method="s") #centering the genotype matrix



##################run rrBLUP REML#####################
init=set_init(pig,df=5,scale=NULL,vare=NULL,pi_snp=1,h2=0.5,c=NULL,model="rrBLUP",centered=T,trait="driploss")
run_para=list(maxiter=100)
priors=list(nu_e=-2,tau2_e=0,nu_s=-2,tau2_s=0) #REML like prior, change nu_e=-1 if you want to use Gelman's prior
update_para=list(scale=T,vare=T)
op<-create.options(model="rrBLUP",method="EM",ante=FALSE,priors=NULL,init=init,
    update_para=update_para,run_para=run_para,save.at="rrBLUP",cv=NULL,print_mcmc=NULL,convcrit=1E-4)

rr<-bafit(dataobj=pig,op=op,trait="driploss") #about 5 minutes

prr1<-get_pvalues(rr)
prr2<-get_pvalues(rr,type="fixed")

manhattan_plot(prr1,pig$map,threshold = 0.001,main="rrBLUP random effects test")
manhattan_plot(prr2,pig$map,threshold = 0.001,main="rrBLUP fixed effects test")


##################run BayesA EM#####################
init=set_init(pig,df=5,scale=rr$hyper_est[2],vare=rr$hyper_est[1],g=rr$ghat,b=rr$betahat,
              pi_snp=1,h2=0.5,model="BayesA",centered=T,trait="driploss",from="rrBLUP")
run_para=list(maxiter=100)
priors=list(nu_e=-1,tau2_e=0,nu_s=-1,tau2_s=0) 
update_para=list(df=FALSE,scale=TRUE,vare=T,pi=FALSE)
op<-create.options(model="BayesA",method="EM",ante=FALSE,priors=priors,init=init,D="P",
    update_para=update_para,run_para=run_para,save.at="BayesA",cv=NULL,print_mcmc=NULL,convcrit=1E-4)

ba_em<-bafit(dataobj=pig,op=op,trait="driploss")
ba_em

plot(ba_em$ghat,rr$ghat)
abline(a=0,b=1)

pba<-get_pvalues(ba_em,type="random")

manhattan_plot(pba,pig$map,threshold = 0.05,main="BayesA random effect test")




##################run BayesC EM#####################
init=set_init(pig,df=5,scale=rr$hyper_est[2],vare=rr$hyper_est[1],g=rr$ghat,beta=rr$betahat,
              pi_snp=0.05,h2=0.5,c=1000,model="SSVS",centered=T,trait="driploss",from="rrBLUP")
run_para=list(maxiter=100)
priors=list(nu_e=-2,tau2_e=0,nu_s=-2,tau2_s=0)
update_para=list(df=FALSE,scale=TRUE,pi=T)
op<-create.options(model="SSVS",method="EM",ante=FALSE,priors=priors,init=init,
    update_para=update_para,run_para=run_para,save.at="SSVS",cv=NULL,print_mcmc=NULL,convcrit=1E-4)

bc_em<-bafit(dataobj=pig,op=op,trait="driploss")
#bc_em<-BayesE(dataobj=pig,op=op,trait="driploss")

pbc<-get_pvalues(bc_em,type="random")

manhattan_plot(pbc,pig$map,threshold = 0.05,main="SSVS random effect test")

##################run BayesB MCMC###################
init=set_init(pig,df=5,pi_snp=0.05,h2=0.5,c=NULL,model="BayesB",centered=T,trait="driploss")
run_para=list(niter=2000,burnIn=100,skip=10)
print_mcmc=list(piter=100,time_est=T,print_to="screen")
update_para=list(df=F,scale=TRUE,pi=TRUE)
op<-create.options(model="BayesB",method="MCMC",ante=FALSE,priors=NULL,init=init,
update_para=update_para,run_para=run_para,save.at="BayesB",cv=NULL,print_mcmc=print_mcmc)

bb<-bafit(dataobj=pig,op=op,trait="driploss")

postprob_plot(bb$prob,pig$map)

##################run BayesC MCMC###################
init=set_init(pig,df=5,pi_snp=0.05,h2=0.5,c=1000,model="SSVS",centered=T,trait="driploss")
run_para=list(niter=20000,burnIn=10000,skip=10)
print_mcmc=list(piter=2000)
update_para=list(df=F,scale=T,pi=T)
priors=list(shape_scale=5,rate_scale=0.1,alphapi=1,betapi=9,nu_e=-2,tau2_e=0)
op<-create.options(model="SSVS",method="MCMC",
                   ante=FALSE,priors=priors,init=init,update_para=update_para,
                   run_para=run_para,save.at="SSVS",cv=NULL,print_mcmc=print_mcmc)

bc<-bafit(dataobj=pig,op=op,trait="driploss")
#bc<-BayesM(dataobj=pig,op=op,trait="driploss")
postprob_plot(bc$phisave,pig$map)

