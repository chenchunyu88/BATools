#Vignette demo
rm(list=ls())
library(BATools)
library(regress)
data("MSUPRP_sample")
summary(MSUPRP_sample)



pheno<-MSUPRP_sample$pheno
geno<-MSUPRP_sample$geno[,1:20000]
ped<-MSUPRP_sample$pedigree
map=MSUPRP_sample$map

sex<-ped$sex
sex<-as.factor(sex)
x<-model.matrix( ~ sex -1,contrasts.arg=list(sex=contrasts(sex, contrasts=F)))
colnames(x)<-c("female","male")
#x<-model.matrix(~1,data=sex)
rownames(x)<-ped$ID

trait="driploss"
y=na.omit(pheno[,trait,1])

X=x
Z=std_geno(geno,method="s") #centering the genotype matrix
nx=rownames(X)
ng=rownames(Z)
np=names(y)
idx <- Reduce(intersect, list(nx,ng,np))
X=as.matrix(X[idx,],ncol=dimX)
y=y[idx]
Z=Z[idx,]


######################BayesA####################
###############Setting up options###############
init=list(df=5,scale=0.01,pi=1) #You can use the set_init function if you like
run_para=list(niter=2000,burnIn=1000,skip=10)
print_mcmc=list(piter=1000)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-create.options(model="BayesA",method="MCMC",ante=FALSE,priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="BayesA",cv=NULL,print_mcmc=print_mcmc)

##################run BayesA MCMC###################

system.time({
ba<-bafit(op=op,y=y,Z=Z,X=X)})
ba
baplot(y=y,Z=Z,BAout=ba,type="pre")
baplot(y=y,Z=Z,BAout=ba,type="trace",op=op)

library(BGLR)
nIter=2000;
burnIn=1000;
thin=3;
saveAt='';
S0=NULL;
weights=NULL;
R2=0.5;
ETA<-list(list(X=Z,model='BayesA'))
system.time({
fit_BA=BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
})
plot(fit_BA$yHat,y)



##################run BayesA EM#####################
library(rrBLUP)
rr<-mixed.solve(y=y,Z=Z,X=X)
init=list(df=5,scale=rr$Vu,vare=rr$Ve,g=rr$u,beta=rr$beta,pi=1)
run_para=list(maxiter=100)
#print_mcmc=list(piter=2000)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-create.options(model="BayesA",method="EM",ante=FALSE,priors=NULL,init=init,D="V",
                   update_para=update_para,run_para=run_para,save.at="BayesA",cv=NULL,print_mcmc=NULL,convcrit=1E-4)

ba_em<-bafit(op=op,y=y,Z=Z,X=X)
ba_em
#
system.time({ba_ema<-BayesEA(op=op,y=y,Z=Z,X=X)})
baplot(y=y,Z=Z,BAout=ba_em,type="pre")

plot(ba$ghat,ba_em$ghat,xlab="MCMC",ylab="EM",main="BayesA MCMC v.s. EM")
abline(a=0,b=1)



