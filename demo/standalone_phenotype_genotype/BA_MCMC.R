#demo for BayesA using MCMC
rm(list=ls())
library(BATools)
library(regress)
data("MSUPRP_sample")
summary(MSUPRP_sample)



pheno<-MSUPRP_sample$pheno
geno<-MSUPRP_sample$geno
ped<-MSUPRP_sample$pedigree
map=MSUPRP_sample$map

dim(map)
dim(geno)

index_chr17=which(map$chr==17)

geno17=geno[,rownames(map)[index_chr17]]

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
init=set_init(y=y,Z=Z,df=5,pi_snp=1,h2=0.5,c=NULL,model="BayesA",centered=T)
#or set your own starting values using list(df=5,scale=0.01,pi=1) 
run_para=list(niter=20000,burnIn=10000,skip=10)
print_mcmc=list(piter=2000)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-create.options(model="BayesA",method="MCMC",ante=FALSE,priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="BayesA",cv=NULL,print_mcmc=print_mcmc)

##################run BayesA MCMC###################


ba<-bafit(op=op,y=y,Z=Z,X=X)
ba
baplot(y=y,Z=Z,BAout=ba,type="pre")
baplot(y=y,Z=Z,BAout=ba,type="trace",op=op)

##################run BayesA EM#####################
init=list(df=5,scale=rr$hyper_est[2],vare=rr$hyper_est[1],g=rr$ghat,b=rr$betahat,pi=1)
run_para=list(maxiter=100)
#print_mcmc=list(piter=2000)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-create.options(model="BayesA",method="EM",ante=FALSE,priors=NULL,init=init,D="V",
    update_para=update_para,run_para=run_para,save.at="BayesA",cv=NULL,print_mcmc=NULL,convcrit=1E-4)

ba_em<-bafit(op=op,y=y,Z=Z,X=X)
ba_em
#ba_em<-BayesAe(dataobj=pig,op=op,trait="driploss")
baplot(y=y,Z=Z,BAout=ba_em,type="pre")

plot(ba$ghat,ba_em$ghat,xlab="MCMC",ylab="EM",main="BayesA MCMC v.s. EM")
abline(a=0,b=1)



##################run BayesC MCMC###################
init=list(df=5,scale=0.2,pi=0.1,c=1000)
run_para=list(niter=2000,burnIn=1000,skip=1)
print_mcmc=list(piter=200)
update_para=list(df=F,scale=T,pi=T)
priors=list(shape_scale=5,rate_scale=0.1,alphapi=1,betapi=9)
op<-create.options(model="SSVS",method="MCMC",
ante=FALSE,priors=NULL,init=init,update_para=update_para,
run_para=run_para,save.at="SSVS",cv=NULL,print_mcmc=print_mcmc)

bc<-bafit(op=op,y=y,Z=Z,X=X)


##################run BayesC EM#####################
init=list(df=5,scale=rr$hyper_est[2],vare=rr$hyper_est[1],g=rr$ghat,b=rr$betahat,pi=0.1,c=1000)
run_para=list(maxiter=100)
update_para=list(df=FALSE,scale=TRUE,pi=T)
op<-create.options(model="SSVS",method="EM",ante=FALSE,priors=NULL,init=init,
    update_para=update_para,run_para=run_para,save.at="SSVS",cv=NULL,print_mcmc=NULL,convcrit=1E-4)

bc_em<-bafit(op=op,y=y,Z=Z,X=X)

##################run BayesB MCMC###################
init=list(df=5,scale=0.01,pi=0.11)
run_para=list(niter=2000,burnIn=100,skip=10)
print_mcmc=list(piter=100,time_est=T,print_to="screen")
update_para=list(df=F,scale=TRUE,pi=TRUE)
op<-create.options(model="BayesB",method="MCMC",ante=FALSE,priors=NULL,init=init,
update_para=update_para,run_para=run_para,save.at="BayesB",cv=NULL,print_mcmc=print_mcmc)

bb<-bafit(op=op,y=y,Z=Z,X=X)





