#This code demonstrate GBLUP model
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

init=set.init(~driploss,data=PigPheno,geno=genoNew,~id,df=5,pi_snp=1,h2=0.5,c=NULL,model="ssGBLUP",centered=TRUE)
#or set your own starting values using 
#init=list(df=5,scale=0.01,pi=1) 
run_para=list(maxiter=100)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-set.options(model="ssGBLUP",method="REML",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="ssGBLUP",print_mcmc=NULL,ssGBLUPvar = "homVAR")

###Tested it's the same with other REML packages using the default settings
ssgblup<-baFit(driploss~sex+car_wt,data=PigPheno,geno=genoNew ,genoid = ~id,ped=PigPed,options = op,map=PigMap,GWA="Win")
ssgblup
par(mfrow=c(1,2))
man_plot_pvalue(ssgblup)
man_plot_pvalue(ssgblup,type="Win")
#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(~driploss,data = PigPheno,k=5)
op$ssGBLUPvar="homVAR"
cvssgblup2<-baFit(driploss~sex,data=PigPheno,geno=genoNew ,genoid = ~id,ped=PigPed,options = op, train=~cv1)
par(mfrow=c(1,1))
cvssgblup2

op$init$scale=0.2568498
op$init$vare=0.3741784
op$ssGBLUPvar="hetVAR"
cvssgblup1<-baFit(driploss~sex,data=PigPheno,geno=genoNew ,genoid = ~id,ped=PigPed,options = op, train=~cv1)
cvssgblup1
par(mfrow=c(1,1))
baplot(cvssgblup1)
