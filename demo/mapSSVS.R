#This code demonstrate MAP SSVS model
rm(list=ls())
library(BATools)
data("Pig")

#Standardize genotype matrix
geno=std_geno(PigM,method="s",freq=PigAlleleFreq)
#use MCMC BayesA results as starting value
#or set your own starting values using 
#init=list(df=5,scale=0.01) 
demo(GBLUP)
init=set.init(~driploss,data=PigPheno,geno=geno,~id,
              df=5,scale=gblup$hyper_est[2],vare = gblup$hyper_est[1],g=gblup$ghat,
              beta=gblup$betahat,pi_snp=0.001,post_prob = NULL,h2=0.5,c=1000,model="SSVS",centered=TRUE,from="GBLUP")
run_para=list(maxiter=100)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-set.options(model="SSVS",method="MAP",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="mapSSVS",print_mcmc=NULL)
mapSSVS<-baFit(driploss~sex,data=PigPheno,geno=geno,genoid = ~id,options = op,map = PigMap,GWA="Win")
mapSSVS
par(mfrow=c(1,2))
man_plot_pvalue(mapSSVS)
man_plot_pvalue(mapSSVS,type="Win")

#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(~driploss,data = PigPheno,k=5)
cvmapSSVS<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1)
cvmapSSVS
par(mfrow=c(1,1))
baplot(cvmapSSVS)




