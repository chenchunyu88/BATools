#This code demonstrate MAP SSVS model
rm(list=ls())
devtools::install("../BATools")
library(BATools)
data("Pig")

#Standardize genotype matrix
geno=std_geno(PigM,method="s",freq=PigAlleleFreq)
#use MCMC BayesA results as starting value
#or set your own starting values using 
#init=list(df=5,scale=0.01) 
demo(GBLUP)

init=set_init("driploss",data=PigPheno,geno=geno,"id",
              df=5,scale=gblup$hyper_est[2],vare = gblup$hyper_est[1],g=gblup$ghat,
              beta=gblup$betahat,pi_snp=0.001,h2=0.5,c=1000,model="SSVS",centered=TRUE,from="GBLUP")
run_para=list(maxiter=100)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-create.options(model="SSVS",method="MAP",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="mapSSVS",print_mcmc=NULL,D="P")
mapSSVS<-baTest(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op)
mapSSVS
#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(data = PigPheno,k=5,"driploss")
cvmapSSVS<-baTest(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1)
cvmapSSVS
plot(cvmapSSVS)




