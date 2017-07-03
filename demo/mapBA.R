#This code demonstrate MAP BayesA model
rm(list=ls())

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
              beta=gblup$betahat,pi_snp=1,h2=0.5,c=NULL,model="BayesA",centered=TRUE,from="GBLUP")

run_para=list(maxiter=100)
update_para=list(df=FALSE,scale=TRUE,pi=FALSE)
op<-create.options(model="BayesA",method="MAP",priors=NULL,init=init,
                   update_para=update_para,run_para=run_para,save.at="mapBayesA",print_mcmc=NULL,D="P")

mapBA<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op,map = PigMap,GWA="Win")
mapBA
par(mfrow=c(1,2))
man_plot_pvalue(mapBA)
man_plot_pvalue(mapBA,type="Win")
#### Cross-validation using BATools
set.seed(1234)
PigPheno=createCV(data = PigPheno,k=5,"driploss")
cvmapBA<-baFit(driploss~sex,data=PigPheno,geno=geno ,genoid = ~id,options = op, train=~cv1)
cvmapBA
par(mfrow=c(1,1))
plot(cvmapBA)



