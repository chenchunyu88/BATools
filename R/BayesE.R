#######################################################################################
#######################################################################################
#dataobj=pig;op=op;y=NULL;Z=NULL;G=NULL;X=NULL;trait="driploss"
#Bayes EM approach
BayesE = function(dataobj=NULL,op=NULL,y=NULL,Z=NULL,X=NULL,trait=NULL)  
#  startpi is defined by the proportion of markers that ARE associate with genes
{   # START OF FUNCTION
	
	set.seed(op$seed)  
	# Parameters
	if(!is.null(dataobj))
	{
		Z=dataobj$geno
		if(is.null(trait)){
			y=na.omit(dataobj$pheno[,1,1])	
		}else{
			y=na.omit(dataobj$pheno[,trait,1])
		}	
	}

	whichNa=which(is.na(y))
	if(length(whichNa)>0)
	{
		Z0=Z
		y0=y
		y=y[-whichNa]
		Z=Z[-whichNa,]
	}

  	if(op$poly)
	{
		Ainv<-dataobj$Ainv
		Zu<-diag(rep(1,nanim))
	  	pU<-dim(Zu)[1]
		Zu2<-apply(Zu,2,crossprod) #colSums(Zu*Zu)
		ueff<-array(0,pU)
		names(ueff)<-rownames(Ainv)
		varu<-op$init$varu
	}

  # input data

	if(!is.null(dataobj))
	{
	  	X = dataobj$fixed
		dimX = dim(X)[2]
		if(length(whichNa)>0)
		{
			X0=X
			X=X[-whichNa,]
		}
	}else{
		if(is.null(X)){
			X = matrix(1,length(y),1)
			dimX = 1
		}else dimX=dim(X)[2]
	}

    if(!is.null(dataobj)){
	    nx=rownames(X)
	    ng=rownames(Z)
	    np=names(y)
	    idx <- Reduce(intersect, list(nx,ng,np))
	    X=as.matrix(X[idx,],ncol=dimX)
	    y=y[idx]
	    Z=Z[idx,]
	}else idx=NULL
	nrecords = length(y)
    nSNP=dim(Z)[2]  #number of SNP
    nanim=dim(Z)[1]; #number of animals
	
  	
	def=op$init$df
	if(is.null(op$init$vare)) {vare=as.numeric(var(y,na.rm=TRUE))} else {vare=op$init$vare}
	if(!is.null(op$init$scale)) scalea  = op$init$scale else{
		#h2=0.5
	    #sumMeanZSq = sum((apply(Z,2L,mean))^2)
		#z2=apply(Z,2L,function(x) sum(x^2))
		#MSz=z2/nanim-sumMeanZSq
		#scalea=var(y,na.rm=TRUE)*h2/MSz	
		#try to use de los campos paper as starting values
		scalea=vare/100	
	}
	lambda=vare/scalea
	if(is.null(op$init$g)) SNPeff = rep(0,nSNP) else SNPeff=op$init$g
	if(is.null(op$init$b)) fixedeff = rep(0,dimX) else fixedeff=op$init$b 
	if(is.null(op$init$phi_est)) phi_est = rep(0,nSNP) else phi_est=op$init$phi_est
	if(!is.null(op$init$pi)) pi_snp=op$init$pi
	if(!is.null(op$init$c)) c=op$init$c
	
	nu_e=op$priors$nu_e
	tau2_e=op$priors$tau2_e
	nu_s=op$prior$nu_s
	tau2_s=op$prior$tau2_s
	
	maxiter=op$run_para$maxiter
	seed=op$seed
	set.seed(seed)

	W = cbind(X,Z)
	#WWdiag = apply(W,2,crossprod) # diagonals are important
	dimW = dim(W)[2]  # number of columns = # of location parameters

	# initialize a few arrays
	theta       = array(0,dimW)
	# so this is the residual.
	ycorr = y - as.vector(W%*%theta)   # adjust obs by everything before starting iterations

	XX = crossprod(X)
	XZ=  crossprod(X,Z)
	ZX=  crossprod(Z,X)
	ZZ=  crossprod(Z)
	Wy=crossprod(W,y)
	thetakeep = array(0,ncol(W))
	convcrit = op$convcrit
	convcurr = 1E10

	tscale=rep(0,2)
	tvare=rep(0,2)
	tpi=rep(0,2)
	tcon=rep(0,2)
	rankX=as.numeric(rankMatrix(X))
	derivAI=matrix(0,2,1)
	informAI=matrix(0,2,2)
	iter=0
	#vare=rb$Ve
	#scalea=rb$Vu
	#lambda=1000000
	##########start BayesEM###########
	while (abs(convcurr) > convcrit)
	{
  	  	SNPeff0=SNPeff
	  	iter = iter+1
		if(op$model=="SSVS")
	  	{
			h1=dnorm(SNPeff,mean=0,sd=sqrt(scalea))
			h0=dnorm(SNPeff,mean=0,sd=sqrt(scalea/c))
			phi_est=pi_snp/((h0/h1)*(1-pi_snp)+pi_snp)
			Dinv=((1-phi_est)*c+phi_est)
	  	}
		
		if(op$model=="BayesA")
	  	{
		    if(op$D=="V") Dinv=(as.numeric((def+1)/(def + SNPeff*SNPeff/scalea)))
		  	else Dinv=(as.numeric((def-1)/(def + SNPeff*SNPeff/scalea)))
	  	}
		
		if(op$model=="rrBLUP") Dinv=(length(Z[1,]))

		
		ZZ_G=ZZ
		diag(ZZ_G)=diag(ZZ)+Dinv*as.numeric(lambda)
	  	coeff=rbind( cbind(XX,XZ),
	               	 cbind(ZX,ZZ_G))
					 
		if(op$update_para$scale || op$update_para$vare){
		  	C=solve(coeff)
		  	theta=C%*%Wy
		  	SNPeff=theta[-(1:dimX)]
		  	ycorr=y-W%*%theta
		
			if(op$model=="rrBLUP" && iter==1)
		  	{
		      vare=as.numeric(crossprod(y,ycorr)/(nrecords-rankX))
		      scalea=vare/lambda   
		  	}
		
	  	   	if(op$update_para$pi){
				alpha_pi=1
		   		beta_pi=9
		   		pi_snp=(sum(phi_est)+alpha_pi-1)/(alpha_pi+beta_pi+nSNP-2)
			}
		  	Cgg=C[(dimX+1):(dimX+nSNP),(dimX+1):(dimX+nSNP)]*vare
		
		  	fsigma2e=ycorr/vare
		  	WCW=W%*%C%*%t(W)
		  	Pfsigma2e=(fsigma2e-WCW%*%fsigma2e)/vare

		  	fsigma2u=Z%*%SNPeff/scalea
		  	Pfsigma2u=(fsigma2u-WCW%*%fsigma2u)/vare

		  	informAI[1,1]=t(fsigma2e)%*%Pfsigma2e+nu_e*tau2_e/(vare^3)-(nu_e+2)/(2*vare^2)
		  	informAI[1,2]=t(fsigma2e)%*%Pfsigma2u
		  	informAI[2,1]=t(fsigma2u)%*%Pfsigma2e
		  	informAI[2,2]=t(fsigma2u)%*%Pfsigma2u+nu_s*tau2_s/(scalea^3)-(nu_s+2)/(2*scalea^2)
		  	informAI=informAI/2
	  	
			if(op$model=="rrBLUP"){
			    traceCgg=sum(diag(Cgg))	
			    derivAI[1]=-0.5*((nrecords-rankX)/vare-(nSNP-traceCgg/scalea)/vare-crossprod(ycorr)/(vare^2))+nu_e*tau2_e/(2*vare^2)-(nu_e+2)/(2*vare)
			    derivAI[2]=-0.5*(nSNP/scalea-traceCgg/(scalea^2)-crossprod(SNPeff)/(scalea^2))+nu_s*tau2_s/(2*scalea^2)-(nu_s+2)/(2*scalea)
			}else{
				traceCgg=sum(Dinv%*%diag(Cgg))
				derivAI[1]=-0.5*((nrecords-rankX)/vare-(nSNP-traceCgg/scalea)/vare-crossprod(ycorr)/(vare^2)) +nu_e*tau2_e/(2*vare^2)-(nu_e+2)/(2*vare)
				derivAI[2]=-0.5*(nSNP/scalea-traceCgg/(scalea^2)-t(SNPeff)%*%Dinv%*%SNPeff/(scalea^2)) +nu_s*tau2_s/(2*scalea^2)-(nu_s+2)/(2*scalea)
			}
		
		

		  	vardiff=solve(informAI)%*%derivAI
		  	if(op$update_para$vare){
			    if((vardiff[1]+vare)<=0)
				{
					vare=vare/2
				}else{
					vare=vare+as.numeric(vardiff[1])
				}
		  	}
			if(op$update_para$scale)
			{
			  	if((scalea+vardiff[2])<=0)
			  	{	
			    	scalea=scalea/2
			  	}else{
			    	scalea=scalea+as.numeric(vardiff[2])
			  	}	
			}
		}else{
			theta=solve(coeff,Wy)
			SNPeff=theta[-(1:dimX)]
		}

	  	if(iter%%4==0){
			scalea=as.numeric(scalea-(scalea-tscale[iter-1])^2/(scalea-2*tscale[iter-1]-tscale[iter-2]))
			vare=as.numeric(vare-(vare-tvare[iter-1])^2/(vare-2*tvare[iter-1]-tvare[iter-2]))
			if(op$model=="SSVS") pi_snp=as.numeric(pi_snp-(pi_snp-tpi[iter-1])^2/(pi_snp-2*tpi[iter-1]-tpi[iter-2]))
		}

	  	lambda = as.numeric(vare/scalea);
	  	#gamma = 1/lambda;

		if(op$model=="SSVS") cat("SSVS EM iter=",iter,"\n")
		if(op$model=="BayesA") cat("BayesA EM iter=",iter,"\n")
		if(op$model=="rrBLUP")  cat("rrBLUP iter=",iter,"\n")
		cat ("Residual Variance is ",vare,sep="")
		cat (" Scale is ",scalea,sep="")
		if(op$model=="BayesC") cat (" pi is ",pi_snp,sep="")
		cat ("\n")
		
		

	  	tscale[iter]=scalea
	  	tvare[iter]=vare
		if(op$model=="SSVS") tpi[iter]=pi_snp
	  	thetakeep = cbind(thetakeep,theta)            #keep iterate
	  	if(iter>1) {
			if(op$update_para$scale || op$update_para$vare) convcurr=sqrt(sum(vardiff^2)/(vare^2+scalea^2))
			else convcurr=crossprod(SNPeff-SNPeff0)/crossprod(SNPeff)
		}
	  	tcon[iter]=convcurr
		cat("Convergence criteria is ",convcrit," and current value is ",convcurr,"\n",sep="")

	}
	if(op$model=="SSVS") cat("\nSSVS converged after ",iter," iterations at ",convcurr,"\n",sep="")
	if(op$model=="BayesA") cat("\nBayesA converged after ",iter," iterations at ",convcurr,"\n",sep="")
	if(op$model=="rrBLUP") cat("\nrrBLUP converged after ",iter," iterations at ",convcurr,"\n",sep="")
	
	save(SNPeff,iter,tscale,tvare,tpi,thetakeep,theta,file = paste(op$save.at,op$seed,".RData",sep=""))
##########End BayesC###########

	yhat=X%*%theta[1:dimX]+Z%*%SNPeff
	
	betahat=theta[1:dimX]
	if(op$update_para$scale || op$update_para$vare) sdbeta=sqrt(diag(C[1:dimX,1:dimX]))
		else {Cgg=NULL;sdbeta=NULL}
	names(betahat)=colnames(X)
	
	
	if(op$model=="SSVS"){
		hyper_est=c(vare,scalea,pi_snp)
		names(hyper_est)=c("vare","varg","pi")
	}
	if(op$model=="BayesA"){
		hyper_est=c(vare,scalea)
		names(hyper_est)=c("vare","varg")
	}
	if(op$model=="rrBLUP"){
		hyper_est=c(vare,scalea)
		names(hyper_est)=c("vare","varg")
	}
	
  	#if(length(whichNa)>0)
  	#{
	#	Z=Z0
	#	if(!is.null(dataobj)) X=X0
  	#}
	
 	if(op$poly)
	{
		BAout<-list(betahat=betahat,ghat=SNPeff, yhat=theta[1]+Z%*%SNPeff+Zu%*%SNPeffBLUP_u,uhat=SNPeffBLUP_u,hypers=c(vare,scalea),Ginv=Ainv,Cuu=Cgg,sdbeta=sdbeta)
	}else BAout<-list(betahat=betahat,ghat=SNPeff, yhat=yhat,hyper_est=hyper_est,Cgg=Cgg,pi_snp=pi_snp,phi_est=phi_est,idx=idx,trait=trait,iter=iter,sdbeta=sdbeta)
  
  	class(BAout)="ba"
  	return(BAout)	
}  # END OF FUNCTION

#######################################################################################
#######################################################################################

