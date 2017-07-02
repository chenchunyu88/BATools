#######################################################################################
#######################################################################################
#Bayes EM approach using animal effect model
aBayesE= function(op=NULL,y=NULL,Z=NULL,X=NULL,vtrain=NULL,GWA=NULL,map=NULL) 
#  startpi is defined by the proportion of markers that ARE associate with genes
{   # START OF FUNCTION
	
  pi_math = 3.14159265359	
  set.seed(op$seed)  
  whichNa=which(vtrain==FALSE)
  Z0=Z
  y0=y
  X0=X
  if(length(whichNa)>0)
  {
    y=y[-whichNa]
    Z=Z[-whichNa,]
    X=X[-whichNa,]
  }
  
	if(GWA=="Win") {
	  win=c()
	  for(j in 1:max(map$idw)){
	    win=append(win,sum(map$idw==j))
	  }
	  chr=map %>% distinct(chr,idw) %>% select(chr) %>% t %>% as.numeric
	}else{
	  win=NULL
	  chr=NULL
	}
  
  
  nx=rownames(X)
  ng=rownames(Z)
  np=names(y)
  idx <- Reduce(intersect, list(nx,ng,np))
  X=X[idx,]
  y=y[idx]
  Z=Z[idx,]
  
  dimX=dim(X)[2]
  nSNP=dim(Z)[2]  #number of SNP
  nanim=dim(Z)[1]; #number of animals
  nrecords = length(y)
  
	
  	
	def=op$init$df
	if(is.null(op$init$vare)) {vare=as.numeric(var(y,na.rm=TRUE))} else {vare=op$init$vare}
	if(!is.null(op$init$scale)) scalea  = op$init$scale else{
		scalea=vare/100	
	}
	lambda=vare/scalea
	if(is.null(op$init$g)) SNPeff = rep(0,nSNP) else SNPeff=op$init$g
	if(is.null(op$init$beta)) fixedeff = rep(0,dimX) else fixedeff=op$init$beta 
	if(is.null(op$init$post_prob)) phi_est = rep(0,nSNP) else phi_est=op$init$post_prob
	if(!is.null(op$init$pi)) pi_snp=op$init$pi
	if(!is.null(op$init$c)) c=op$init$c
	
	nu_e=op$priors$nu_e
	tau2_e=op$priors$tau2_e
	nu_s=op$prior$nu_s
	tau2_s=op$prior$tau2_s
	
	maxiter=op$run_para$maxiter
	seed=op$seed
	set.seed(seed)
	
	Za=diag(1,nanim)
	Wa = cbind(X,Za)
	#WWdiag = apply(W,2,crossprod) # diagonals are important
	dimW = dim(Wa)[2]  # number of columns = # of location parameters

	# initialize a few arrays
	theta       = array(0,dimW)
	# so this is the residual.
	ycorr = y - as.vector(Wa%*%theta)   # adjust obs by everything before starting iterations

	XX = crossprod(X)
	XZa=  crossprod(X,Za)
	ZaX=  t(XZa)
	ZZa=  crossprod(Za)
	Wya=crossprod(Wa,y)
	thetakeep_a = array(0,ncol(Wa))
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
	##########start BayesEM###########
	while (abs(convcurr) > convcrit)
	{
  	  	SNPeff0=SNPeff
	  	iter = iter+1
		if(op$model=="SSVS")
	  	{
			h1=dnorm(SNPeff,mean=0,sd=sqrt(scalea))
			h0=dnorm(SNPeff,mean=0,sd=sqrt(scalea/c))
			if(iter!=1 && is.null(phi_est)) phi_est=pi_snp/((h0/h1)*(1-pi_snp)+pi_snp)
			Dinv=as.numeric((1-phi_est)*c+phi_est)
	  	}
		
		if(op$model=="BayesA")
	  	{
		    if(op$D=="V") Dinv=(as.numeric((def+1)/(def + SNPeff*SNPeff/scalea)))
		  	else Dinv=(as.numeric((def-1)/(def + SNPeff*SNPeff/scalea)))
	  	}
		
		if(op$model %in% c("rrBLUP","GBLUP")) Dinv=rep(1,length(Z[1,]))

		D=as.numeric(1/Dinv)
		
		ZD=.Call("BATools_set_ZD",Z,D)
		A=tcrossprod(ZD,Z)
		Ainv=solve(A)
		ZZ_G=ZZa+Ainv*as.numeric(lambda)
		coeff=rbind( cbind(XX,XZa),
		             cbind(ZaX,ZZ_G))
	
					 
		if(op$update_para$scale || op$update_para$vare){
		  	C=solve(coeff)
		  	theta=C%*%Wya
		  	a=theta[-(1:dimX)]
		  	ycorr=y-Wa%*%theta
		
			if(op$model %in% c("rrBLUP","GBLUP") && iter==1)
		  	{
		      vare=as.numeric(crossprod(y,ycorr)/(nrecords-rankX))
		      scalea=vare/lambda   
		  	}
		
	  		if(op$update_para$pi){
				  alpha_pi=1
		   		beta_pi=9
		   		pi_snp=(sum(phi_est)+alpha_pi-1)/(alpha_pi+beta_pi+nSNP-2)
			}
		  	Caa=C[(dimX+1):(dimX+nanim),(dimX+1):(dimX+nanim)]*vare
		
		  	fsigma2e=ycorr/vare
		  	WCW=Wa%*%C%*%t(Wa)
		  	Pfsigma2e=(fsigma2e-WCW%*%fsigma2e)/vare

		  	fsigma2u=Za%*%a/scalea
		  	Pfsigma2u=(fsigma2u-WCW%*%fsigma2u)/vare

		  	informAI[1,1]=t(fsigma2e)%*%Pfsigma2e+nu_e*tau2_e/(vare^3)-(nu_e+2)/(2*vare^2)
		  	informAI[1,2]=t(fsigma2e)%*%Pfsigma2u
		  	informAI[2,1]=t(fsigma2u)%*%Pfsigma2e
		  	informAI[2,2]=t(fsigma2u)%*%Pfsigma2u+nu_s*tau2_s/(scalea^3)-(nu_s+2)/(2*scalea^2)
		  	informAI=informAI/2
	  	
			traceCaa=sum(rowSums(Ainv*Caa))
		
			derivAI[1]=-0.5*((nanim-rankX)/vare-(nanim-traceCaa/scalea)/vare-crossprod(ycorr)/(vare^2))+nu_e*tau2_e/(2*vare^2)-(nu_e+2)/(2*vare)
			derivAI[2]=-0.5*(nanim/scalea-traceCaa/(scalea^2)-t(a)%*%Ainv%*%a/(scalea^2)) +nu_s*tau2_s/(2*scalea^2)-(nu_s+2)/(2*scalea)
	

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
			theta=solve(coeff,Wya)
			a=theta[-(1:dimX)]
		}
		  SNPeff=D*(t(Z)%*%(Ainv%*%a))
	  	if(iter%%4==0){
			scalea=as.numeric(scalea-(scalea-tscale[iter-1])^2/(scalea-2*tscale[iter-1]-tscale[iter-2]))
			vare=as.numeric(vare-(vare-tvare[iter-1])^2/(vare-2*tvare[iter-1]-tvare[iter-2]))
			if(op$model=="SSVS") pi_snp=as.numeric(pi_snp-(pi_snp-tpi[iter-1])^2/(pi_snp-2*tpi[iter-1]-tpi[iter-2]))
		}

	  	lambda = as.numeric(vare/scalea);
	  	#gamma = 1/lambda;

		if(op$model=="SSVS") cat("SSVS MAP iter=",iter,"\n")
		if(op$model=="BayesA") cat("BayesA MAP iter=",iter,"\n")
		if(op$model %in% c("rrBLUP","GBLUP"))  cat("GBLUP iter=",iter,"\n")
		cat ("Residual Variance is ",vare,sep="")
		cat (" Scale is ",scalea,sep="")
		if(op$model=="BayesC") cat (" pi is ",pi_snp,sep="")
		cat ("\n")
		
		

	  	tscale[iter]=scalea
	  	tvare[iter]=vare
		if(op$model=="SSVS") tpi[iter]=pi_snp
	  	thetakeep_a = cbind(thetakeep_a,theta)            #keep iterate
	  	if(iter>1) {
			if(op$update_para$scale || op$update_para$vare) convcurr=sqrt(sum(vardiff^2)/(vare^2+scalea^2))
			else convcurr=crossprod(SNPeff-SNPeff0)/crossprod(SNPeff)
		}
	  	tcon[iter]=convcurr
		cat("Convergence criteria is ",convcrit," and current value is ",convcurr,"\n",sep="")

	}
	if(op$model=="SSVS") cat("\nSSVS converged after ",iter," iterations at ",convcurr,"\n",sep="")
	if(op$model=="BayesA") cat("\nBayesA converged after ",iter," iterations at ",convcurr,"\n",sep="")
	if(op$model%in%c("rrBLUP","GBLUP")) cat("\n GBLUP converged after ",iter," iterations at ",convcurr,"\n",sep="")
	
	
	yhat=X0%*%theta[1:dimX]+Z0%*%SNPeff
	
	betahat=theta[1:dimX]
	if(op$update_para$scale || op$update_para$vare) sdbeta=sqrt(diag(C[1:dimX,1:dimX]))
		else {Cgg=NULL;sdbeta=NULL}
	
	names(betahat)=colnames(X)
	

	if(op$model=="SSVS"){
		hyper_est=c(vare,scalea,pi_snp)
		names(hyper_est)=c("vare","scale","pi")
	}
	if(op$model=="BayesA"){
		hyper_est=c(vare,scalea)
		names(hyper_est)=c("vare","scale")
	}
	if(op$model %in% c("rrBLUP","GBLUP")){
		hyper_est=c(vare,scalea)
		names(hyper_est)=c("vare","scale")
	}
	
	if(GWA!="No"){
		cat("Start calculate GWA results\n")
		
		if(op$model=="BayesA"){
		  meanvarg=(SNPeff^2+def*scalea)/(def+1)
		  Dinv=as.numeric(1/meanvarg*(1-2*SNPeff^2/(def+1)/meanvarg))
		}
	
		if(op$model=="SSVS"){
		  h1=dnorm(SNPeff,mean=0,sd=sqrt(scalea))
		  h0=dnorm(SNPeff,mean=0,sd=sqrt(scalea/c))
		  tau=as.numeric(pi_snp/((h0/h1)*(1-pi_snp)+pi_snp))
		  Dinv=as.numeric((tau+c*(1-tau))/scalea-SNPeff^2*tau*(1-tau)*(1-c)^2/scalea^2)
	  
		}
		if(op$model %in% c("rrBLUP","GBLUP")) Dinv=rep(1,length(Z[1,]))
		
		D=1/Dinv
  
		ZD=set_ZD(Z,D)
	  	A=tcrossprod(ZD,Z)
		Ainv=tryCatch({
			solve(A)
		}, error=function(e){
			solve(A+diag(0.0001,nanim))
		})

	  	if(op$model %in% c("rrBLUP","GBLUP")) ZZ_G=ZZa+Ainv*lambda else ZZ_G=ZZa+Ainv*vare
				
	 	coeff=rbind(cbind(XX,XZa),
	              cbind(ZaX,ZZ_G))

	    C=solve(coeff)
  
	    Caa=C[(dimX+1):(dimX+nanim),(dimX+1):(dimX+nanim)]*vare
  
	  	if(op$model %in% c("rrBLUP","GBLUP")) varahat=A*as.numeric(scalea)-Caa else varahat=A-Caa
		
		AiVAi=Ainv%*%varahat%*%Ainv

		AiVAiZD=AiVAi%*%ZD
		
		varg=rep(0,nSNP)
		
		for(i in 1:nSNP) varg[i]=sum(ZD[,i]*AiVAiZD[,i])
			
		if(op$model %in% c("rrBLUP","GBLUP")){
			zscore=SNPeff/sqrt(abs(varg))
			pvalue=2*pnorm(-abs(zscore))	
			Cgg=scalea-varg	
		}else{
			Cgg=D-varg
			zscore=SNPeff/sqrt(abs(Cgg))
			pvalue=2*pnorm(-abs(zscore)) 
		}
		
		if(GWA=="Win"){
			
			calc.win.r<-function(win){
				Cha=rep(0,length(win))
	
				current=1

				for(i in 1:length(win))
				{
					gw=SNPeff[current:(current+win[i]-1)]
					if(win[i]>1){
						t1=Z[,current:(current+win[i]-1)]
						t2=AiVAiZD[,current:(current+win[i]-1)]
						Cw=-crossprod(t1,t2) #mat_mut(t1,t2) #
						diag(Cw)=Cgg[current:(current+win[i]-1)]
		
						CwInv=tryCatch({
								solve(Cw)
							}, error=function(e){
								solve(Cw+diag(diag(Cw)*0.00000001))
						})
						Cha[i]=gw%*%CwInv%*%gw
					}else{
						Cha[i]=gw*gw/Cgg[current]
					}
					current=current+win[i]
				}
				Cha
			}

			calc.win.f<-function(win){
	
				Chf_RR=rep(0,length(win))
				current=1

				for(i in 1:length(win))
				{
					gw=SNPeff[current:(current+win[i]-1)]
					if(win[i]>1){
						t1=Z[,current:(current+win[i]-1)]
						t2=AiVAiZD[,current:(current+win[i]-1)]
						Cw=crossprod(t1,t2) #mat_mut(t1,t2) #


						diag(Cw)=varg[current:(current+win[i]-1)]
						CwInv=tryCatch({
								solve(Cw)
							}, error=function(e){
								solve(Cw+diag(diag(Cw)*0.00000001))
						})
						Chf_RR[i]=gw%*%CwInv%*%gw

					}else{
						Chf_RR[i]=gw*gw/varg[current]
			
					}
					current=current+win[i]
				}
				Chf_RR
			}
			
			if(op$model %in% c("rrBLUP","GBLUP")){
				wCh=calc.win.f(win)
				Wpvalue=pchisq(abs(wCh), win, lower.tail = FALSE)
			}else{
				wCh=calc.win.r(win)
				Wpvalue=pchisq(abs(wCh), win, lower.tail = FALSE)
			}
		}else{
			wCh=NULL
			Wpvalue=NULL
		}
		
	}else{
		zscore=NULL
		pvalue=NULL
		Cgg=NULL
		varg=NULL
		wCh=NULL
		Wpvalue=NULL
	}
	save(SNPeff,iter,tscale,tvare,tpi,thetakeep_a,varg,Cgg,theta,file = paste(op$save.at,op$seed,".RData",sep=""))
	
   BAout<-list(betahat=betahat,ghat=SNPeff, yhat=yhat,y=y0,train=vtrain,
	   hyper_est=hyper_est,Cgg=Cgg,varg=varg,pi_snp=pi_snp,phi_est=phi_est,
	   idx=idx,iter=iter,sdbeta=sdbeta,model=op$model,df=def,GWA=GWA,win=win,pvalue=pvalue,
	   zscore=zscore,Wpvalue=Wpvalue,wCh=wCh,map=map,Wchr=chr)
  
  	class(BAout)="ba"
  	return(BAout)	
}  # END OF FUNCTION

#######################################################################################
#######################################################################################

