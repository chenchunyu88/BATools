#######################################################################################
#######################################################################################
#Bayes EM approach using animal effect model
ssBayesE= function(op=NULL,y=NULL,M=NULL,X=NULL,vtrain=NULL,GWA=NULL,map=NULL,Ainv=NULL) 
#  startpi is defined by the proportion of markers that ARE associate with genes
{   # START OF FUNCTION
	M2=M
	g_id=match(rownames(M2),rownames(Ainv))
	g_names=rownames(Ainv)[g_id]
	ng_names=rownames(Ainv)[-g_id]
	
	
	Ainv11=Ainv[-g_id,-g_id]
	Ainv12<-Ainv[-g_id,g_id]
	
	rhs=(-1.0)*Ainv12%*%M2
	M1hat=solve(Ainv11,rhs)
	
  pi_math = 3.14159265359	


  M=rbind(M1hat,M2)

  
  Z1=diag(1,dim(M1hat)[1])
  Z2=diag(1,dim(M2)[1])
  Z=diag(1,nrow(M))
  U=matrix(0,dim(Z)[1],dim(Z1)[2])
  U[1:dim(Z1)[2],1:dim(Z1)[2]]=Z1
  
  X1=X[ng_names,]
  X2=X[g_names,] 
  X=rbind(X1,X2)
  
  y1=y[ng_names]
  y2=y[g_names]
  y=c(y1,y2)

  vtrain1=vtrain[ng_names]
  vtrain2=vtrain[g_names]
  vtrain=c(vtrain1,vtrain2)

  set.seed(op$seed)  
  whichNa=which(vtrain==FALSE)
  M0=M
  y0=y
  X0=X
  Z0=Z
  M10=M1hat
  M20=M2
  U0=U
  if(length(whichNa)>0)
  {
    y=y[-whichNa]
    Z=Z[-whichNa,]
    X=X[-whichNa,]
    M<-M[-whichNa,]
    U<-U[-whichNa,]
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
  
  
  #nx=rownames(X)
  #ng=rownames(M)
  #np=names(y)
  #idx <- Reduce(intersect, list(nx,ng,np))
  #X=X[idx,]
  #y=y[idx]
  #Z=Z[idx,]
  
  dimX=dim(X)[2]
  nSNP=dim(M)[2]  #number of SNP
  nanim=dim(M)[1]; #number of animals
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
	
	W=cbind(X,Z,U)

	#WWdiag = apply(W,2,crossprod) # diagonals are important
	dimW = dim(W)[2]  # number of columns = # of location parameters

	# initialize a few arrays
	theta       = array(0,dimW)
	# so this is the residual.
	ycorr = y - as.vector(W%*%theta)   # adjust obs by everything before starting iterations
  

	Wy=crossprod(W,y)
	
	XX=crossprod(X)
	XZ=crossprod(X,Z)
	XU=crossprod(X,U)
	ZX=crossprod(Z,X)
	ZZ=crossprod(Z)
	MMstar=tcrossprod(M0)+diag(0.001,nrow=dim(M0)[1])
	if(!is.null(vtrain)) nanim=nrow(MMstar)
	ZU=crossprod(Z,U)
	UX=crossprod(U,X)
	UZ=crossprod(U,Z)
	UU=crossprod(U)
	rankX=as.numeric(rankMatrix(X))
	Ginv=solve(MMstar)
	
	vara=varg=scalea
	if(op$ssGBLUPvar=="hetVAR") varg=op$init$varGenetic 
	
	thetakeep = array(0,ncol(W))
	convcrit = op$convcrit
	convcurr = 1E10

	tvara=rep(0,2)
	tvarg=rep(0,2)
	tvare=rep(0,2)
	tcon=rep(0,2)
	tlog=rep(0,2)
	q1=nrow(M1hat)
	A11=Ainv11
	derivAI=matrix(0,2,1)
	if(op$ssGBLUPvar=="hetVAR") 	derivAI=matrix(0,3,1)
	informAI=matrix(0,2,2)
	if(op$ssGBLUPvar=="hetVAR") 	informAI=matrix(0,3,3)
	iter=0
	
	##########start BayesEM###########
	while (abs(convcurr) > convcrit && iter<op$run_para$maxiter)
	{
	 
  	  SNPeff0=SNPeff
	  	iter = iter+1

	  	r1=cbind(XX,XZ,XU)
	  	r2=cbind(ZX,ZZ+Ginv*vare/vara,ZU)
	    r3=cbind(UX,UZ,UU+A11*vare/varg) 
	  	coeff=rbind(r1,r2,r3)
	  


	
					 
		  if(op$update_para$scale || op$update_para$vare){
		    C=qr.solve(coeff)
		    theta=C%*%Wy
		    beta=theta[1:dimX]
		    g=theta[(dimX+1):(dimX+nanim)]
		    e=theta[-(1:(dimX+nanim))]
		    ycorr=y-W%*%theta
		
		  	Caa=C[(dimX+1):(dimX+nanim),(dimX+1):(dimX+nanim)]*vare
		  	if(q1==0) Cgg=matrix(0,0,0) else Cgg=C[(dimX+nanim+1):dim(C)[1],(dimX+nanim+1):dim(C)[1]]*vare
		  	
		
		  	WCW=W%*%C%*%t(W)
		  	fsigma2e=ycorr/as.numeric(vare)
		  	Pfsigma2e=(fsigma2e-W%*%C%*%t(W)%*%fsigma2e)/vare
		  	
		  	traceCgg=sum(diag(A11%*%Cgg))
		  	traceCaa=sum(diag(Ginv%*%Caa))
		  	if(op$ssGBLUPvar=="hetVAR") {

		  	  fsigma2a=Z%*%g/vara
		  	  Pfsigma2a=(fsigma2a-WCW%*%fsigma2a)/vare
		  	  
		  	  fsigma2g=U%*%e/varg
		  	  Pfsigma2g=(fsigma2g-WCW%*%fsigma2g)/vare
		  	  
		  	  informAI[1,1]=t(fsigma2e)%*%Pfsigma2e-1/(vare^2)
		  	  informAI[1,2]=t(fsigma2e)%*%Pfsigma2a
		  	  informAI[1,3]=t(fsigma2e)%*%Pfsigma2g  
		  	  informAI[2,1]=t(fsigma2a)%*%Pfsigma2e
		  	  informAI[2,2]=t(fsigma2a)%*%Pfsigma2a-1/(vara^2)
		  	  informAI[2,3]=t(fsigma2a)%*%Pfsigma2g
		  	  informAI[3,1]=t(fsigma2g)%*%Pfsigma2e
		  	  informAI[3,2]=t(fsigma2g)%*%Pfsigma2a
		  	  informAI[3,3]=t(fsigma2g)%*%Pfsigma2g-1/(varg^2) 
		  	  
		  	  derivAI[1]=-0.5*((nrecords-rankX)/vare-(nanim-traceCaa/vara+q1-traceCgg/varg)/vare-crossprod(ycorr)/(vare^2))-1/(vare*2)
		  	  derivAI[2]=-0.5*(nanim/vara-traceCaa/(vara^2)-crossprod(g,Ginv%*%g)/(vara^2))-1/(vara*2)  
		  	  derivAI[3]=-0.5*(q1/varg-traceCgg/(varg^2)-crossprod(e,A11%*%e)/(varg^2)) -1/(varg*2)
		  	  
		  	}else{
		  	  fsigma2a=Z%*%g/vara+U%*%e/vara
		  	  Pfsigma2a=(fsigma2a-W%*%C%*%t(W)%*%fsigma2a)/vare
		  	  
		  	  informAI[1,1]=t(fsigma2e)%*%Pfsigma2e+nu_e*tau2_e/(vare^3)-(nu_e+2)/(2*vare^2)
		  	  informAI[1,2]=t(fsigma2e)%*%Pfsigma2a
		  	  informAI[2,1]=t(fsigma2a)%*%Pfsigma2e
		  	  informAI[2,2]=t(fsigma2a)%*%Pfsigma2a+nu_s*tau2_s/(vara^3)-(nu_s+2)/(2*vara^2)
	
		  	  derivAI[1]=-0.5*((nrecords-rankX)/vare-(nanim-traceCaa/vara+q1-traceCgg/vara)/vare-crossprod(ycorr)/(vare^2))#-1/(vare*2)
		  	  derivAI[2]=-0.5*((nanim+q1)/vara-(traceCaa+traceCgg)/(vara^2)-(crossprod(g,Ginv%*%g)+crossprod(e,A11%*%e))/(vara^2))#-1/(vara*2)
		  	  
		  	}
		  	
		  
		  	informAI=informAI/2
	  	
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
			  	if((vara+vardiff[2])<=0)
			  	{	
			  	  vara=vara/2
			  	}else{
			    	vara=vara+as.numeric(vardiff[2])
			  	}	
			  
			    if(op$ssGBLUPvar!="hetVAR") varg=vara else{
			      if((varg+vardiff[3])<=0)
			      {	
			        varg=varg/2
			      }else{
			        varg=varg+as.numeric(vardiff[3])
			      }
			    }
			}
		  
		}else{
			theta=solve(coeff,Wy)
			a=theta[-(1:dimX)]
		}
		 
	  	if(iter%%4==0){
			    tmp=as.numeric(vara-(vara-tvara[iter-1])^2/(vara-2*tvara[iter-1]-tvara[iter-2]))
			    vara=ifelse(tmp>0,tmp,vara)
			    if(op$ssGBLUPvar!="hetVAR") varg=vara else {
			      tmp=as.numeric(varg-(varg-tvarg[iter-1])^2/(varg-2*tvarg[iter-1]-tvarg[iter-2]))
			      varg=ifelse(tmp>0,tmp,varg)
			    }
			    vare=as.numeric(vare-(vare-tvare[iter-1])^2/(vare-2*tvare[iter-1]-tvare[iter-2]))
			}

	  	#lambda = as.numeric(vare/vara);
	  	#gamma = 1/lambda;

		  
		  if(op$model %in% c("rrBLUP","ssGBLUP"))  cat("ssGBLUP iter=",iter,"\n")
		  cat ("Residual Variance is ",vare,sep="")
		  cat (", Marker variance is ",vara,sep="")
		  cat (", Genetic variance is ",varg,sep="")
		  cat ("\n")
		
		

	  	tvara[iter]=vara
	  	tvarg[iter]=varg
	  	tvare[iter]=vare
	
	  	thetakeep = cbind(thetakeep,theta)            #keep iterate
	  	if(iter>1) {
			if(op$update_para$scale || op$update_para$vare) convcurr=sqrt(sum(vardiff^2)/(vare^2+vara^2))
			else convcurr=crossprod(SNPeff-SNPeff0)/crossprod(SNPeff)
		}
	  	tcon[iter]=convcurr
		cat("Convergence criteria is ",convcrit," and current value is ",convcurr,"\n",sep="")

	}
	if(abs(convcurr) > convcrit) stop("\n ssGBLUP NOT converge after ",iter,
	                                  " iterations at ",convcurr,"\n",sep="") else{
	                                    cat("\n ssGBLUP  converged after ",iter," iterations at ",convcurr,"\n",sep="")
	                                    
	                                  }
	
	G2=tcrossprod(M2)
	names(g)=names(y0)
	SNPeff=crossprod(M2,solve(G2))%*%g[g_names]
	

	yhat=as.numeric(X0%*%theta[1:dimX]+M0%*%SNPeff+U0%*%e)
	names(yhat)=rownames(X0)
	betahat=theta[1:dimX]
	if(op$update_para$scale || op$update_para$vare) sdbeta=sqrt(diag(C[1:dimX,1:dimX])) else {Cgg=NULL;sdbeta=NULL}
	
	names(betahat)=colnames(X)
	

	
		hyper_est=c(vare,vara,varg)
		names(hyper_est)=c("vare","varMarker","varGenetics")
	
	
	if(GWA!="No"){
		cat("Start calculate GWA results\n")
		
	
	  gwas<-get_pvalue(U=g[g_names],V=c(vare,vara),G=G2,Z=M2,X=X2)
	  
	  est<-t(gwas[[1]])
	  colnames(est)=c("ghat","var_ghat")
    
	 
		
		varSNP=est[,2]
		
	
			

			zscore=est[,1]/sqrt(abs(varSNP))
			pvalue=2*pnorm(-abs(zscore))	
			Cgg=vara-varSNP	

		
		if(GWA=="Win"){
			
		  
		  calc.win<-function(win,AiVAiZ,Z,SNPeff,varg,Vu,M){
		    Chf_RR=rep(0,length(win))
		    current=1
		    Cgg=Vu-varg
		    for(i in 1:length(win))
		    {
		      idx=current:(current+win[i]-1)
		      
		      if(win[i]>1){
		        gw=SNPeff[idx]
		        t1=Z[,idx]
		        t2=AiVAiZ[,idx]
		        Cw=crossprod(t1,t2) #mat_mut(t1,t2) #
		        diag(Cw)=varg[idx]
		        CwInv=tryCatch({
		          solve(Cw)
		        }, error=function(e){
		          solve(Cw+diag(0.00001,nrow(Cw)))
		        })
		        Chf_RR[i]=gw%*%CwInv%*%gw
		        
		        if(abs(Chf_RR[i])>100){
		          subM=M[,idx]+1
		          if(length(idx)>2){
		            #subm=new("SnpMatrix",subM)
		            #LD <- ld(subm, depth = dim(subM)[2] - 1, stats = "R.squared", symmetric = T)
		            LD<-cor(subM)
		            diag(LD)=0
		            highLD_pairs=which(LD>0.99,arr.ind=TRUE)
		            if(length(highLD_pairs)!=0){
		              highLD_exclude=unique(apply(highLD_pairs,1,max))
		              if(length(highLD_exclude)>0) idx=idx[-highLD_exclude]
		            }
		          }else{
		            r  = snpgdsLDpair(subM[,1]-1,subM[,2]-1,method="r")[1]
		            if(r>.99) idx=idx[1]
		          }
		          if(length(idx)>1){
		            gw=SNPeff[idx]
		            t1=Z[,idx]
		            t2=AiVAiZ[,idx]
		            Cw=crossprod(t1,t2) #mat_mut(t1,t2) #
		            diag(Cw)=varg[idx]
		            CwInv=tryCatch({
		              solve(Cw)
		            }, error=function(e){
		              solve(Cw+diag(0.00001,nrow(Cw)))
		            })
		            Chf_RR[i]=gw%*%CwInv%*%gw
		            Cw=-Cw
		            diag(Cw)=Cgg[idx]
		          }else{
		            Chf_RR[i]=SNPeff[idx]*SNPeff[idx]/varg[idx]
		          }
		        }
		      }else{
		        gw=SNPeff[idx]
		        Chf_RR[i]=gw*gw/varg[current]
		      }
		      current=current+win[i]
		    }
		    Chf_RR
		  }
			
			
				wCh=calc.win(win,gwas[[2]],M2,gwas[[1]][1,],gwas[[1]][2,],vara,M2)
				Wpvalue=pchisq(abs(wCh), win, lower.tail = FALSE)
			
		}else{
			wCh=NULL
			Wpvalue=NULL
		}
		
	}else{
		zscore=NULL
		pvalue=NULL
		Cgg=NULL
		varSNP=NULL
		wCh=NULL
		Wpvalue=NULL
	}
	save(SNPeff,iter,tvarg,tvare,tvara,thetakeep,varSNP,Cgg,theta,vtrain,file = paste(op$save.at,op$seed,".RData",sep=""))
	
   BAout<-list(betahat=betahat,ghat=SNPeff, yhat=yhat,y=y0,y1,y2,train=vtrain,
	   hyper_est=hyper_est,Cgg=Cgg,varSNP=varSNP,pi_snp=pi_snp,phi_est=phi_est,
	   iter=iter,sdbeta=sdbeta,model=op$model,df=def,GWA=GWA,win=win,pvalue=pvalue,
	   zscore=zscore,Wpvalue=Wpvalue,wCh=wCh,map=map,Wchr=chr)
  
  	class(BAout)="ba"
  	return(BAout)	
}  # END OF FUNCTION

#######################################################################################
#######################################################################################

