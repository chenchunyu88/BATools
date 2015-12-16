#dataobj=pig;op=op;y=NULL;Z=NULL;yNa=NULL;trait="driploss"
BayesM = function(dataobj=NULL,op=NULL,y=NULL,Z=NULL,X=NULL,trait=NULL,yNa=NULL) 
#  startpi is defined by the proportion of markers that ARE associate with genes
{   # START OF FUNCTION
	pi_math = 3.14159265359	
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
		varusave=array(0,op$run_para$niter/op$run_para$skip)
		varu<-op$init$varu
		# This is for computing mean and variance of each SNP polygenenic effect.
		Mulast =    array(0,pU)
		Qulast =    array(0,pU)
		Mucurrent = array(0,pU)
		Qucurrent = array(0,pU)
	}
	# Vectors for saved samples of hyperparameters
	if (op$update_para$df)    {defsave   =  array(0,op$run_para$niter/op$run_para$skip) }
	if (op$update_para$scale)  {scalesave =  array(0,op$run_para$niter/op$run_para$skip) }
	if (op$update_para$pi)     {pisave    =  array(0,op$run_para$niter/op$run_para$skip) } 
	varesave = array(0,op$run_para$niter/op$run_para$skip)

	# To monitor MH acceptance ratios on degrees of freedom parameter
	alphadef_save = array(0,(op$run_para$niter-op$run_para$burnIn)/op$run_para$skip)
	# TO MONITOR MH ACCEPTANCE RATIOS ON SCALE PARAMETER (WHERE NECESSARY) POST BURN-IN
	alphascalea_save = array(0,(op$run_para$niter-op$run_para$burnIn)/op$run_para$skip)
	# input data
	
	if(!is.null(dataobj))
	{
	  	X = dataobj$fixed
		dimX= dim(X)[2]
	if(length(whichNa)>0)
	{
		X0=X
		X=X[-whichNa,]
	}
	}else{
	  X = matrix(1,length(y),1)
		dimX = 1
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
    nSNP=dim(Z)[2]  #number of SNP
    nanim=dim(Z)[1]; #number of animals
	nrecords = length(y)
  	G1=nSNP

  	W = as.matrix(cbind(X,Z))
  	Wstacked=as.vector(W);

	# specify prior degrees of belief on genetic variance
	# specify  prior scale on genetic variance
	# inital values   

	def=op$init$df			 # Starting Bayes A df 
	scalea  = op$init$scale    	 # Starting scale parameter
	Pi_SNP = op$init$pi                   # Starting pi parameter   (probability of being a gene!)
	cdef       = op$priors$cdef		 # cdef is scale parameter for proposal density on df
	cscalea     = op$priors$cscalea
	alphacheck = 0
	alphaQTLvar= 0			 # to monitor BayesB QTLvar rejection rates 
	alphatot1= 0
	if (op$update_para$scale & op$update_para$df) alphatot2= 0
	WWdiag     = apply(W,2,crossprod) # diagonals are important

	# used to monitor Metropolis Hastings acceptance ratio for df
	dimW     = dim(W)[2]  # number of columns = # of location parameters
	# initialize a few arrays
	theta      = array(0,dimW)
	SNPeff     = array(0,nSNP)
	if(op$model=="SSVS"){
	    phi=array(0,nSNP)
	    hratio=array(0,nSNP)
	    loghratio=array(0,nSNP)
	    phisave=array(0,nSNP)
	    hratiosave=c()
		  c=op$init$c
	}
	varq       = array(scalea,nSNP)  # initialize each SNP variance to overall scale parameter
	SNPvar     = array(0,nSNP)       # used to take running sum of SNP specific variances
	SNPvarlast = array(0,nSNP)
	postprob   = array(0,nSNP)       # indicator variable for inclusion of SNP in model
	                               # later used to determine posterior probability...unique to BayesB
	alphametrovar = array(0,nSNP)    # used to monitor MH acceptance rates for included genes.  

	# This is for computing mean and variance of each SNP effect.
	Mlast =    array(0,dimW)
	Qlast =    array(0,dimW)
	Mcurrent = array(0,dimW)
	Qcurrent = array(0,dimW)
	storeDbarj<-c()

	# adjust y
	ycorr = y - as.vector(W%*%theta)   # adjust obs by everything before starting iterations
	# so this is the residual.

	# use this to set up the starting value for the residual variance
	vare = crossprod(ycorr)/(nrecords-dimX)

	starttime<-proc.time()
	timeSNP = 0
	timetotal=0

	# mcmc sampling
	for (iter in 1:op$run_para$niter) 
	{  #"START" OF MCMC SAMPLING

     	itersave = iter-op$run_para$burnIn
     	if (itersave ==1) 
     	{
        	alphametrovar = array(0,nSNP)  # reinitialize to 0 post burn-in
        	if(op$model=="BayesB") postprob = array(0,nSNP) # reinitialize to 0 post burn-in
     	}

        #####################################################################################
        # sample intercept
        #####################################################################################
     	for(i in 1:dimX)
        {
   			ycorr    = ycorr + W[,i]*theta[i]
   			rhs=t(W[,i])%*%ycorr/vare# 1'(y-W*theta)/sigmae
   			lhs=WWdiag[i]/vare
   			invLhs=1.0/lhs		
           	meancond = rhs*invLhs                          
           	theta[i] = rnorm(1,meancond,sqrt(invLhs))       
           	ycorr    = ycorr - W[,i]*theta[i]  # take it off again to create new "residual"	
        }	

     	#####################################################################################	
     	# sample variance & effect for each SNP (BayesA)	
     	####################################################################################
		if (Pi_SNP == 1)  
		{ # "START" OF BAYESA
			#####################################################################################
			# sample variance for each SNP (Bayes A)	
			#####################################################################################			
			if(op$model=="BayesA") varq = (scalea*def + SNPeff*SNPeff)/rchisq(nSNP,def+1) else varq=rep(scalea,nSNP)
			#####################################################################################
			# sample effect for each SNP (Bayes A)	
			#####################################################################################
			startSNP = proc.time()
			BayesC<- .Call("BayesACL",nSNP,dimX,nanim,Wstacked,WWdiag,theta,ycorr,varq,vare)
			theta=BayesC[[1]]
			ycorr=BayesC[[2]]
			SNPeff=theta[-(1:dimX)]
			G1=nSNP
		}   # "END" FOR BAYESA

		#####################################################################################
		# sample effect for each SNP (Bayes B)	
		#####################################################################################
		if (Pi_SNP < 1 )  
		{   # "START" OF BAYESB
			startSNP = proc.time()
			if(op$model=="BayesB"){
				G1=0
				#####################################################################################
				# sample effect and variance for each SNP (Bayes B)	
				#####################################################################################
				BayesB<-.Call("BayesBC",nSNP,dimX,nanim,Wstacked,WWdiag,theta,ycorr,varq,vare,postprob,scalea,def,alphametrovar,iter,Pi_SNP,G1) 
				theta=BayesB[[1]]
				ycorr=BayesB[[2]]
				SNPeff=theta[-(1:dimX)]
				postprob=BayesB[[3]]
				if(itersave<=0) postprob   = array(0,nSNP)  
				alphametrovar=BayesB[[4]]
				G1=BayesB[[5]][1] #add G1
				alpha1=BayesB[[6]][1]
				alphacheck=BayesB[[7]][1]	
			}
			if(op$model=="SSVS"){
				#####################################################################################
				# sample effect and variance for each SNP (Bayes C)	
				#####################################################################################
		        BayesC<- .Call("BayesCC",nSNP,dimX,nanim,Wstacked,WWdiag, theta, ycorr,rep(scalea,nSNP),vare,Pi_SNP,phi,hratio,c)
		        theta=BayesC[[1]]
		        ycorr=BayesC[[2]]
		        SNPeff=theta[-(1:dimX)]
		        phi=BayesC[[3]]
		        if(itersave>0) phisave=phisave+phi
			}
	 		 	
		} # "END" OF BAYESB
		#sample polygenenic effects
		if(op$poly)
		{
			Avaru=Ainv/as.numeric(varu)
			startPoly=proc.time()
			res1<-.Call("sampleeff",pU,nanim,Zu,Zu2,ueff,ycorr,vare,Avaru)
			ueff=res1[[1]]
			ycorr=res1[[2]]
			endPoly=proc.time()
		    varu<-( t(ueff)%*%Ainv%*%ueff )/rchisq(1,length(ueff)-1) #sample sigmau
		}
		
		if (itersave>1) 
		{     # "START" OF COMPUTE POSTERIOR MEANS & VARIANCES
			Mcurrent = Mlast + (theta-Mlast)/itersave
			Qcurrent = Qlast + (itersave-1)*((theta-Mlast)^2)/itersave
			Mlast = Mcurrent
			Qlast = Qcurrent 

			SNPvar=SNPvarlast + (varq-SNPvarlast)/itersave
			SNPvarlast=SNPvar
			if(op$poly)
			{
				Mucurrent = Mulast + (ueff-Mulast)/itersave
				Qucurrent = Qulast + (itersave-1)*((ueff-Mulast)^2)/itersave
				Mulast =    Mucurrent
				Qulast =    Qucurrent	
			}
		}  # "END" OF COMPUTE POSTERIOR MEANS & VARIANCES
		if (itersave==1) {
			Mlast = theta;SNPvarlast=varq;
			if(op$poly)
			{
				Mulast=ueff;
			}
		} 

		#####################################################################################	
		#   sample residual variance
		#####################################################################################
		vare = ( crossprod(ycorr) +op$priors$tau2_e )/rchisq(1,nrecords+op$priors$nu_e)   # prior on vare used
		#  scale parameter should be e`e + Se,
		#  degrees of freedom should be n + ve

    	#####################################################################################
     	# sample degrees of freedom
     	#####################################################################################
     	if (op$update_para$df)  
		{  #"START" OF SAMPLE THE DEGREES OF FREEDOM
			SNPeff2=SNPeff*SNPeff
			for (cycle in 1:10){
		        mhdf=metropDf(def, scalea, G1=nSNP, SNPeff2, priordf = 1, cdef , max = Inf,pi_math=pi_math)
		        def=as.numeric(mhdf[[1]])
		        alpha1=as.numeric(mhdf[[2]])
		        alphatot1=alphatot1+alpha1
				if (cycle==10 & iter <= op$run_para$burnIn)
				{
				 	alphacheck1 = alphatot1/10;
				 	if (alphacheck1 > .70) { cdef = cdef*1.2} 
					else if(alphacheck1  < .20) {cdef = cdef*.7}     # tune Metropolis sampler #
					alphatot1 = 0;
				}
				else if (cycle==10 & iter > op$run_para$burnIn) 
				{
				 	alphacheck1 = alphatot1/10;
				 	alphadef_save[itersave/op$run_para$skip] = alphacheck1
				 	alphatot1 = 0
				}
			}
     	}   # "END" OF SAMPLE THE DEGREES OF FREEDOM  

        if(op$model=="SSVS"){
	        #####################################################################################
	        #  sample scale#
	        #  based on Gelman Uniform(0,A) prior on scale parameter for "large" A.
	        #####################################################################################
	        if (op$update_para$scale) 
	         {  #"START" OF SAMPLE THE SCALE PARAMETER 
	              #shape_g = 0.5*(nSNP+5)
	   	   	   	  tmp=SNPeff*SNPeff/((1-phi)/c+phi)
	              #scale_g = 0.5*(sum(tmp)+0.04)
	              scalea =  (sum(tmp))/rchisq(1,nSNP)
	   	    	  #rgamma(1,shape=shape_g,scale=scale_g)
	         } #"END" OF SAMPLE THE SCALE PARAMETER 
        }else if(op$model=="rrBLUP"){
	        if (op$update_para$scale) 
	         {  #"START" OF SAMPLE THE SCALE PARAMETER 
	              #shape_g = 0.5*(nSNP+5)
	   	   	   	  tmp=SNPeff*SNPeff
	              #scale_g = 0.5*(sum(tmp)+0.04)
	              scalea =  (sum(tmp))/rchisq(1,nSNP)
	   	    	  #rgamma(1,shape=shape_g,scale=scale_g)
	         } #"END" OF SAMPLE THE SCALE PARAMETER 
        }else if(op$model%in%c("BayesA","BayesB")){
	        #####################################################################################
	        #  sample scale#
	        #  based on Gelman Uniform(0,A) prior on scale parameter for "large" A.
	        #####################################################################################
	        if (op$update_para$scale & !(op$update_para$df)) 
	         {  #"START" OF SAMPLE THE SCALE PARAMETER 
	              shape_g = 0.5*(G1*def)+op$priors$shape_scale
	   		      scale_g = 0.5*def*sum(1/varq[which(varq>0)])+op$priors$rate_scale
	              scalea = rgamma(1,shape=shape_g,rate=scale_g)
	         } #"END" OF SAMPLE THE SCALE PARAMETER   
		     
	    	#####################################################################################
	     	#  sample scale usning UNIMH#
	     	#####################################################################################
			if (op$update_para$scale & op$update_para$df)       	  
			{  #"START" OF SAMPLE THE SCALE PARAMETER 
		 		SNPeff2=SNPeff*SNPeff
				for (cycle in 1:10){
					mhscale=metropScale(def, scalea, G1=G1, SNPeff2, cscalea,op$priors$shape_scale, op$priors$rate_scale, pi_math=pi_math)
					scalea=as.numeric(mhscale[[1]])
					alpha2=as.numeric(mhscale[[2]])
					alphatot2=alphatot2+alpha2
					if (cycle==10 & iter <= op$run_para$burnIn)
					{
						alphacheck2 = alphatot2/10;
						if (alphacheck2 > .70) { cscalea = cscalea*1.2}
						else if(alphacheck2  < .20) {cscalea = cscalea*.7}     # tune Metropolis sampler #
						alphatot2 = 0;
					}
					else if (cycle==10 & iter > op$run_para$burnIn) 
					{
						alphacheck2 = alphatot2/10;
						alphascalea_save[itersave/op$run_para$skip] = alphacheck2
						alphatot2 = 0
					}
				}
	      	} #"END" OF SAMPLE THE SCALE PARAMETER 
		}
		#####################################################################################
		#      sample pi
		#####################################################################################
		if (op$update_para$pi) 
		{ 
			if(op$model=="BayesB") 	Pi_SNP = rbeta(1,op$priors$alphapi+G1,op$priors$betapi+nSNP-G1) 
			if(op$model=="SSVS") {
		 	   m1=sum(phi)
		       Pi_SNP=rbeta(1,op$priors$alphapi+m1,op$priors$betapi+nSNP-m1)
			}	
		} #"END" OF SAMPLE THE PI PARAMETER 

		###########################################################################################
		#store MCMC samples at every "skip"
		###########################################################################################
		if (iter%%op$run_para$skip == 0) 
		{
			varesave[iter/op$run_para$skip] = vare
			if(op$update_para$df) {defsave[iter/op$run_para$skip] = def}
			if(op$update_para$scale) {scalesave[iter/op$run_para$skip] = scalea}
			if(op$update_para$pi) {pisave[iter/op$run_para$skip] = Pi_SNP}
			if(op$poly){varusave[iter/op$run_para$skip]=varu}
		}

		endSNP = proc.time()
		timeSNP = timeSNP + (endSNP[1]-startSNP[1]) 
		endtime = proc.time()
		timetotal = endtime[1]-starttime[1]
		
		if((op$print_mcmc$piter!=0) & iter%%op$print_mcmc$piter==0)
		{
			tmp=paste("iter= ",iter," vare= ",round(vare,6), "scale= ", round(scalea,8),"\ntimepercycle= ", round(timeSNP/iter,3), "estimated time left=", round(timeSNP/iter*(op$run_para$niter-iter),2),"seconds \n",sep=" ")
			if(op$print_mcmc$print_to=="screen")
			{
				cat(tmp)
			}else{
			write(tmp, file = paste("log",op$seed,".txt",sep=""),append = TRUE)
			}
		}
       if (iter%%1000==0) 
       {
			#############################################################################################
			#    PROCESSING MCMC SAMPLES
			#############################################################################################
			hyperparameters = as.matrix(varesave)
			colnames(hyperparameters) = "varesave"
			if (op$update_para$scale) {hyperparameters = cbind(hyperparameters,scalesave) }
			if (op$update_para$df) {hyperparameters = cbind(hyperparameters,defsave) }
			if (op$update_para$pi) {hyperparameters = cbind(hyperparameters,pisave) }
			if(op$poly) {hyperparameters = cbind(hyperparameters,varusave)}	
			rownames(hyperparameters) = seq(op$run_para$skip,op$run_para$niter,op$run_para$skip)
			save(iter,hyperparameters,file=paste(op$save.at, "Hyperparameters",op$seed,".RData",sep="")) 
			SNPtimepercycle = timeSNP/iter
			timepercycle=timetotal/iter
			save(iter,timepercycle,SNPtimepercycle,file=paste(op$save.at, "time",op$seed,".RData",sep=""))

       	}# "END" OF PROCESSING MCMC SAMPLING
               
		if (iter%%1000==0 & itersave>1 ) 
		{
		#############################################################################################
		#    PROCESSING MCMC SAMPLES
		#############################################################################################
			meanmu=Mcurrent[1:dimX]
			names(meanmu)=colnames(X)
			meang= Mcurrent[(dimX+1):dimW]
			names(meang)=colnames(Z)
			varg = Qcurrent[(dimX+1):dimW]/itersave
			names(varg)=colnames(Z)
			postprob_save   = postprob/itersave
			names(postprob_save)=colnames(Z)
 	        save(iter,alphadef_save,alphascalea_save,file=paste(op$save.at, "alphadef",op$seed,".RData",sep=""))
 			if(op$poly)
 			{
 				meanu<-Mucurrent
 				pvaru<-Qucurrent/itersave
 				if(op$model=="BayesB") save(iter,meanmu,meang,postprob_save,varg,SNPvar,meanu,pvaru,file = paste(op$save.at,"EffectsResults",op$seed,".RData",sep=""))
 					else save(iter,meanmu,meang,varg,SNPvar,meanu,pvaru,file = paste(op$save.at,"EffectsResults",op$seed,".RData",sep=""))
 			}else{
 				if(op$model=="BayesB") save(iter,meanmu,meang,postprob_save,varg,SNPvar,file = paste(op$save.at,"EffectsResults",op$seed,".RData",sep=""))
 					else save(iter,meanmu,meang,varg,SNPvar,file = paste(op$save.at,"EffectsResults",op$seed,".RData",sep=""))
 			}
		} # "END" OF PROCESSING MCMC SAMPLING

  	} # "END" OF MCMC SAMPLING
	#if(length(whichNa)>0)
	#{
	#Z=Z0
	#if(!is.null(dataobj)) X=X0
	#}
	if(op$model=="SSVS"){
	  phisave=phisave/itersave
	  names(phisave)=colnames(Z)
	}
	sample_idx=c((op$run_para$burnIn/op$run_para$skip+1):(op$run_para$niter/op$run_para$skip))
	hyper_est=apply(hyperparameters[sample_idx,],2,mean)
	names_hypers=c("vare")
	if (op$update_para$scale) {names_hypers = c(names_hypers,"scale") }
	if (op$update_para$df) {names_hypers = c(names_hypers,"df") }
	if (op$update_para$pi) {names_hypers = c(names_hypers,"pi") }
	if(op$poly) {hyperparameters = c(names_hypers,"varu")}	
	colnames(hyperparameters)=names_hypers
	names(hyper_est)=names_hypers  
	#if(op$model%in%c("BayesB","SSVS")) colnames(hyperparameters)=c("vare","varg","pi") else colnames(hyperparameters)=c("vare","varg")
	#if(op$model%in%c("BayesB","SSVS")) names(hyper_est)=c("vare","varg","pi") else names(hyper_est)=c("vare","varg")
  	if(op$poly)
  	{
		if(op$model=="BayesB") BAout<-list(betahat=meanmu,ghat=meang,yhat=X%*%meanmu+Z%*%meang+Zu%*%meanu,prob=postprob_save,uhat=meanu,eff_sample=effectiveSize(hyperparameters),hypers=hyperparameters,idx=idx,trait=trait,hyper_est=hyper_est)
		else if (op$model=="SSVS") BAout<-list(betahat=meanmu,ghat=meang, yhat=X%*%meanmu+Z%*%meang+Zu%*%meanu,uhat=meanu,eff_sample=effectiveSize(hyperparameters),hypers=hyperparameters,phisave=phisave,idx=idx,trait=trait,hyper_est=hyper_est)
		else BAout<-list(betahat=meanmu,ghat=meang, yhat=X%*%meanmu+Z%*%meang+Zu%*%meanu,uhat=meanu,eff_sample=effectiveSize(hyperparameters),hypers=hyperparameters,idx=idx,trait=trait,hyper_est=hyper_est)
	}else {
		if(op$model=="BayesB") BAout<-list(betahat=meanmu,ghat=meang, yhat=X%*%meanmu+Z%*%meang,prob=postprob_save,eff_sample=effectiveSize(hyperparameters),hypers=hyperparameters,idx=idx,trait=trait,hyper_est=hyper_est)
		else if (op$model=="SSVS") BAout<-list(betahat=meanmu,ghat=meang, yhat=X%*%meanmu+Z%*%meang,eff_sample=effectiveSize(hyperparameters),hypers=hyperparameters,phisave=phisave,idx=idx,trait=trait,hyper_est=hyper_est)
		else BAout<-list(betahat=meanmu,ghat=meang, yhat=X%*%meanmu+Z%*%meang,eff_sample=effectiveSize(hyperparameters),hypers=hyperparameters,idx=idx,trait=trait,hyper_est=hyper_est)
	}
	class(BAout)="ba"
   	return(BAout)
}  # END OF FUNCTION

#######################################################################################
#######################################################################################
