##calculate the variance explained by each window
calc.varw<-function(win,Z,SNPeff){
	len_win=length(win)
	varw=rep(0,len_win)
	start=1
	nanim=dim(Z)[1]
	for(k in 1:len_win){
		end=start+win[k]-1
		if(start==end) gw=Z[,start]*as.numeric(SNPeff[start]) else gw=Z[,start:end]%*%as.numeric(SNPeff[start:end])
		varw[k]=sum(gw*gw)/nanim-(sum(gw)/nanim)^2
		start=end+1
	}
	varw
}

##calculate the posterior probability of each window contains at least one non-zero effect for variable selection methods 
calc.phiw<-function(win,phi){
	start=1
	lenw=length(win)
	phiw=rep(0,lenw)
	for(k in 1:lenw)
	{
		end=start+win[k]-1
		phiw[k]=(sum(phi[start:end])>=1)
		if(k!=lenw) start=end+1
	}
	phiw
}


#Main function for Bayesian whole genome regression when all individuals are genotyped
ssBayesM <- function(op=NULL,y=NULL,M=NULL,X=NULL,vtrain=NULL,GWA=NULL,map=NULL,Ainv=NULL) 
{   # START OF FUNCTION
  M2=M
  g_id=match(rownames(M2),rownames(Ainv))
  g_names=rownames(Ainv)[g_id]
  ng_names=rownames(Ainv)[-g_id]
  
  Ainv=as.matrix(Ainv)
  
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
	
	if(dim(X)[2]==1) X1=as.matrix(X[ng_names,],ncol=1) else X1=X[ng_names,]
	if(dim(X)[2]==1) X2=as.matrix(X[g_names,],ncol=1) else X2=X[g_names,] 
	X=rbind(X1,X2)
	
	y1=y[ng_names]
	y2=y[g_names]
	y=c(y1,y2)
	
	nanim.ungeno=nrow(M1hat)
	
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
	
	# Vectors for saved samples of hyperparameters
	if (op$update_para$df)    {defsave   =  array(0,op$run_para$niter/op$run_para$skip) }
	 {scalesave =  array(0,op$run_para$niter/op$run_para$skip) }
	if (op$update_para$pi)     {pisave    =  array(0,op$run_para$niter/op$run_para$skip) } 
	varesave = array(0,op$run_para$niter/op$run_para$skip)

	# To monitor MH acceptance ratios on degrees of freedom parameter
	alphadef_save = array(0,(op$run_para$niter-op$run_para$burnIn)/op$run_para$skip)
	# TO MONITOR MH ACCEPTANCE RATIOS ON SCALE PARAMETER (WHERE NECESSARY) POST BURN-IN
	alphascalea_save = array(0,(op$run_para$niter-op$run_para$burnIn)/op$run_para$skip)
	# input data
	

	
	
	#nx=rownames(X)
	#ng=rownames(Z)
	#np=names(y)
	#idx <- Reduce(intersect, list(nx,ng,np))
	#X=X[idx,]
	#y=y[idx]
	#Z=Z[idx,]
	
	dimX=dim(X)[2]
  nSNP=dim(M)[2]  #number of SNP
  nanim=dim(M)[1]; #number of animals
	nrecords = length(y)
  G1=nSNP
  
  nanim.ungeno=dim(M1hat)[1]
  
  dimAinv=dim(Ainv11)[1]
  Mlast2    =  array(0,dimAinv);  Qlast2    = array(0,dimAinv);              # for save imputation residuals
  Mcurrent2 =  array(0,dimAinv);  Qcurrent2 = array(0,dimAinv);
  epsilon     = array(0,dimAinv)                                                  # epsilon is the imputation residuals
  sigma_g=op$init$varGenetic
  sigma_gsave    = array(0,op$run_para$niter/op$run_para$skip)

  W=as.matrix(cbind(X,M))
  
  Z1Z1diag=apply(Z1,2,crossprod)


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
	if(op$model=="ssSSVS"){
	    hratio=array(0,nSNP)
	    loghratio=array(0,nSNP)
	    hratiosave=c()
		  c=op$init$c
	}
	varq       = array(scalea,nSNP)  # initialize each SNP variance to overall scale parameter
	SNPvar     = array(0,nSNP)       # used to take running sum of SNP specific variances
	SNPvarlast = array(0,nSNP)
	postprob   = array(0,nSNP)       # indicator variable for inclusion of SNP in model
                                    # later used to determine posterior probability...unique to BayesB and SSVS
	postprob_save=array(0,nSNP)
	if(GWA=="Win" && op$model=="ssBayesB")
	{
	  p1=array(0,nSNP)
	  p2=array(0,nSNP)
	  pw=array(0,nSNP)
	}
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
	if(is.null(op$init$vare)) vare = crossprod(ycorr)/(nrecords-dimX) else vare=op$init$vare

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
        	if(op$model=="ssBayesB") postprob = array(0,nSNP) # reinitialize to 0 post burn-in
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
			if(op$model=="ssBayesA") varq = (scalea*def + SNPeff*SNPeff)/rchisq(nSNP,def+1) else varq=rep(scalea,nSNP)
			#####################################################################################
			# sample effect for each SNP (Bayes A)	
			#####################################################################################
			startSNP = proc.time()
			BayesC<- .Call("BayesACL",nSNP,dimX,nanim,W,WWdiag,theta,ycorr,varq,vare)
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
			if(op$model=="ssBayesB"){
				G1=0
				#####################################################################################
				# sample effect and variance for each SNP (Bayes B)	
				#####################################################################################
				BayesB<-.Call("BayesBC",nSNP,dimX,nanim,W,WWdiag,theta,ycorr,varq,vare,postprob,scalea,def,alphametrovar,iter,Pi_SNP,G1) 
				theta=BayesB[[1]]
				ycorr=BayesB[[2]]
				SNPeff=theta[-(1:dimX)]
				p1=BayesB[[3]]
				postprob=BayesB[[3]]
			  postprob_save=postprob_save+postprob 
				postprob   = array(0,nSNP) 
				alphametrovar=BayesB[[4]]
				G1=BayesB[[5]][1] #add G1
				alpha1=BayesB[[6]][1]
				alphacheck=BayesB[[7]][1]	
			}
			if(op$model=="ssSSVS"){
				#####################################################################################
				# sample effect and variance for each SNP (Bayes C)	
				#####################################################################################
		        BayesC<- .Call("BayesCC",nSNP,dimX,nanim,W,WWdiag, theta, ycorr,rep(scalea,nSNP),vare,Pi_SNP,postprob,hratio,c)
		        theta=BayesC[[1]]
		        ycorr=BayesC[[2]]
		        SNPeff=theta[-(1:dimX)]
		        postprob=BayesC[[3]]
		        if(itersave>0) postprob_save=postprob_save+postprob
			}
	 		 	
		} # "END" OF BAYESB/SSVS
    
     	####sample imputation like updating polygenenic effects
     	Avaru=Ainv11/as.numeric(sigma_g)  # A/varu easier to pass to c function
     	
     	res1<-.Call("sampleeff",dimAinv,nanim.ungeno,Z1,Z1Z1diag,epsilon,ycorr,vare,Avaru)
     	epsilon=res1[[1]]
     	ycorr=res1[[2]]
     	
     	sigma_g<-( t(epsilon)%*%Ainv11%*%epsilon +op$priors$def_sigmag0*op$priors$scale_sigmag0)/rchisq(1,nanim.ungeno+op$priors$def_sigmag0) #sample sigmau
     	
     	
		if (itersave>1) 
		{     # "START" OF COMPUTE POSTERIOR MEANS & VARIANCES
			Mcurrent = Mlast + (theta-Mlast)/itersave
			Qcurrent = Qlast + (itersave-1)*((theta-Mlast)^2)/itersave
			Mlast = Mcurrent
			Qlast = Qcurrent 
			
			#imputation residual
			Mcurrent2 = Mlast2 + (epsilon - Mlast2)/itersave;
			Qcurrent2 = Qlast2 + (itersave - 1)*((epsilon - Mlast2)^2)/itersave;
			Mlast2 = Mcurrent2; Qlast2 = Qcurrent2;
			

			SNPvar=SNPvarlast + (varq-SNPvarlast)/itersave
			SNPvarlast=SNPvar
			if (iter%%op$run_para$skip == 0) {
			  if(op$model %in% c("ssrrBLUP","ssBayesA")){
			    count=count+(SNPeff>0)
			    if(GWA=="Win"){
			      g_total<-M%*%SNPeff
			      var_total=sum(g_total*g_total)/nanim-(sum(g_total)/nanim)^2
			      varw=calc.varw(win,M,SNPeff)
			      qw=varw/var_total
			      countqw=countqw+(qw>0.01)
			    }
			  }
			 if(op$model %in% c("ssSSVS","ssBayesB") && GWA=="Win"){
			   if(op$model=="ssSSVS") Wpostprob=Wpostprob+calc.phiw(win,postprob) else Wpostprob=Wpostprob+calc.phiw(win,p1)
			   n_saved_sample=n_saved_sample+1
			   #cat("\n p2",p2,"\n Wpb",Wpostprob,"\n")
			  } 
			} 
		}  # "END" OF COMPUTE POSTERIOR MEANS & VARIANCES
		if (itersave==1) {
			Mlast = theta;SNPvarlast=varq;Mlast2 = epsilon
		
			  if(op$model %in% c("ssrrBLUP","ssBayesA")){
			    count=rep(0,nSNP)
			    if(GWA=="Win") countqw=rep(0,length(win))
			  }
			  if(op$model %in% c("ssSSVS","ssBayesB") && GWA=="Win"){
			    Wpostprob=rep(0,length(win))
			    n_saved_sample=0
			  } 
			
		} 

		#####################################################################################	
		#   sample residual variance
		#####################################################################################
     if(op$update_para$vare) vare = ( crossprod(ycorr) +op$priors$tau2_e )/rchisq(1,nrecords+op$priors$nu_e)   # prior on vare used
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

        if(op$model=="ssSSVS"){
	        #####################################################################################
	        #  sample scale#
	        #  based on Gelman Uniform(0,A) prior on scale parameter for "large" A.
	        #####################################################################################
	        if (op$update_para$scale) 
	         {  #"START" OF SAMPLE THE SCALE PARAMETER 
	              #shape_g = 0.5*(nSNP+5)
	   	   	   	  tmp=SNPeff*SNPeff/((1-postprob)/c+postprob)
	              #scale_g = 0.5*(sum(tmp)+0.04)
	              scalea =  (sum(tmp))/rchisq(1,nSNP)
	   	    	  #rgamma(1,shape=shape_g,scale=scale_g)
	         } #"END" OF SAMPLE THE SCALE PARAMETER 
        }else if(op$model=="ssrrBLUP"){
	        if (op$update_para$scale) 
	         {  #"START" OF SAMPLE THE SCALE PARAMETER 
	              #shape_g = 0.5*(nSNP+5)
	   	   	   	  tmp=SNPeff*SNPeff
	              #scale_g = 0.5*(sum(tmp)+0.04)
	              scalea =  (sum(tmp))/rchisq(1,nSNP)
	   	    	  #rgamma(1,shape=shape_g,scale=scale_g)
	         } #"END" OF SAMPLE THE SCALE PARAMETER 
        }else if(op$model%in%c("ssBayesA","ssBayesB")){
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
			if(op$model=="ssBayesB") 	Pi_SNP = rbeta(1,op$priors$alphapi+G1,op$priors$betapi+nSNP-G1) 
			if(op$model=="ssSSVS") {
		 	   m1=sum(postprob)
		       Pi_SNP=rbeta(1,op$priors$alphapi+m1,op$priors$betapi+nSNP-m1)
			}	
		} #"END" OF SAMPLE THE PI PARAMETER 

		###########################################################################################
		#store MCMC samples at every "skip"
		###########################################################################################
		if (iter%%op$run_para$skip == 0) 
		{
			varesave[iter/op$run_para$skip] = vare
			scalesave[iter/op$run_para$skip] = scalea
			sigma_gsave[iter/op$run_para$skip]=sigma_g
			
			if(op$update_para$df) {defsave[iter/op$run_para$skip] = def}
			
			if(op$update_para$pi) {pisave[iter/op$run_para$skip] = Pi_SNP}
			
		}

		endSNP = proc.time()
		timeSNP = timeSNP + (endSNP[1]-startSNP[1]) 
		endtime = proc.time()
		timetotal = endtime[1]-starttime[1]
		
		if((op$print_mcmc$piter!=0) & iter%%op$print_mcmc$piter==0)
		{
		  if(op$model %in% c("ssrrBLUP","ssBayesA")){
		    tmp=paste("iter= ",iter," vare= ",round(vare,6), " Marker scale= ", round(scalea,6),"Genetic Variance= ",round(sigma_g,8) ,
		              "\ntimepercycle= ", round(timeSNP/iter,3), "estimated time left=", 
		              round(timeSNP/iter*(op$run_para$niter-iter),2),"seconds \n",sep=" ")
		  }
		  if(op$model %in% c("ssSSVS","ssBayesB")){
		    tmp=paste("iter= ",iter," vare= ",round(vare,6), "Marker scale= ", round(scalea,6),"pi= ",round(Pi_SNP,6),"Genetic Variance= ",round(sigma_g,8) ,
		              "\ntimepercycle= ", round(timeSNP/iter,3), "estimated time left=", 
		              round(timeSNP/iter*(op$run_para$niter-iter),2),"seconds \n",sep=" ")
		  } 
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
      hyperparameters=c()
      counth=0
       hyperparameters = cbind(hyperparameters,varesave);counth=counth+1
			 hyperparameters = cbind(hyperparameters,scalesave) ;counth=counth+1
			 hyperparameters = cbind(hyperparameters,sigma_gsave) ;counth=counth+1
			if (op$update_para$df) {hyperparameters = cbind(hyperparameters,defsave);counth=counth+1}
			if (op$update_para$pi) {hyperparameters = cbind(hyperparameters,pisave) ;counth=counth+1}
			if(counth>0) rownames(hyperparameters) = seq(op$run_para$skip,op$run_para$niter,op$run_para$skip)
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
		  if(op$model %in% c("ssrrBLUP","ssBayesA")){
		    tmp=count/(itersave/run_para$skip)
		    bpvalue=2*apply(cbind(tmp,1-tmp),1,min)
		    if(GWA=="Win") Wprob=countqw/(itersave/run_para$skip) else Wprob=NULL
		  }
		  if(op$model %in% c("ssSSVS","ssBayesB")){
		    #if(op$model=="SSVS") postprob_est   = postprob_save/itersave else postprob_est=postprob/itersave
		    postprob_est   = postprob_save/itersave
		    if(GWA=="Win") Wpostprob_est =Wpostprob/n_saved_sample else Wpostprob_est=NULL
		    names(postprob_est)=colnames(Z)
		  } 
		  meanmu=Mcurrent[1:dimX]
			names(meanmu)=colnames(X)
			meang= Mcurrent[(dimX+1):dimW]
			names(meang)=colnames(M2)
			
			meane=Mcurrent2
			names(meane)=rownames(M1hat)
			
			varg = Qcurrent[(dimX+1):dimW]/itersave
			names(varg)=colnames(M2)
			
 	    save(iter,alphadef_save,alphascalea_save,file=paste(op$save.at, "alphadef",op$seed,".RData",sep=""))
 			
 			if(op$model %in% c("ssBayesB","ssSSVS")) save(iter,meanmu,meang,meane,postprob_est,varg,SNPvar,file = paste(op$save.at,"EffectsResults",op$seed,".RData",sep=""))
 				else save(iter,meanmu,meang,meane,varg,SNPvar,bpvalue,Wprob,file = paste(op$save.at,"EffectsResults",op$seed,".RData",sep=""))
 			
		} # "END" OF PROCESSING MCMC SAMPLING

  } # "END" OF MCMC SAMPLING



	names_hypers=c()
	 names_hypers=c(names_hypers,"vare")
	 names_hypers = c(names_hypers,"scale")
	 names_hypers = c(names_hypers,"varGenetic")
	if (op$update_para$df) {names_hypers = c(names_hypers,"df") }
	if (op$update_para$pi) {names_hypers = c(names_hypers,"pi") }

	sample_idx=c((op$run_para$burnIn/op$run_para$skip+1):(op$run_para$niter/op$run_para$skip))
	if(length(names_hypers)>1) hyper_est=apply(hyperparameters[sample_idx,],2,mean)
	else if (length(names_hypers)==1) hyper_est=mean(hyperparameters[sample_idx,])
	else hyper_est=c(vare,scale)
	
	colnames(hyperparameters)=names_hypers
	if(length(names_hypers)>1) names(hyper_est)=names_hypers  


	if(op$model %in% c("ssBayesB","ssSSVS")) BAout<-list(bhat=meanmu,ahat=meang, meane ,yhat=X0%*%meanmu+M0%*%meang+U0%*%meane,
	                                                 prob=postprob_est,Wprob=Wpostprob_est,eff_sample=effectiveSize(hyperparameters),
	                                                 hypers=hyperparameters,hyper_est=hyper_est,train=vtrain,y=y0,map=map,win=win,Wchr=chr,GWA=GWA,model=op$model)
    else BAout<-list(bhat=meanmu,ahat=meang, ehat=meane,yhat=X0%*%meanmu+M0%*%meang+U0%*%meane,
                     bpvalue=bpvalue,Wprob=Wprob,eff_sample=effectiveSize(hyperparameters),hypers=hyperparameters,
                     hyper_est=hyper_est,train=vtrain,y=y0,map=map,win=win,Wchr=chr,GWA=GWA,model=op$model)
	
	class(BAout)="ba"
   	return(BAout)
}  # END OF FUNCTION

#######################################################################################
#######################################################################################
