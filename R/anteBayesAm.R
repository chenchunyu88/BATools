###################################################################################################
###############             anteBayesA                      ########################################
###################################################################################################
anteBayesAm =function(op=NULL,y=NULL,Z=NULL,X=NULL,vtrain=NULL,GWA=NULL,map=NULL,ante_p=NULL)  
#  startpi is defined by the proportion of markers that ARE associate with genes
{   # START OF FUNCTION

  set.seed(op$seed)  
  pi_math = 3.14159265359	
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
  
  # Vectors for saved samples of hyperparameters
  if (op$update_para$df)    {defsave   =  array(0,op$run_para$niter/op$run_para$skip) }
  if (op$update_para$scale)  {scalesave =  array(0,op$run_para$niter/op$run_para$skip) }
  if (op$update_para$pi)     {pisave    =  array(0,op$run_para$niter/op$run_para$skip) } 
  varesave = array(0,op$run_para$niter/op$run_para$skip)
  
  if (op$update_para$mut) {mutsave   = array(0,op$run_para$niter/op$run_para$skip) }
  if (op$update_para$vart) { vartsave  = array(0,op$run_para$niter/op$run_para$skip)  }

  # To monitor MH acceptance ratios on degrees of freedom parameter
  alphadef_save = array(0,(op$run_para$niter-op$run_para$burnIn))
  # TO MONITOR MH ACCEPTANCE RATIOS ON SCALE PARAMETER (WHERE NECESSARY) POST BURN-IN
  alphascalea_save = array(0,(op$run_para$niter-op$run_para$burnIn)/op$run_para$skip)
  
  nx=rownames(X)
  ng=rownames(Z)
  np=names(y)
  idx <- Reduce(intersect, list(nx,ng,np))
  if(dim(X)[2]==1) X=as.matrix(X[idx,],ncol=1) else X=X[idx,]
  y=y[idx]
  Z=Z[idx,]
  
  dimX=dim(X)[2]
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
  mut = op$init$mut
  vart =op$init$vart
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
  assoct     = array(0,nSNP)
  SNPg      = array(0,nSNP) 
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
  Tlast= array(0,nSNP)
  Tcurrent= array(0,nSNP)
  TQlast= array(0,nSNP)
  TQcurrent= array(0,nSNP)

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
     # sample variance & effect for each SNP	
     ####################################################################################
     if (Pi_SNP == 1)  
     { # "START" OF anteBAYESA
       #####################################################################################
       # sample conditional variance for each SNP (anteBayes A)	
       #####################################################################################
       varq = (scalea*def + SNPeff*SNPeff)/rchisq(nSNP,def+1)  # joint sample

       #####################################################################################
       # sample conditional effect for each SNP (anteBayes A)	
       #####################################################################################
        startSNP = proc.time()
		    BayesC<- .Call("anteBayesAC",assoct,nSNP,dimX,nanim,Wstacked,WWdiag, theta, ycorr,varq,vare,SNPeff,SNPg,mut,vart)
		    theta=BayesC[[1]]
		    ycorr=BayesC[[2]]
		    varq=BayesC[[3]]
		    vare=BayesC[[4]]
		    SNPeff=BayesC[[5]]
		    SNPg=BayesC[[6]]
		    assoct=BayesC[[7]]
		    if(!is.null(ante_p)) assoct[ante_p]=0
		    G1=nSNP
      } # "END" FOR anteBAYESA


     

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

      
     #####################################################################################
     #  sample scale#
     #  based on Gelman Uniform(0,A) prior on scale parameter for "large" A.
     #####################################################################################
     if (op$update_para$scale && !(op$update_para$df))
      {  #"START" OF SAMPLE THE SCALE PARAMETER 
           shape_g = 0.5*(nSNP*def+1)+op$priors$shape_scale
           scale_g = 0.5*def*sum(1/varq)+op$priors$rate_scale
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
      #####################################################################################
      #  sample mean and variance for association parameters
      #####################################################################################
      if (op$update_para$mut)
      {
          ######sample mut#######
          sigma2hat_t=1/(1/op$priors$sigma2_m_t+length(which(assoct!=0))/vart) 
          muhat_t=((mean(assoct[which(assoct!=0)])*length(which(assoct!=0))/vart)+(op$priors$mu_m_t/op$priors$sigma2_m_t))*sigma2hat_t  
          mut = rnorm(1,mean=muhat_t,sd=sqrt(sigma2hat_t))   
      }
      ########sample vart#######################
      if(op$update_para$vart)
      {
	      vart = ( t(assoct[which(assoct!=0)]-mut)%*%(assoct[which(assoct!=0)]-mut))/rchisq(1,(length(assoct[which(assoct!=0)])-1))
      }

      if (itersave>1) 
      {     # "START" OF COMPUTE POSTERIOR MEANS & VARIANCES
           Mcurrent = Mlast + (theta-Mlast)/itersave
           Qcurrent = Qlast + (itersave-1)*((theta-Mlast)^2)/itersave
           Mlast = Mcurrent
           Qlast = Qcurrent 

           SNPvar=SNPvarlast + (varq-SNPvarlast)/itersave
           SNPvarlast=SNPvar

           Tcurrent=Tlast+ (assoct-Tlast)/itersave
           Tlast=Tcurrent

	         TQcurrent= TQlast+(itersave-1)*((assoct-Tlast)^2)/itersave
	         TQlast=TQcurrent
	         
	         ###Do some GWA
	         if (iter%%op$run_para$skip == 0) {
	           count=count+(SNPeff>0)
	           if(GWA=="Win"){
	             g_total<-Z%*%SNPeff
	             var_total=sum(g_total*g_total)/nanim-(sum(g_total)/nanim)^2
	             varw=calc.varw(win,Z,SNPeff)
	             qw=varw/var_total
	             countqw=countqw+(qw>0.01)
	           }
	         }
	         
     }  # "END" OF COMPUTE POSTERIOR MEANS & VARIANCES
     if (itersave==1) {
		    Mlast = theta;SNPvarlast=varq;Tlast= assoct;TQlast=rep(vart,length(assoct))
		    count=rep(0,nSNP)
		    if(GWA=="Win") countqw=rep(0,length(win))
	   }
       ###########################################################################################
       #store MCMC samples at every "skip"
       ###########################################################################################
       if (iter%%op$run_para$skip  == 0) 
        {
            varesave[iter/op$run_para$skip ] = vare
            if(op$update_para$df) {defsave[iter/op$run_para$skip ] = def}
            if(op$update_para$scale) {scalesave[iter/op$run_para$skip ] = scalea}
            if(op$update_para$pi) {pisave[iter/op$run_para$skip ] = Pi_SNP}
            if(op$update_para$mut){mutsave[iter/op$run_para$skip]=mut} 
	          if(op$update_para$vart){vartsave[iter/op$run_para$skip]=vart}
        }

       endSNP = proc.time()
       timeSNP = timeSNP + (endSNP[1]-startSNP[1]) 

       endtime = proc.time()
       timetotal = endtime[1]-starttime[1]
      if((op$print_mcmc$piter!=0) & iter%%op$print_mcmc$piter==0)
	    {
		    tmp=paste("iter= ",iter," vare= ",round(vare,6), "scale= ", round(scalea,6),
		              "\ntimepercycle= ", round(timeSNP/iter,3), "estimated time left=", round(timeSNP/iter*(op$run_para$niter-iter),2),"\n",sep=" ")
		
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
          #if (op$update_para$scale) {hyperparameters = cbind(hyperparameters,scalesave);counth=counth+1 }
          if (op$update_para$df) {hyperparameters = cbind(hyperparameters,defsave);counth=counth+1 }
          if (op$update_para$pi) {hyperparameters = cbind(hyperparameters,pisave);counth=counth+1 }
          if (op$update_para$mut) {hyperparameters = cbind(hyperparameters,mutsave);counth=counth+1}
 	        if (op$update_para$vart) {hyperparameters = cbind(hyperparameters,vartsave) ;counth=counth+1}
          if(counth)rownames(hyperparameters) = seq(op$run_para$skip,op$run_para$niter,op$run_para$skip)
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
          tmp=count/(itersave/run_para$skip)
          bpvalue=2*apply(cbind(tmp,1-tmp),1,min)
          if(GWA=="Win") Wprob=countqw/(itersave/run_para$skip) else Wprob=NULL
         
         
          meanmu= Mcurrent[1:dimX]
          names(meanmu)=colnames(X)
          meang= Mcurrent[(dimX+1):dimW]
          varg = Qcurrent[(dimX+1):dimW]/itersave
          names(meang)=names(varg)=colnames(Z)
          
          meant= Tcurrent[1:(nSNP-1)]
          vart=TQcurrent[1:(nSNP-1)]/itersave
          names(meant)=names(vart)=names(meang)[-1]
          save(iter,meanmu,meang,varg,SNPvar,meant,vart,bpvalue,Wprob,file = paste(op$save.at,"EffectsResults",op$seed,".RData",sep=""))
       

          #  monitor MH acceptance rates over cycle for the degrees of freedom term.
          save(iter,alphadef_save,alphascalea_save,file=paste(op$save.at, "alphadef",op$seed,".RData",sep=""))

       } # "END" OF PROCESSING MCMC SAMPLING

    } # "END" OF MCMC SAMPLING
    names_hypers=c()
    names_hypers=c(names_hypers,"vare")
    names_hypers = c(names_hypers,"scale") 
    if (op$update_para$df) {names_hypers = c(names_hypers,"df") }
    if (op$update_para$pi) {names_hypers = c(names_hypers,"pi") }
    if (op$update_para$mut) names_hypers = c(names_hypers,"mut")
    if (op$update_para$vart) names_hypers = c(names_hypers,"vart")
    sample_idx=c((op$run_para$burnIn/op$run_para$skip+1):(op$run_para$niter/op$run_para$skip))
    
    if(length(names_hypers)>1) hyper_est=apply(hyperparameters[sample_idx,],2,mean)
      else if (length(names_hypers)==1) hyper_est=mean(hyperparameters[sample_idx,])
      else hyper_est=c(vare,scale)
    
    colnames(hyperparameters)=names_hypers
    if(length(names_hypers)>1) names(hyper_est)=names_hypers  
    
    BAout<-list(bhat=meanmu,ahat=meang, yhat=X0%*%meanmu+Z0%*%meang,
                bpvalue=bpvalue,Wprob=Wprob,ante_t=meant,eff_sample=effectiveSize(hyperparameters),
                hypers=hyperparameters,idx=idx,hyper_est=hyper_est,train=vtrain,y=y0,map=map,win=win,
                Wchr=chr,GWA=GWA,model=op$model)
    class(BAout)="ba"
    return(BAout)
}  # END OF FUNCTION