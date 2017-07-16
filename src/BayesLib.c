#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <math.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>
#include <stdio.h>
#define N 2

SEXP anteBayesAC(SEXP cassoct, SEXP cnSNP, SEXP cdimX, SEXP cnanim, SEXP cW, SEXP cWWdiag, SEXP ctheta, SEXP cycorr, SEXP cvarq, SEXP cvarE, SEXP cSNPeff, SEXP cSNPg,SEXP cmut, SEXP cvart)
{
	double *xj, *W, *WWdiag, *theta, *ycorr, *varq, *SNPeff,*assoct,*SNPg;
	double rhs=0,varE,invlhs=0,mean=0,lhs=0,sd,mut,vart;
	int SNP,i, nSNP, nanim,dimX;
    SEXP list;
	
  	GetRNGstate();
	
	nSNP=INTEGER_VALUE(cnSNP);
	nanim=INTEGER_VALUE(cnanim);
	dimX=INTEGER_VALUE(cdimX);

	varE=NUMERIC_VALUE(cvarE); 
	mut=NUMERIC_VALUE(cmut);
	vart=NUMERIC_VALUE(cvart);

	PROTECT(cassoct=AS_NUMERIC(cassoct));
	assoct=NUMERIC_POINTER(cassoct); 
		
  	PROTECT(cW=AS_NUMERIC(cW));
	W=NUMERIC_POINTER(cW); 
        
    PROTECT(cWWdiag=AS_NUMERIC(cWWdiag));
    WWdiag=NUMERIC_POINTER(cWWdiag); 

    PROTECT(ctheta=duplicate(AS_NUMERIC(ctheta)));
    theta=NUMERIC_POINTER(ctheta);  
	
    PROTECT(cycorr=duplicate(AS_NUMERIC(cycorr)));
    ycorr=NUMERIC_POINTER(cycorr);

    PROTECT(cvarq=AS_NUMERIC(cvarq));
    varq=NUMERIC_POINTER(cvarq);
	
	PROTECT(cSNPeff=AS_NUMERIC(cSNPeff));
	SNPeff=NUMERIC_POINTER(cSNPeff);


	PROTECT(cSNPg=AS_NUMERIC(cSNPg));
	SNPg=NUMERIC_POINTER(cSNPg);

	

    xj=(double *) R_alloc(nanim,sizeof(double));

		
	for(SNP=0; SNP<nSNP; SNP++)
	{		
		rhs=0;	
		for(i=0; i<nanim; i++)
		{	
			xj[i]=W[i+(SNP+dimX)*nanim];
			ycorr[i]=ycorr[i]+SNPg[SNP]*xj[i];	
			rhs=rhs+xj[i]*ycorr[i]/varE;
		}
		
		if(SNP==(nSNP-1))
		{
			lhs=WWdiag[SNP+dimX]/varE+1.0/varq[SNP];
		}
		if(SNP<(nSNP-1))
		{
			lhs=WWdiag[SNP+dimX]/varE+1.0/varq[SNP]+R_pow_di(assoct[SNP],2)/varq[SNP+1];
		}

		if(SNP==0)
		{
			rhs=rhs+assoct[SNP]*SNPg[SNP+1]/varq[SNP+1];
		}
		if(SNP==(nSNP-1))
		{
			rhs=rhs+assoct[SNP-1]*SNPg[SNP-1]/varq[SNP];
		}
		if(SNP>0 && SNP<(nSNP-1))
		{
			rhs=rhs+assoct[SNP]*SNPg[SNP+1]/varq[SNP+1]+assoct[SNP-1]*SNPg[SNP-1]/varq[SNP];
		}
		
		invlhs=1.0/lhs;
		mean=invlhs*rhs;
		sd=sqrt(invlhs);
		SNPg[SNP]= rnorm(mean,sd);		
		for(i=0; i<nanim; i++)
		{
			ycorr[i]=ycorr[i]-SNPg[SNP]*xj[i];
		}
		theta[SNP+dimX]=SNPg[SNP];
	}  
	//Rprintf("rhs %f \n",rhs);
	//Rprintf("invlhs %f \n",invlhs);

	for(SNP=0; SNP<nSNP; SNP++)
	{
		if (SNP==0) { SNPeff[SNP]=SNPg[SNP];}
        if (SNP>0) { SNPeff[SNP]=SNPg[SNP]-assoct[SNP-1]*SNPg[SNP-1];}
    }

    for(SNP=0; SNP<(nSNP-1); SNP++)
	{
        rhs=SNPg[SNP]*SNPg[SNP+1]/varq[SNP+1]+mut/vart;  
        lhs=R_pow_di(SNPg[SNP],2)/varq[SNP+1]+1/vart;
        invlhs=1/lhs;
        mean=invlhs*rhs;
		sd=sqrt(invlhs);
        assoct[SNP]=rnorm(mean,sd);
    }
		
    // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 7));
	// attaching theta vector to list:
	SET_VECTOR_ELT(list, 0, ctheta);
	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 1, cycorr);
	
	SET_VECTOR_ELT(list, 2, cvarq);
	
	SET_VECTOR_ELT(list, 3, cvarE);
	// attaching SNPeff vector to list:
	SET_VECTOR_ELT(list, 4, cSNPeff);

	SET_VECTOR_ELT(list, 5, cSNPg);

	SET_VECTOR_ELT(list, 6, cassoct);


  	PutRNGstate();

  	UNPROTECT(9);
  	return(list);
	
}



SEXP anteBayesBC(SEXP cassoct, SEXP cnSNP, SEXP cdimX, SEXP cnanim, SEXP cW, SEXP cM, SEXP cWWdiag, SEXP ctheta, SEXP cycorr, SEXP cvarq, SEXP cvarE, SEXP cSNPeff, SEXP cSNPg, SEXP cpostprob,SEXP cscalea,SEXP cdef, SEXP calphametrovar, SEXP cpi,SEXP cG1, SEXP cmut, SEXP cvart)
{
	double *xj, *W,*M, *theta, *ycorr, *varq, *SNPeff,*assoct,*SNPg,*alphametrovar;
	double rhs=0,varE,invlhs=0,mean=0,lhs=0,sd,mut,vart,def,scalea,pi,MpM;
	double u,v_no, v_yes,logDataNullModel,logDataOld,alphacheck=0,varCandidate=0,logDataNew,acceptProb;
	int SNP,i, nSNP ,*postprob, nanim,*G1,dimX, mhiter;
    SEXP list,alphacheck1;
	
  	GetRNGstate();
	
	nSNP=INTEGER_VALUE(cnSNP);
	nanim=INTEGER_VALUE(cnanim);
	dimX=INTEGER_VALUE(cdimX);

	varE=NUMERIC_VALUE(cvarE); 
	mut=NUMERIC_VALUE(cmut);
	vart=NUMERIC_VALUE(cvart);
	
	def=NUMERIC_VALUE(cdef);
	scalea=NUMERIC_VALUE(cscalea);
	
	pi=NUMERIC_VALUE(cpi);
	
	PROTECT(cG1=AS_INTEGER(cG1));
	G1=INTEGER_POINTER(cG1);

	PROTECT(cassoct=AS_NUMERIC(cassoct));
	assoct=NUMERIC_POINTER(cassoct); 
		
  	PROTECT(cW=AS_NUMERIC(cW));
	W=NUMERIC_POINTER(cW); 
	
	PROTECT(cM=AS_NUMERIC(cM));
	M=NUMERIC_POINTER(cM);
        
    //PROTECT(cWWdiag=AS_NUMERIC(cWWdiag));
    //WWdiag=NUMERIC_POINTER(cWWdiag); 

    PROTECT(ctheta=duplicate(AS_NUMERIC(ctheta)));
    theta=NUMERIC_POINTER(ctheta);  
	
    PROTECT(cycorr=duplicate(AS_NUMERIC(cycorr)));
    ycorr=NUMERIC_POINTER(cycorr);

    PROTECT(cvarq=AS_NUMERIC(cvarq));
    varq=NUMERIC_POINTER(cvarq);
	
	PROTECT(cSNPeff=AS_NUMERIC(cSNPeff));
	SNPeff=NUMERIC_POINTER(cSNPeff);


	PROTECT(cSNPg=AS_NUMERIC(cSNPg));
	SNPg=NUMERIC_POINTER(cSNPg);
	
	PROTECT(cpostprob=AS_INTEGER(cpostprob));
	postprob=INTEGER_POINTER(cpostprob);

	PROTECT(calphametrovar=AS_NUMERIC(calphametrovar));
	alphametrovar=NUMERIC_POINTER(calphametrovar);
	
	PROTECT(alphacheck1 = allocVector(REALSXP,1));
	

    xj=(double *) R_alloc(nanim,sizeof(double));

		
	for(SNP=(nSNP-1); SNP>=0; SNP--)
	{		
	
		for(i=0; i<nanim; i++)
		{	
			xj[i]=M[i+SNP*nanim];
			ycorr[i]+=SNPeff[SNP]*xj[i];				
		}		
		if(SNP==(nSNP-1))
		{
			for(i=0; i<nanim; i++)
			{
				M[i+SNP*nanim]=W[i+(SNP+dimX)*nanim];
			}
		}
	
		if(SNP<(nSNP-1))
		{
			for(i=0; i<nanim; i++)
			{
				M[i+SNP*nanim]=M[i+(SNP+1)*nanim]*assoct[SNP]+W[i+(SNP+dimX)*nanim];
			}
		}
		rhs=0;
			
		for(i=0; i<nanim; i++)
		{	
			rhs+=M[i+SNP*nanim]*ycorr[i];			
		}		
		
		MpM=0;
		for(i=0; i<nanim; i++)
		{	
			MpM+=M[i+SNP*nanim]*M[i+SNP*nanim];			
		}		

		v_no=MpM*varE;
		v_yes=MpM*MpM*varq[SNP]+v_no;
		logDataNullModel=-0.5*(log(v_no)+R_pow_di(rhs,2)/v_no);
		if(varq[SNP]>0.0)
		{ 
			logDataOld   = -0.5*(log(v_yes) + R_pow_di(rhs,2)/v_yes); 
		}
		else
		{
			logDataOld = logDataNullModel;
		}
		alphacheck=0;
		for (mhiter=0;mhiter<10;mhiter++)
		{
			u = runif(0.0,1.0);
			varCandidate=0;
			if(u<pi)
			{
				varCandidate=scalea*def/rchisq(def);
			}
			if(varCandidate>0.0)
			{
				v_yes=MpM*MpM*varCandidate+v_no;
				logDataNew= -0.5*(log(v_yes) + rhs*rhs/v_yes);
			}
			else
			{
				logDataNew=logDataNullModel;
			}
			
			acceptProb=exp(logDataNew-logDataOld);
			u = runif(0.0,1.0);
			if(u<acceptProb)
			{
				varq[SNP]=varCandidate;
				logDataOld=logDataNew;
			}
			alphacheck+=acceptProb;		
		}
		alphametrovar[SNP]+=alphacheck/10;
		
		if(varq[SNP]>0)
		{
			G1[0]++;
			postprob[SNP]++;
			lhs=MpM/varE+1.0/varq[SNP];
			invlhs=1/lhs;
			mean=invlhs*rhs/varE;
			sd=sqrt(invlhs);
			SNPeff[SNP]=rnorm(mean,sd);
			for(i=0; i<nanim; i++)
			{	
				ycorr[i]-=	M[i+SNP*nanim]*SNPeff[SNP];	
			}			
		}
		else if(varq[SNP]==0)
		{
			SNPeff[SNP]=0.0;
		}		
	}	
	
	for (SNP=0;SNP<nSNP;SNP++)
	{
		if (SNP==0) { SNPg[SNP]=SNPeff[SNP];}
		if (SNP > 0) { SNPg[SNP]=assoct[SNP-1]*SNPg[SNP-1]+SNPeff[SNP];}
		theta[dimX+SNP] = SNPg[SNP];
    }
	
    for(SNP=0; SNP<(nSNP-1); SNP++)
	{
        if(varq[SNP+1]!=0)
		{
			rhs=SNPg[SNP]*SNPg[SNP+1]/varq[SNP+1]+mut/vart;  
			lhs=(SNPg[SNP]*SNPg[SNP])/varq[SNP+1]+1/vart;
			invlhs=1/lhs;
			mean=invlhs*rhs;
			sd=sqrt(invlhs);
			assoct[SNP]=rnorm(mean,sd);
		}
		if (varq[SNP+1]==0 && SNPg[SNP]!=0)
		{
            assoct[SNP]=SNPg[SNP+1]/SNPg[SNP];
		}
		if (varq[SNP+1]==0 && SNPg[SNP]==0)
		{
			assoct[SNP]=rnorm(mut,sqrt(vart));
		}		
    }
		
    // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 14));
	
	SET_VECTOR_ELT(list, 0, ctheta);
	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 1, cycorr);
	
	SET_VECTOR_ELT(list, 2, cvarq);
	
	SET_VECTOR_ELT(list, 3, cvarE);
	// attaching SNPeff vector to list:
	SET_VECTOR_ELT(list, 4, cSNPeff);
	
	SET_VECTOR_ELT(list, 5, cpostprob);
	
	SET_VECTOR_ELT(list, 6, calphametrovar);

	SET_VECTOR_ELT(list, 7, cG1);
	
	SET_VECTOR_ELT(list, 8, cdef);
	
	SET_VECTOR_ELT(list, 9, cscalea);
	
	SET_VECTOR_ELT(list, 10, cpi);
	
	REAL(alphacheck1)[0]=alphacheck;
	
	
	SET_VECTOR_ELT(list, 11, alphacheck1);	

	SET_VECTOR_ELT(list, 12, cSNPg);

	SET_VECTOR_ELT(list, 13, cassoct);


  	PutRNGstate();

  	UNPROTECT(13);
  	return(list);
	
}


SEXP BayesAC(SEXP cnSNP, SEXP cdimX, SEXP cnanim, SEXP cW, SEXP cWWdiag, SEXP ctheta, SEXP cycorr, SEXP cvarq, SEXP cvarE)
{
	double *xj, *W, *WWdiag, *theta, *ycorr, *varq;
	double rhs,varE,lhs;
	int j,i, nSNP, nanim,dimX;
    SEXP list;
	
  	GetRNGstate();
	
	  nSNP=INTEGER_VALUE(cnSNP);
	  nanim=INTEGER_VALUE(cnanim);
	  dimX=INTEGER_VALUE(cdimX);

	  varE=NUMERIC_VALUE(cvarE); 
		
  	PROTECT(cW=AS_NUMERIC(cW));
	  W=NUMERIC_POINTER(cW); 
        
    PROTECT(cWWdiag=AS_NUMERIC(cWWdiag));
    WWdiag=NUMERIC_POINTER(cWWdiag); 

    PROTECT(ctheta=duplicate(AS_NUMERIC(ctheta)));
    theta=NUMERIC_POINTER(ctheta);  
	
    PROTECT(cycorr=duplicate(AS_NUMERIC(cycorr)));
    ycorr=NUMERIC_POINTER(cycorr);

    PROTECT(cvarq=AS_NUMERIC(cvarq));
    varq=NUMERIC_POINTER(cvarq);	

    xj=(double *) R_alloc(nanim,sizeof(double));

			
	  for(j=0; j<nSNP; j++)
	  {		
		  rhs=0;
		for(i=0; i<nanim; i++)
		{	
			xj[i]=W[i+(j+dimX)*nanim];
			ycorr[i]=ycorr[i]+theta[j+dimX]*xj[i];
			rhs+=xj[i]*ycorr[i];
		}
		
		rhs=rhs/varE;
		lhs=WWdiag[j+dimX]/varE+1.0/varq[j];
		theta[j+dimX]=rhs/lhs + sqrt(1.0/lhs)*rnorm(0,1) ;

		for(i=0; i<nanim; i++)
		{
			ycorr[i]=ycorr[i]-theta[j+dimX]*xj[i];
		}
		
	}  
		
    // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 2));
	// attaching theta vector to list:
	SET_VECTOR_ELT(list, 0, ctheta);
	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 1, cycorr);

  	PutRNGstate();

  	UNPROTECT(6);
  	return(list);
	
}


SEXP BayesBC(SEXP cnSNP, SEXP cdimX, SEXP cnanim, SEXP cW, SEXP cWWdiag, SEXP ctheta, SEXP cycorr, SEXP cvarq, SEXP cvarE, SEXP cpostprob,SEXP cscalea,SEXP cdef, SEXP calphametrovar,SEXP citer, SEXP cpi, SEXP cG1)
{
	int  *postprob,inc=1;
	double *W, *WWdiag,*xj,*theta, *ycorr, *varq,*alphametrovar;
	double rhs,varE,invlhs,mean,lhs,sd,v_no,v_yes,logDataNullModel;
	double logDataOld,alphacheck=0,def,scalea,acceptProb=0;
	double logDataNew,varCandidate,u,pi,b;
	int SNP,nSNP,nanim,*G1,dimX,mhiter;
    	SEXP list,alphacheck1,alpha;
	
  	GetRNGstate();
	
	nSNP=INTEGER_VALUE(cnSNP);
	nanim=INTEGER_VALUE(cnanim);
	dimX=INTEGER_VALUE(cdimX);
	def=NUMERIC_VALUE(cdef);
	//iter=INTEGER_VALUE(citer);
	//G1=INTEGER_VALUE(cG1);
	
	scalea=NUMERIC_VALUE(cscalea);
	varE=NUMERIC_VALUE(cvarE); 
	pi=NUMERIC_VALUE(cpi);
	

		//G1=INTEGER_VALUE(cG1);
		
	PROTECT(cG1=AS_INTEGER(cG1));
	G1=INTEGER_POINTER(cG1); 	
	
  	PROTECT(cW=AS_NUMERIC(cW));
	W=NUMERIC_POINTER(cW); 
        
    PROTECT(cWWdiag=AS_NUMERIC(cWWdiag));
    WWdiag=NUMERIC_POINTER(cWWdiag); 

    PROTECT(ctheta=duplicate(AS_NUMERIC(ctheta)));
    theta=NUMERIC_POINTER(ctheta);  
	
    PROTECT(cycorr=duplicate(AS_NUMERIC(cycorr)));
    ycorr=NUMERIC_POINTER(cycorr);

    PROTECT(cvarq=AS_NUMERIC(cvarq));
    varq=NUMERIC_POINTER(cvarq);
	
	PROTECT(cpostprob=AS_INTEGER(cpostprob));
	postprob=INTEGER_POINTER(cpostprob);
	
	PROTECT(calphametrovar=AS_NUMERIC(calphametrovar));
	alphametrovar=NUMERIC_POINTER(calphametrovar);
	
	PROTECT(alpha = allocVector(REALSXP,1));
	
	PROTECT(alphacheck1 = allocVector(REALSXP,1));

    //xj=(int *) R_alloc(nanim,sizeof(int));

		
	for(SNP=0; SNP<nSNP; SNP++)
	{
		b=theta[SNP+dimX];
		xj=W+(SNP+dimX)*nanim;
		F77_NAME(daxpy)(&nanim,&b,xj,&inc,ycorr,&inc);
		rhs=F77_NAME(ddot)(&nanim,xj,&inc,ycorr,&inc);
		/*rhs=0;
		for(i=0; i<nanim; i++)
		{	
			xj[i]=W[i+(SNP+dimX)*nanim];
			ycorr[i]=ycorr[i]+theta[SNP+dimX]*xj[i];
			rhs=rhs+xj[i]*ycorr[i];
		}
		*/
		v_no=WWdiag[SNP+dimX]*varE;
		v_yes=R_pow_di(WWdiag[SNP+dimX],2)*varq[SNP]+v_no; //R_pow_di(x,2)=x^2

		logDataNullModel =   -0.5*(log(v_no)+R_pow_di(rhs,2)/v_no);
		
		if (varq[SNP]>0.0)
		{
			logDataOld =  -0.5*(log(v_yes)+R_pow_di(rhs,2)/v_yes);
		}
		else
		{
			logDataOld=logDataNullModel;
		}
		
		alphacheck=0;
		for (mhiter=0; mhiter<10;mhiter++) 
		{  
			u = runif(0.0,1.0);
			varCandidate=0;
			if(u<pi)
			{
				varCandidate=scalea*def/(rchisq(def));
			}
			if(varCandidate>0.0)
			{
				v_yes =R_pow_di(WWdiag[dimX+SNP],2)*varCandidate + v_no;
				logDataNew = -0.5*(log(v_yes)+R_pow_di(rhs,2)/v_yes);
			}
			else
			{
				logDataNew = logDataNullModel;					
			}

			acceptProb = exp(logDataNew-logDataOld);
			u=runif(0,1);
			if(u<acceptProb)
			{
				varq[SNP]=varCandidate;
                logDataOld=logDataNew;
			}	
			alphacheck=alphacheck+acceptProb;
			//if(varq[SNP]<1E-8) varq[SNP]=0.0;
				
		}
		alphametrovar[SNP] = alphametrovar[SNP] + alphacheck/10;
		if(varq[SNP]>0.0)
		{
			G1[0]=G1[0]+1;
			postprob[SNP]=postprob[SNP]+1;
	
			lhs=WWdiag[SNP+dimX]/varE+1.0/varq[SNP];
			invlhs=1.0/lhs;
			mean=invlhs*rhs/varE;
			sd=sqrt(invlhs);
			theta[SNP+dimX]=rnorm(mean,sd);
			b=-theta[SNP+dimX];
			F77_NAME(daxpy)(&nanim,&b,xj,&inc,ycorr,&inc);
			/*
			for(i=0; i<nanim; i++)
			{		
				ycorr[i]=ycorr[i]-theta[SNP+dimX]*xj[i];
			} 
			*/
		}
		else if (varq[SNP]==0.0)
		{
			theta[dimX+SNP]=0.0;
		}
	}

	
				
		
	//Rprintf("G1 in C %d",G1[0]); 
		
    // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 7));
	// attaching theta vector to list:
	SET_VECTOR_ELT(list, 0, ctheta);
	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 1, cycorr);

	SET_VECTOR_ELT(list, 2, cpostprob);
	
	SET_VECTOR_ELT(list, 3, calphametrovar);
	
	SET_VECTOR_ELT(list, 4, cG1);
	
	REAL(alpha)[0]=acceptProb;
	
	REAL(alphacheck1)[0]=alphacheck;
	
	SET_VECTOR_ELT(list, 5, alpha);
	
	SET_VECTOR_ELT(list, 6, alphacheck1);
	
	
  	PutRNGstate();

  	UNPROTECT(11);
  	return(list);
}


SEXP BayesACL(SEXP cnSNP, SEXP cdimX, SEXP cnanim, SEXP cW, SEXP cWWdiag, SEXP ctheta, SEXP cycorr, SEXP cvarq, SEXP cvarE)
{
	double *xj, *W, *WWdiag, *theta, *ycorr, *varq;
	double rhs,varE,lhs,b;
	int j, nSNP, nanim,dimX,inc=1;
    SEXP list;
	
  	GetRNGstate();
	
	nSNP=INTEGER_VALUE(cnSNP);
	nanim=INTEGER_VALUE(cnanim);
	dimX=INTEGER_VALUE(cdimX);

	varE=NUMERIC_VALUE(cvarE); 
		
  	PROTECT(cW=AS_NUMERIC(cW));
	W=NUMERIC_POINTER(cW); 
        
    PROTECT(cWWdiag=AS_NUMERIC(cWWdiag));
    WWdiag=NUMERIC_POINTER(cWWdiag); 

    PROTECT(ctheta=duplicate(AS_NUMERIC(ctheta)));
    theta=NUMERIC_POINTER(ctheta);  
	
    PROTECT(cycorr=duplicate(AS_NUMERIC(cycorr)));
    ycorr=NUMERIC_POINTER(cycorr);

    PROTECT(cvarq=AS_NUMERIC(cvarq));
    varq=NUMERIC_POINTER(cvarq);
	

    //xj=(double *) R_alloc(nanim,sizeof(double));

			
	for(j=0; j<nSNP; j++)
	{		
		b=theta[j+dimX];
		xj=W+(j+dimX)*nanim;
		F77_NAME(daxpy)(&nanim,&b,xj,&inc,ycorr,&inc);
		rhs=F77_NAME(ddot)(&nanim,xj,&inc,ycorr,&inc);
		//for(i=0; i<nanim; i++)
		//{	
		//	xj[i]=W[i+(j+dimX)*nanim];
		//	ycorr[i]=ycorr[i]+theta[j+dimX]*xj[i];
		//	rhs+=xj[i]*ycorr[i];
		//}
		rhs=rhs/varE;
		lhs=WWdiag[j+dimX]/varE+1.0/varq[j];
		theta[j+dimX]=rhs/lhs + sqrt(1.0/lhs)*rnorm(0,1) ;
		b=-theta[j+dimX];
		F77_NAME(daxpy)(&nanim,&b,xj,&inc,ycorr,&inc);
		//for(i=0; i<nanim; i++)
		//{
		//	ycorr[i]=ycorr[i]-theta[j+dimX]*xj[i];
		//}
		
	}  
		
    // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 2));
	// attaching theta vector to list:
	SET_VECTOR_ELT(list, 0, ctheta);
	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 1, cycorr);

  	PutRNGstate();

  	UNPROTECT(6);
  	return(list);
	
}



SEXP BayesCC(SEXP cnSNP, SEXP cdimX, SEXP cnanim, SEXP cW, SEXP cWWdiag, SEXP ctheta, SEXP cycorr, SEXP cvarq, SEXP cvarE,SEXP cpi,SEXP cphi,SEXP chratio,SEXP ca)
{
	double *xj, *W, *WWdiag, *theta, *ycorr, *varq,*phi,*hratio;
	double rhs,varE,lhs,pi;
	int j,i, nSNP, nanim,dimX;
	double h1,h0,varG,a; //for phi
    SEXP list;
	
  	GetRNGstate();
	
	nSNP=INTEGER_VALUE(cnSNP);
	nanim=INTEGER_VALUE(cnanim);
	dimX=INTEGER_VALUE(cdimX);

	varE=NUMERIC_VALUE(cvarE); 
	pi=NUMERIC_VALUE(cpi);
	a=NUMERIC_VALUE(ca);
		
  	PROTECT(cW=AS_NUMERIC(cW));
	  W=NUMERIC_POINTER(cW); 
        
    PROTECT(cWWdiag=AS_NUMERIC(cWWdiag));
    WWdiag=NUMERIC_POINTER(cWWdiag); 

    PROTECT(ctheta=duplicate(AS_NUMERIC(ctheta)));
    theta=NUMERIC_POINTER(ctheta);  
	
    PROTECT(cycorr=duplicate(AS_NUMERIC(cycorr)));
    ycorr=NUMERIC_POINTER(cycorr);

    PROTECT(cvarq=AS_NUMERIC(cvarq));
    varq=NUMERIC_POINTER(cvarq);
		
	PROTECT(cphi=AS_NUMERIC(cphi));
	phi=NUMERIC_POINTER(cphi);
		
	PROTECT(chratio=AS_NUMERIC(chratio));
	hratio=NUMERIC_POINTER(chratio);


    xj=(double *) R_alloc(nanim,sizeof(double));


	varG=varq[0];
		
	for(j=0; j<nSNP; j++)
	{		
		rhs=0;
		//h1=log1p(WWdiag[j+dimX]*varG/varE)-log1p(WWdiag[j+dimX]*varG/varE/100);
		//h2=1/(WWdiag[j+dimX]/varE+100/varG)-1/(WWdiag[j+dimX]/varE+1/varG);
		//Rprintf("h1=%f\n",h1);
		//Rprintf("h2=%f\n",h2);
		for(i=0; i<nanim; i++)
		{	
			xj[i]=W[i+(j+dimX)*nanim];
			ycorr[i]=ycorr[i]+theta[j+dimX]*xj[i];
			rhs+=xj[i]*ycorr[i];
			//Rprintf("rhs=%f\n",rhs);
		}
		h1=Rf_dnorm4(theta[j+dimX],0,sqrt(varG),0);
		h0=Rf_dnorm4(theta[j+dimX],0,sqrt(varG/a),0);
		//logh=0.5*(h1+h2*R_pow_di(rhs/varE,2));
		//Rprintf("rhs=%f\n",rhs);
		//Rprintf("vare=%f\n",varE);
		
		//Rprintf("h2=%f\n",h2);
		//hratio[j]=R_pow(M_E,0.5*logh);
		//Rprintf("hratio=%f\n",hratio);
		hratio[j]=h0/h1;
		phi[j]=rbinom(1,pi/(hratio[j]*(1-pi)+pi));
		//phi[j]=rbinom(1,h1*pi);
		//Rprintf("phi=%f\n",phi);
		
		rhs=rhs/varE;
		lhs=WWdiag[j+dimX]/varE+1.0/(varG*((1-phi[j])/a+phi[j]));
		theta[j+dimX]=rhs/lhs + sqrt(1.0/lhs)*rnorm(0,1) ;

		for(i=0; i<nanim; i++)
		{
			ycorr[i]=ycorr[i]-theta[j+dimX]*xj[i];
		}
		
	}  
		
    // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 4));
	// attaching theta vector to list:
	SET_VECTOR_ELT(list, 0, ctheta);
	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 1, cycorr);
	
	SET_VECTOR_ELT(list, 2, cphi);
	
	SET_VECTOR_ELT(list, 3, chratio);
  	PutRNGstate();

  	UNPROTECT(8);
  	return(list);
	
}

/*
SEXP BayesCCm(SEXP cnSNP, SEXP cdimX, SEXP cnanim, SEXP cW, SEXP cWWdiag, SEXP ctheta, SEXP cycorr, SEXP cvarq, SEXP cvarE,SEXP cpi,SEXP cphi,SEXP chratio,SEXP ca)
{
	//Marginalized version
	double *xj, *W, *WWdiag, *theta, *ycorr, *varq,*phi,*hratio;
	double rhs,varE,lhs,pi;
	int j,i, nSNP, nanim,dimX;
	double h1,h0,varG,a; //for phi
    SEXP list;
	
  	GetRNGstate();
	
	nSNP=INTEGER_VALUE(cnSNP);
	nanim=INTEGER_VALUE(cnanim);
	dimX=INTEGER_VALUE(cdimX);

	varE=NUMERIC_VALUE(cvarE); 
	pi=NUMERIC_VALUE(cpi);
	a=NUMERIC_VALUE(ca);
		
  	PROTECT(cW=AS_NUMERIC(cW));
	W=NUMERIC_POINTER(cW); 
        
    PROTECT(cWWdiag=AS_NUMERIC(cWWdiag));
    WWdiag=NUMERIC_POINTER(cWWdiag); 

    PROTECT(ctheta=duplicate(AS_NUMERIC(ctheta)));
    theta=NUMERIC_POINTER(ctheta);  
	
    PROTECT(cycorr=duplicate(AS_NUMERIC(cycorr)));
    ycorr=NUMERIC_POINTER(cycorr);

    PROTECT(cvarq=AS_NUMERIC(cvarq));
    varq=NUMERIC_POINTER(cvarq);
		
	PROTECT(cphi=AS_NUMERIC(cphi));
	phi=NUMERIC_POINTER(cphi);
		
	PROTECT(chratio=AS_NUMERIC(chratio));
	hratio=NUMERIC_POINTER(chratio);


    xj=(double *) R_alloc(nanim,sizeof(double));
	


	varG=varq[0];
		
	for(j=0; j<nSNP; j++)
	{		
		rhs=0;
		//h1=log1p(WWdiag[j+dimX]*varG/varE)-log1p(WWdiag[j+dimX]*varG/varE/100);
		//h2=1/(WWdiag[j+dimX]/varE+100/varG)-1/(WWdiag[j+dimX]/varE+1/varG);
		//Rprintf("h1=%f\n",h1);
		//Rprintf("h2=%f\n",h2);
		Zy[j]=0;
		for(i=0; i<nanim; i++)
		{	
			xj[i]=W[i+(j+dimX)*nanim];
			ycorr[i]=ycorr[i]+theta[j+dimX]*xj[i];
			rhs+=xj[i]*ycorr[i];
			//Rprintf("rhs=%f\n",rhs);
		}
		//h1=Rf_dnorm4(theta[j+dimX],0,sqrt(varG),0);
		//h0=Rf_dnorm4(theta[j+dimX],0,sqrt(varG/a),0);
		//logh=0.5*(h1+h2*R_pow_di(rhs/varE,2));
		//Rprintf("rhs=%f\n",rhs);
		//Rprintf("vare=%f\n",varE);
		
		//Rprintf("h2=%f\n",h2);
		//hratio[j]=R_pow(M_E,0.5*logh);
		//Rprintf("hratio=%f\n",hratio);
		logh=R_log(1+WWdiag[j+dimX]*varG/varE)-R_log(1+WWdiag[j+dimX]*varG/varE/a);
		logh+=(rhs/varE)^2*(1/(WWdiag[j+dimX]/varE+a/varG)-1/(WWdiag[j+dimX]/varE-1/varG));
		logh/=2;
		hratio[j]=exp(logh);
		phi[j]=rbinom(1,pi/(hratio[j]*(1-pi)+pi));
		//phi[j]=rbinom(1,h1*pi);
		//Rprintf("phi=%f\n",phi);
		
		rhs=rhs/varE;
		lhs=WWdiag[j+dimX]/varE+1.0/(varG*((1-phi[j])/a+phi[j]));
		theta[j+dimX]=rhs/lhs + sqrt(1.0/lhs)*rnorm(0,1) ;

		for(i=0; i<nanim; i++)
		{
			ycorr[i]=ycorr[i]-theta[j+dimX]*xj[i];
		}
		
	}  
		
    // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 4));
	// attaching theta vector to list:
	SET_VECTOR_ELT(list, 0, ctheta);
	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 1, cycorr);
	
	SET_VECTOR_ELT(list, 2, cphi);
	
	SET_VECTOR_ELT(list, 3, chratio);
  	PutRNGstate();

  	UNPROTECT(8);
  	return(list);
	
}

*/

SEXP sampleeff(SEXP cp, SEXP cn, SEXP cZ, SEXP cZ2, SEXP ceff, SEXP cycorr, SEXP cvarE, SEXP cAvaru)
{
	double *xj, *Z, *Z2, *eff, *ycorr,*Avaru,*aj;
	double rhs,varE,lhs;
	int j,i,p,n;
        SEXP list;
	
  	GetRNGstate();
	
	p=INTEGER_VALUE(cp);
	n=INTEGER_VALUE(cn);

	varE=NUMERIC_VALUE(cvarE); 
		
  	PROTECT(cZ=AS_NUMERIC(cZ));
	Z=NUMERIC_POINTER(cZ); 
        
	PROTECT(cZ2=AS_NUMERIC(cZ2));
	Z2=NUMERIC_POINTER(cZ2); 

	PROTECT(cycorr=AS_NUMERIC(cycorr));
	ycorr=NUMERIC_POINTER(cycorr);

	//PROTECT(cvareff=AS_NUMERIC(cvareff));
	//vareff=NUMERIC_POINTER(cvareff);
	
	PROTECT(ceff=AS_NUMERIC(ceff));
	eff=NUMERIC_POINTER(ceff);
	
	PROTECT(cAvaru=AS_NUMERIC(cAvaru));
	Avaru=NUMERIC_POINTER(cAvaru);
	

    	xj=(double *) R_alloc(n,sizeof(double));
	aj=(double *) R_alloc(p,sizeof(double));
	/*
	Rprintf("p=%d\n",p);
	Rprintf("nanim=%d\n",n); //Avaru[1:5,1:5]
	Rprintf("vare=%f\n",varE);
	for(j=0; j<p; j++)
	{
		Rprintf("eff[%d]=%f\n",j,eff[j]);
		Rprintf("ycorr[%d]=%f\n",j,ycorr[j]);
		Rprintf("Z2[%d]=%f\n",j,Z2[j]);	

	}
	Rprintf("Avaru[%d]=%f\n",31,Avaru[31]);	
	*/		
	for(j=0; j<p; j++)
	{
		//Rprintf("eff[j]=%f\n",eff[j]);		
		rhs=0;
		for(i=0; i<n; i++)
		{	
			xj[i]=Z[i+j*n];
		
			ycorr[i]=ycorr[i]+eff[j]*xj[i];
			//Rprintf("xj[%d]=%f\n",i,xj[i]);
			//Rprintf("aj[%d]=%f\n",i,aj[i]);
			//Rprintf("Avaru[%d]=%f\n",i,aj[i+j*n]);
			rhs+=xj[i]*ycorr[i]/varE;
		}
		for(i=0; i<p; i++)
		{
				aj[i]=Avaru[i+j*p];
				if(i!=j) rhs-=aj[i]*eff[i];
		}
		//Rprintf("rhs=%f\n",rhs);
		//rhs+=aj[j]*eff[j];
		//Rprintf("rhs=%f\n",rhs);
		lhs=Z2[j]/varE+aj[j];
		//Rprintf("aj[%d]=%f\n",j,aj[j]);
		//Rprintf("ajj[%d]=%f\n",j,aj[j+j*n]);
		//Rprintf("lhs=%f\n",lhs);
		eff[j]=rhs/lhs + sqrt(1.0/lhs)*rnorm(0,1) ;
		//Rprintf("eff[j]=%f\n",eff[j]);
		for(i=0; i<n; i++)
		{
			ycorr[i]=ycorr[i]-eff[j]*xj[i];
		}
		
	}  
		
    // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 2));
	// attaching theta vector to list:
	SET_VECTOR_ELT(list, 0, ceff);
	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 1, cycorr);

  	PutRNGstate();

  	UNPROTECT(6);
  	return(list);
	
}




SEXP GBLUP_GSRU(SEXP cnSNP, SEXP cdimX,SEXP cdimQ, SEXP cnanim, SEXP cW, SEXP cWWdiag,SEXP cWWdiagadd, SEXP ctheta, SEXP cycorr, SEXP cvarq, SEXP cvarE)
{
	double *xj, *newtheta, *W, *WWdiag, *WWdiagadd,*theta, *ycorr, *varq;
	double rhs,varE,lhs;
	int j,i, nSNP, nanim,dimX,dimQ;
    SEXP list;
	
  	GetRNGstate();
	
	nSNP=INTEGER_VALUE(cnSNP);
	nanim=INTEGER_VALUE(cnanim);
	dimX=INTEGER_VALUE(cdimX);
	dimQ=INTEGER_VALUE(cdimQ);

	varE=NUMERIC_VALUE(cvarE);

  	PROTECT(cW=AS_NUMERIC(cW));
	W=NUMERIC_POINTER(cW);

    PROTECT(cWWdiag=AS_NUMERIC(cWWdiag));
    WWdiag=NUMERIC_POINTER(cWWdiag); 
    
    PROTECT(cWWdiagadd=AS_NUMERIC(cWWdiagadd));
    WWdiagadd=NUMERIC_POINTER(cWWdiagadd);

    PROTECT(ctheta=duplicate(AS_NUMERIC(ctheta)));
    theta=NUMERIC_POINTER(ctheta);  
	
    PROTECT(cycorr=duplicate(AS_NUMERIC(cycorr)));
    ycorr=NUMERIC_POINTER(cycorr);

    PROTECT(cvarq=AS_NUMERIC(cvarq));
    varq=NUMERIC_POINTER(cvarq);
	

    xj=(double *) R_alloc(nanim,sizeof(double));
    newtheta=(double *) R_alloc(dimX+nSNP+dimQ,sizeof(double));

        for(j=0; j<(dimX+nSNP+dimQ); j++)
	{


                rhs=0;
		for(i=0; i<nanim; i++)
		{	
			xj[i]=W[i+(j)*nanim];
			ycorr[i]=ycorr[i]+theta[j]*xj[i];
			rhs+=xj[i]*ycorr[i];
		}

		lhs = (WWdiag[j] + WWdiagadd[j]);
		newtheta[j] = (rhs+WWdiag[j]*theta[j])/lhs + sqrt(varE/lhs)*rnorm(0,1) ;
		theta[j] = newtheta[j];

		for(i=0; i<nanim; i++)
		{
			ycorr[i]=ycorr[i]-theta[j]*xj[i];
		}
		
	}  
		
        // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 2));
	// attaching theta vector to list:
	SET_VECTOR_ELT(list, 0, ctheta);
	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 1, cycorr);

  	PutRNGstate();

  	UNPROTECT(7);
  	return(list);
	
}
