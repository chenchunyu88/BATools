#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP BayesAslope(SEXP cnSNP, SEXP cnanim, SEXP cW, SEXP cWWdiag, SEXP cycorr, SEXP cvarq, SEXP cvarE, SEXP cSNPeff)
{
	double *xj, *W, *WWdiag, *ycorr, *varq, *SNPeff;
	double rhs,varE,lhs;
	int j,i, nSNP, nanim;
        SEXP list;
	
  	GetRNGstate();
	
	nSNP=INTEGER_VALUE(cnSNP);
	nanim=INTEGER_VALUE(cnanim);
	varE=NUMERIC_VALUE(cvarE); 
		
  	PROTECT(cW=AS_NUMERIC(cW));
	W=NUMERIC_POINTER(cW); 
        
       PROTECT(cWWdiag=AS_NUMERIC(cWWdiag));
       WWdiag=NUMERIC_POINTER(cWWdiag);  
	
       PROTECT(cycorr=duplicate(AS_NUMERIC(cycorr)));
       ycorr=NUMERIC_POINTER(cycorr);

       PROTECT(cvarq=AS_NUMERIC(cvarq));
       varq=NUMERIC_POINTER(cvarq);
	
       PROTECT(cSNPeff=AS_NUMERIC(cSNPeff));
       SNPeff=NUMERIC_POINTER(cSNPeff);

      xj=(double *) R_alloc(nanim,sizeof(double));
		
	for(j=0; j<nSNP; j++)
	{		
		rhs=0;
		for(i=0; i<nanim; i++)
		{	
			xj[i]=W[i+j*nanim];
			ycorr[i]=ycorr[i]+SNPeff[j]*xj[i];
			rhs+=xj[i]*ycorr[i];
		}
		
		rhs=rhs/varE;
		lhs=WWdiag[j]/varE+1.0/varq[j];
		SNPeff[j]=rhs/lhs + sqrt(1.0/lhs)*rnorm(0,1) ;

		for(i=0; i<nanim; i++)
		{
			ycorr[i]=ycorr[i]-SNPeff[j]*xj[i];
		}
		
	}  
		
    // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 4));

	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 0, cycorr);
	
	SET_VECTOR_ELT(list, 1, cvarq);
	
	SET_VECTOR_ELT(list, 2, cvarE);
	// attaching SNPeff vector to list:
	SET_VECTOR_ELT(list, 3, cSNPeff);
  	PutRNGstate();

  	UNPROTECT(6);
  	return(list);
	
}

SEXP BayesAintercept(SEXP cnSNP, SEXP cnanim, SEXP cZ, SEXP cycorr, SEXP cvarq, SEXP cvarE, SEXP cSNPeff, SEXP cphi, SEXP cd)
{
	double *xj, *Z, *ycorr, *varq, *SNPeff, *phi, *d;
	double rhs,diag,varE,lhs;
	int j,i, nSNP, nanim;
       SEXP list;
	
  	GetRNGstate();
	
	nSNP=INTEGER_VALUE(cnSNP);
	nanim=INTEGER_VALUE(cnanim);
	varE=NUMERIC_VALUE(cvarE); 
		
  	PROTECT(cZ=AS_NUMERIC(cZ));
	Z=NUMERIC_POINTER(cZ);  
	
       PROTECT(cycorr=duplicate(AS_NUMERIC(cycorr)));
       ycorr=NUMERIC_POINTER(cycorr);

       PROTECT(cvarq=AS_NUMERIC(cvarq));
       varq=NUMERIC_POINTER(cvarq);
	
       PROTECT(cSNPeff=AS_NUMERIC(cSNPeff));
       SNPeff=NUMERIC_POINTER(cSNPeff);

       PROTECT(cphi=AS_NUMERIC(cphi));
       phi=NUMERIC_POINTER(cphi);

       PROTECT(cd=AS_NUMERIC(cd));
       d=NUMERIC_POINTER(cd);

       xj=(double *) R_alloc(nanim,sizeof(double));

			
	for(j=0; j<nSNP; j++)
	{		
		rhs=0;
              diag=0;
		for(i=0; i<nanim; i++)
		{	
 			xj[i]=Z[i+j*nanim]+Z[i+j*nanim]*d[i]*phi[j];
			ycorr[i]=ycorr[i]+SNPeff[j]*xj[i];
			rhs+=xj[i]*ycorr[i];
                     diag+=xj[i]*xj[i];
		}
		
		rhs=rhs/varE;
              lhs=diag/varE+1.0/varq[j];
		SNPeff[j]=rhs/lhs + sqrt(1.0/lhs)*rnorm(0,1) ;

		for(i=0; i<nanim; i++)
		{
			ycorr[i]=ycorr[i]-SNPeff[j]*xj[i];
		}
		
	}  
		
    // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 4));

	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 0, cycorr);
	
	SET_VECTOR_ELT(list, 1, cvarq);
	
	SET_VECTOR_ELT(list, 2, cvarE);
	// attaching SNPeff vector to list:
	SET_VECTOR_ELT(list, 3, cSNPeff);
  	PutRNGstate();

  	UNPROTECT(7);
  	return(list);
	
}

SEXP BayesAassociate(SEXP cnSNP, SEXP cnanim, SEXP cZ, SEXP cycorr, SEXP cvarE, SEXP cSNPeff, SEXP cbb, SEXP cd, SEXP cml0, SEXP csl02)
{
	double *xj, *Z, *ycorr, *SNPeff, *bb, *d;
	double rhs,diag,varE,lhs,ml0,sl02;
	int j,i, nSNP, nanim;
       SEXP list;
	
  	GetRNGstate();
	
	nSNP=INTEGER_VALUE(cnSNP);
	nanim=INTEGER_VALUE(cnanim);
	varE=NUMERIC_VALUE(cvarE); 
	ml0=NUMERIC_VALUE(cml0); 
	sl02=NUMERIC_VALUE(csl02); 
	
  	PROTECT(cZ=AS_NUMERIC(cZ));
	Z=NUMERIC_POINTER(cZ);  
	
        PROTECT(cbb=duplicate(AS_NUMERIC(cbb)));
        bb=NUMERIC_POINTER(cbb); 

        PROTECT(cycorr=duplicate(AS_NUMERIC(cycorr)));
        ycorr=NUMERIC_POINTER(cycorr);
	
        PROTECT(cSNPeff=AS_NUMERIC(cSNPeff));
        SNPeff=NUMERIC_POINTER(cSNPeff);

        PROTECT(cd=AS_NUMERIC(cd));
        d=NUMERIC_POINTER(cd);

        xj=(double *) R_alloc(nanim,sizeof(double));
	
	for(j=0; j<nSNP; j++)
	{		
		rhs=0;
              diag=0;
		for(i=0; i<nanim; i++)
		{	
 			xj[i]=Z[i+j*nanim]*d[i]*bb[j];
			ycorr[i]=ycorr[i]+SNPeff[j]*xj[i];
			rhs+=xj[i]*ycorr[i];
                     diag+=xj[i]*xj[i];
		}
		
		rhs=rhs/varE+ml0/sl02;
              lhs=diag/varE+1.0/sl02;
		SNPeff[j]=rhs/lhs + sqrt(1.0/lhs)*rnorm(0,1) ;

		for(i=0; i<nanim; i++)
		{
			ycorr[i]=ycorr[i]-SNPeff[j]*xj[i];
		}
		
	}  
		
    // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 3));

	// attaching ycorr vector to list:
	SET_VECTOR_ELT(list, 0, cycorr);
	
	SET_VECTOR_ELT(list, 1, cvarE);
	// attaching SNPeff vector to list:
	SET_VECTOR_ELT(list, 2, cSNPeff);
  	PutRNGstate();

  	UNPROTECT(6);
  	return(list);
	
}


