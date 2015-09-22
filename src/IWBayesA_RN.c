#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <math.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>
#include <stdio.h>
#define N 2

void Multiply_Matrices_2x2_CRC(double *C, SEXP cA, double *B) //multiply 2x2 matrix
{
	double *A; //might be tricky here
	PROTECT(cA=AS_NUMERIC(cA));
	A=NUMERIC_POINTER(cA);	
	C[0] = A[0] * B[0] + A[1] * B[2];
	C[1] = A[0] * B[1] + A[1] * B[3];
	C[2] = A[2] * B[0] + A[3] * B[2];
	C[3] = A[2] * B[1] + A[3] * B[3];
	UNPROTECT(1);	
}

void Multiply_Matrices_2x2(double *C, double *A, double *B) //multiply 2x2 matrix
{
	
   C[0] = A[0] * B[0] + A[1] * B[2];
   C[1] = A[0] * B[1] + A[1] * B[3];
   C[2] = A[2] * B[0] + A[3] * B[2];
   C[3] = A[2] * B[1] + A[3] * B[3];
}

void M_inverse_2x2_sys(double *C, double *A) //inverse 2x2 sysmetric matrix
{
	double determ;
	determ = A[0]*A[3] - A[1]*A[2];
	
	C[0] = A[3]/determ;
	C[3] = A[0]/determ;
	C[1] = -A[1]/determ;
	C[2] = -A[1]/determ;
}

void M_inverse_2x2_sys_CR(double *C, SEXP cA) //inverse 2x2 sysmetric matrix
{
	double determ,*A;
	PROTECT(cA=AS_NUMERIC(cA));
    A=NUMERIC_POINTER(cA);
	//Rprintf("inside M_inverse_2x2_sys_CR matrix %f %f %f %f \n",A[0],A[1],A[2],A[3]);
	determ = A[0]*A[3] - A[1]*A[2];
	C[0] = A[3]/determ;
	C[3] = A[0]/determ;
	C[1] = -A[1]/determ;
	C[2] = -A[1]/determ;
	//Rprintf("C inside M_inverse_2x2_sys_CR matrix %f %f %f %f \n",C[0],C[1],C[2],C[3]);
	UNPROTECT(1);	
}



SEXP M_inverse_2x2_sys_RC(SEXP cC, double *A) //inverse 2x2 sysmetric matrix
{
	double *C,determ;
	PROTECT(cC=AS_NUMERIC(cC));
       C = NUMERIC_POINTER(cC);
 	determ = A[0]*A[3]-A[1]*A[2];
	C[0] = A[3]/determ;
	C[3] = A[0]/determ;
	C[1] = -A[1]/determ;
	C[2] = -A[1]/determ;
	UNPROTECT(1);
	return cC;
}


void chol_2x2(double *C, double *A) //Cholesky Decomposition of 2x2 matrix
{
	C[0] = sqrt(A[0]);
	C[1] = A[1]/sqrt(A[0]);
	C[3] = sqrt(A[3]-R_pow_di(C[1],2));
	C[2] = 0;
}

void Self_Transpose_Multiply(double *C, double *A, int nrows, int ncols ) //calculate t(A) %*% A in R
{
   double *pA;
   double *p_A = A;
   double *pB;
   double *pC;
   double *p_C = C;
   double *p_CC;
   int i,j,k;

   for (i = 0; i < ncols; p_C += ncols, i++) {
      p_A = A + i;
      p_CC = p_C + i;
      pC = p_CC;
      for (j = i; j < ncols; pC++, p_CC += ncols, j++) {
         pA = p_A;
         pB = A + j;
         *pC = 0.0;
         for (k = 0; k < nrows; pA += ncols, pB += ncols, k++) 
            *pC += *pA * *pB;
         *p_CC = *pC;
      }
   }
}

SEXP Add_Matrices_ToExsiting_2x2_RC(SEXP cC, double *A) 
{
	double *C;
	PROTECT(cC=AS_NUMERIC(cC));
    C=NUMERIC_POINTER(cC);
	C[0] += A[0] ;
	C[1] += A[1];
	C[2] += A[2];
	C[3] += A[3];
	UNPROTECT(1);
	return cC;
}

SEXP Add_Matrices_ToExsiting_2x2_RR(SEXP cC, SEXP cA) 
{
	double *C,*A;
	PROTECT(cC=AS_NUMERIC(cC));
    C=NUMERIC_POINTER(cC);
	PROTECT(cA=AS_NUMERIC(cA));
    A=NUMERIC_POINTER(cA);
	C[0] += A[0] ;
	C[1] += A[1];
	C[2] += A[2];
	C[3] += A[3];
	UNPROTECT(2);
	return cC;
}

void Add_Matrices_ToExsiting_2x2(double *C, double *A) 
{
   C[0] += A[0] ;
   C[1] += A[1];
   C[2] += A[2];
   C[3] += A[3];
}

void Add_Matrices_2x2(double *C, double *A, double *B) 
{
   C[0] = A[0] + B[0];
   C[1] = A[1] + B[1];
   C[2] = A[2] + B[2];
   C[3] = A[3] + B[3];
}

SEXP Add_Matrices_2x2_RCC(SEXP cC, double *A, double *B) 
{
	double *C;
	PROTECT(cC=AS_NUMERIC(cC));
    C=NUMERIC_POINTER(cC);
	C[0] = A[0] + B[0];
	C[1] = A[1] + B[1];
	C[2] = A[2] + B[2];
	C[3] = A[3] + B[3];
	UNPROTECT(1);
	return cC;
	
}


 void Multiply_2x2_Matrix_by_Scalar(double *A, double x) //Multiply the 2??2 matrix A by the scalar x, i.e. A ?? xA.
{
   A[0] *= x;
   A[1] *= x;
   A[2] *= x;
   A[3] *= x;
}


void Fill_Matrix_with_Scalar(double *A, double x, int nrows, int ncols)
{
   int i,j;
 
   for (i = 0; i < nrows; i++)
      for ( j = 0; j < ncols; j++) *A++ = x;
}

void Multiply_Matrices(double *C, double *A, int nrows, int ncols,double *B, int mcols) 
{
   double *pA = A;
   double *pB;
   double *p_B;
   //double *pC = C;
   int i,j,k;

   for (i = 0; i < nrows; A += ncols, i++) 
      for (p_B = B, j = 0; j < mcols; C++, p_B++, j++) {
         pB = p_B;
         pA = A;
         *C = 0.0; 
         for (k = 0; k < ncols; pA++, pB += mcols, k++) 
            *C += *pA * *pB;
      }
}


//sample variance-covariance matrix for all loci

SEXP sample_vars(SEXP cb,SEXP cS, SEXP cdef, SEXP cSig, SEXP cinvsigma, SEXP cvar, SEXP cprint,SEXP citer)
{
	//FILE *my_file;
	int nmarkers=length(cb)/2,i,j,p,count=0,iprint,iter;
	double CC[2][2]={{0.0, 0.0},{0.0,0.0}},*b,*Sig,*invsigma,*var,temp1=0,temp2=0,*temp,v,def,*s;
	double Z[2][2]={{0.0, 0.0},{0.0,0.0}},wish[2][2]={{0.0, 0.0},{0.0,0.0}},sigma[2][2]={{0.0, 0.0},{0.0,0.0}},tmp[2][2]={{0.0, 0.0},{0.0,0.0}},t[2][2]={{0.0, 0.0},{0.0,0.0}}, invS[2][2]={{0.0, 0.0},{0.0,0.0}},determ=0; //pp[2][2]={{0.0, 0.0},{0.0,0.0}},,mean[2]={0.0,0.0}
	SEXP list,ctemp;
	GetRNGstate();
	
	def=NUMERIC_VALUE(cdef);


	PROTECT(ctemp=allocVector(REALSXP,2));
    	temp=NUMERIC_POINTER(ctemp);
	PROTECT(cS=AS_NUMERIC(cS));
    	s=NUMERIC_POINTER(cS);

	PROTECT(cb=AS_NUMERIC(cb));
    	b=NUMERIC_POINTER(cb);

	PROTECT(cSig=AS_NUMERIC(cSig));
    	Sig=NUMERIC_POINTER(cSig);

	PROTECT(cinvsigma=AS_NUMERIC(cinvsigma));
    	invsigma=NUMERIC_POINTER(cinvsigma);

	PROTECT(cvar=AS_NUMERIC(cvar));
    	var=NUMERIC_POINTER(cvar);
        iprint=INTEGER_VALUE(cprint);
	iter=INTEGER_VALUE(citer);
	for(j=0;j<(2*nmarkers);j++)
	{
		if((j % 2)!=1)
		{
			s[0]=b[j]*b[j];
			s[3]=b[j+1]*b[j+1];
			s[1]=b[j]*b[j+1];
			s[2]=s[1];
			v=1+def;
			for(i=0;i<4;i++)
			{
				s[i]+=Sig[i]*(def-3);
			}
			p=2;
			M_inverse_2x2_sys_CR((double *)invS, cS);
                        chol_2x2(&CC[0][0],&invS[0][0]);                       //Cholesky Decomposition
			Z[0][0]=sqrt(rchisq(v));
			Z[1][1]=sqrt(rchisq((v - p + 1)));
			Z[0][1]=rnorm(0,1);
			Multiply_Matrices_2x2(&tmp[0][0],&Z[0][0], &CC[0][0]);	
			wish[0][0] = tmp[0][0]*tmp[0][0] + tmp[1][0]*tmp[1][0];
			wish[0][1] = tmp[0][0]*tmp[0][1] + tmp[1][0]*tmp[1][1];
			wish[1][0] = wish[0][1];
			wish[1][1] = tmp[0][1]*tmp[0][1] + tmp[1][1]*tmp[1][1];	
			M_inverse_2x2_sys(&sigma[0][0], &wish[0][0]);
			
			determ = sigma[0][0]*sigma[1][1] - sigma[0][1]*sigma[0][1];
			
			temp1 += log(determ);

			var[count]   = sigma[0][0];
			var[count+1] = sigma[1][1];
			var[count+2] = sigma[0][1];
			count += 3;

			invsigma[j+nmarkers*2*j] = sigma[1][1]/determ;
			invsigma[j+1+nmarkers*2*(j+1)] = sigma[0][0]/determ;
			invsigma[j+nmarkers*2*(j+1)]   =- sigma[0][1]/determ;
			invsigma[j+1+nmarkers*2*j]     =- sigma[0][1]/determ;
			
			t[0][0] = invsigma[j+nmarkers*2*j];
			t[1][1] = invsigma[j+1+nmarkers*2*(j+1)];
			t[0][1] = invsigma[j+nmarkers*2*(j+1)];
			t[1][0] = invsigma[j+1+nmarkers*2*j];
			
			
			Multiply_Matrices_2x2_CRC(&tmp[0][0],cSig, &t[0][0]);
			temp2 = temp2 + tmp[0][0] + tmp[1][1];
			
			
		}
	}
	
	PROTECT(list = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(list, 0, cS);
	SET_VECTOR_ELT(list, 1, cvar);
	temp[0] = temp1;
	temp[1] = temp2;
	SET_VECTOR_ELT(list, 2, ctemp);
	SET_VECTOR_ELT(list, 3, cinvsigma);
		
	PutRNGstate();

  	UNPROTECT(7);
  	return(list);
	
}

SEXP sample_snp(SEXP cnmarkers,SEXP cnanim, SEXP cx, SEXP cb,SEXP cxpx,SEXP cinvsigma, SEXP cycorr, SEXP cvare,  SEXP cprint)
{
	double *ycorr, *x,*b,*xpx,*invsigma,*xj,*xjp1;
	double rhs[2],invLhs[2][2]={{0.0, 0.0},{0.0,0.0}},lhs[2][2]={{0.0 ,0.0},{0.0 ,0.0}},tmp[2][2]={{0.0, 0.0},{0.0,0.0}},pp[2][2]={{0.0, 0.0},{0.0,0.0}};//invS[2][2]={{0.0, 0.0},{0.0,0.0}},,t[2][2]={{0.0, 0.0},{0.0,0.0}}
	double mean[2]={0.0,0.0},L[N][N]={{0.0 ,0.0},{0.0 ,0.0}},m[2]={0.0,0.0},vare;
	int i,j,nmarkers,nanim,iprint;
	SEXP list;
    	
  	GetRNGstate();
	
	nmarkers = INTEGER_VALUE(cnmarkers);
	nanim = INTEGER_VALUE(cnanim);

	vare = NUMERIC_VALUE(cvare);
	
	double xx[2][nanim],tx[2][nanim];

	PROTECT(cx=AS_NUMERIC(cx));
    	x=NUMERIC_POINTER(cx);

	PROTECT(cb=AS_NUMERIC(cb));
    	b=NUMERIC_POINTER(cb);
	
	PROTECT(cxpx=AS_NUMERIC(cxpx));
    	xpx=NUMERIC_POINTER(cxpx);
	
	PROTECT(cinvsigma=AS_NUMERIC(cinvsigma));
    	invsigma=NUMERIC_POINTER(cinvsigma);

	PROTECT(cycorr=AS_NUMERIC(cycorr));
    	ycorr=NUMERIC_POINTER(cycorr);
        iprint=INTEGER_VALUE(cprint);	

	xj=(double *) R_alloc(nanim,sizeof(double));	
    	xjp1=(double *) R_alloc(nanim,sizeof(double));	
	
	/*sample genetic effects*/
	Fill_Matrix_with_Scalar(&tmp[0][0],0,2,2);
	tmp[0][0]=1/vare;
	tmp[1][1]=1/vare;
	
	for(j=0;j<(2*nmarkers);j++)
	{
		if((j % 2)!=1)
		{
			Fill_Matrix_with_Scalar(&xx[0][0],0,2,nanim);
			Fill_Matrix_with_Scalar(&tx[0][0],0,2,nanim);
			
                        for(i=0;i<nanim;i++)
			{
				xj[i]=x[i+j*nanim];
				xjp1[i]=x[i+(j+1)*nanim];
				xx[0][i]=xj[i];
				xx[1][i]=xjp1[i];
				ycorr[i]=ycorr[i]+xj[i]*b[j]+xjp1[i]*b[j+1];
				
			}
			//rhs
			Multiply_Matrices(&tx[0][0],&tmp[0][0],2,2,&xx[0][0],nanim);
			Multiply_Matrices(&rhs[0],&tx[0][0],2,nanim,&ycorr[0],1);
			
			//lhs
			pp[0][0] = xpx[j+nmarkers*2*j];
			pp[0][1] = xpx[j+nmarkers*2*(j+1)];
			pp[1][0] = xpx[j+1+nmarkers*2*j];	
			pp[1][1] = xpx[j+1+nmarkers*2*(j+1)];
			
			Multiply_Matrices_2x2(&lhs[0][0],&tmp[0][0],&pp[0][0]);	
			lhs[0][0] += invsigma[j+nmarkers*2*j];
			lhs[0][1] += invsigma[j+nmarkers*2*(j+1)];
			lhs[1][0] += invsigma[j+1+nmarkers*2*j];
			lhs[1][1] += invsigma[j+1+nmarkers*2*(j+1)];
		
			
			M_inverse_2x2_sys(&invLhs[0][0],&lhs[0][0]);
			
			mean[0] = rhs[0]*invLhs[0][0] + rhs[1]*invLhs[1][0];
			mean[1] = rhs[0]*invLhs[0][1] + rhs[1]*invLhs[1][1];
			
			chol_2x2(&L[0][0],&invLhs[0][0]);
			
			
			m[0]=rnorm(0,1);
			m[1]=rnorm(0,1);
			//m[0]=0.1000258;
			//m[1]=-0.1037045;
			//Rprintf("m matrix %f %f\n",m[0],m[1]);
			
			b[j]   = m[0]*L[0][0] + L[1][0]*m[1] + mean[0];
			b[j+1] = m[0]*L[0][1] + L[1][1]*m[1] + mean[1];
			
	//		if(iprint==1) Rprintf("b12 matrix %E %E marker %d \n",b[j],b[j+1],j);
			
			
			for(i=0;i<nanim;i++)
			{
				ycorr[i]=ycorr[i]-xj[i]*b[j]-xjp1[i]*b[j+1];
				//Rprintf("ycorr %f %d\n",ycorr[i],i);				
			}		
		}
	}
	
    	// Creating a list with 7 vector elements:
	PROTECT(list = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(list, 0, cb);
	SET_VECTOR_ELT(list, 1, cycorr);
	PutRNGstate();

  	UNPROTECT(6);
  	return(list);
	
}

SEXP sample_scale(SEXP cnmarkers,SEXP cS,SEXP cSig, SEXP cdef, SEXP cinvsigma, SEXP cdf0, SEXP cSig0,SEXP cprint)
{
	double *s,*invsigma,*Sig;
	double CC[2][2]={{0.0, 0.0},{0.0,0.0}},Z[2][2]={{0.0, 0.0},{0.0,0.0}},wish[2][2]={{0.0, 0.0},{0.0,0.0}},tmp[2][2]={{0.0, 0.0},{0.0,0.0}},t[2][2]={{0.0, 0.0},{0.0,0.0}}, invS[2][2]={{0.0, 0.0},{0.0,0.0}},Sum[2][2]={{0.0, 0.0},{0.0,0.0}};

	double def,df0,v;
	int j,nmarkers,p,iprint;
	SEXP list;
	
	GetRNGstate();

	nmarkers=INTEGER_VALUE(cnmarkers);
	
	def=NUMERIC_VALUE(cdef);
	df0=NUMERIC_VALUE(cdf0);

    	PROTECT(cS=AS_NUMERIC(cS));
    	s=NUMERIC_POINTER(cS);

	PROTECT(cSig=AS_NUMERIC(cSig));
    	Sig=NUMERIC_POINTER(cSig);

	PROTECT(cinvsigma=AS_NUMERIC(cinvsigma));
    	invsigma=NUMERIC_POINTER(cinvsigma);
        iprint = INTEGER_VALUE(cprint);
	//PROTECT(cstoreSig=AS_NUMERIC(cstoreSig));
    	//storeSig=NUMERIC_POINTER(cstoreSig);
	
    	//sampling scales
  	Fill_Matrix_with_Scalar(&Sum[0][0],0,2,2); //get zero matrix
	
	for(j=0;j<(2*nmarkers);j++)
	{
		if((j % 2)!=1)
		{
			t[0][0]=invsigma[j+nmarkers*2*j];
			t[1][1]=invsigma[j+1+nmarkers*2*(j+1)];
			t[0][1]=invsigma[j+nmarkers*2*(j+1)];
			t[1][0]=invsigma[j+1+nmarkers*2*j];
			Add_Matrices_ToExsiting_2x2(&Sum[0][0],&t[0][0]);
		}
        }
	
	v=df0+def*nmarkers;
	Multiply_2x2_Matrix_by_Scalar(&Sum[0][0],(def-3));

	M_inverse_2x2_sys_CR(&t[0][0],cSig0);
	Multiply_2x2_Matrix_by_Scalar(&t[0][0],df0);
	cS=Add_Matrices_2x2_RCC(cS,&Sum[0][0],&t[0][0]);
	Fill_Matrix_with_Scalar(&tmp[0][0],0,2,2); //Fill matrix with zero
	
	p=2;		
						
	M_inverse_2x2_sys_CR((double *)invS, cS);
	s[0]=invS[0][0];
	s[1]=invS[1][0];
	s[2]=invS[1][0];
	s[3]=invS[1][1];
			
			
	chol_2x2(&CC[0][0],&invS[0][0]); //Cholesky Decomposition
	
	Z[0][0] = sqrt(rchisq(v));
	Z[1][1] = sqrt(rchisq((v - p + 1)));
	Z[0][1] = rnorm(0,1);
	
	Multiply_Matrices_2x2(&tmp[0][0],&Z[0][0], &CC[0][0]);	
	Self_Transpose_Multiply(&wish[0][0], &tmp[0][0], 2, 2);	//crossprod	

	Sig[0] = wish[0][0];
	Sig[1] = wish[0][1];
	Sig[2] = wish[1][0];
	Sig[3] = wish[1][1];
	
	PROTECT(list = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(list, 0, cS);

	SET_VECTOR_ELT(list, 1, cSig);
		
	PutRNGstate();

  	UNPROTECT(4);
  	return(list);
}


