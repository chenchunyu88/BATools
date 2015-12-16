#include <RcppArmadillo.h>
#include <Rcpp.h> 
using namespace Rcpp;
extern "C" SEXP pcg_c( SEXP As, SEXP bs, SEXP Minvs, SEXP maxiters, SEXP tols) {
	
    Rcpp::NumericVector br(bs);                 // creates Rcpp vector from SEXP
    Rcpp::NumericMatrix Ar(As);                 // creates Rcpp matrix from SEXP
	Rcpp::NumericMatrix Minvr(Minvs);
    int n = Ar.nrow(), k = Ar.ncol();

	int maxiter=as<int>(maxiters);

	double tol=as<double>(tols),a,beta;

    arma::mat A(Ar.begin(), n, k, false);       // reuses memory and avoids extra copy

	arma::mat Minv(Minvr.begin(),n,n,false);

    arma::colvec b(br.begin(), br.size(), false);

	arma::colvec x(br.size()),Ap(br.size()),r1(br.size()),z1(br.size());
	
	x.fill(0);
	arma::colvec r=b-A*x;
	arma::colvec z=Minv*r;
	arma::colvec p(z);
	int iter=0;
	double sumr2=dot(r,r);
	while(sumr2>tol && iter<maxiter){
		iter++;
		Ap=A*p;
		a=dot(r,z)/dot(p,Ap);
		x+=a*p;
		r1=r-a*Ap;
		z1=Minv*r1;
		beta=dot(z1,r1)/dot(z,r);
		p=z1+beta*p;
		z=z1;
		r=r1;
		sumr2=dot(r,r);
	}
	return Rcpp::wrap(x);
}

/*
extern "C" SEXP set_coeff( SEXP XXs, SEXP XZs, SEXP ZXs, SEXP ZZs) {
	
    Rcpp::NumericMatrix XXr(XXs);                 // creates Rcpp matrix from SEXP
	Rcpp::NumericMatrix XZr(XZs);                 // creates Rcpp matrix from SEXP
	Rcpp::NumericMatrix ZXr(ZXs);                 // creates Rcpp matrix from SEXP
	Rcpp::NumericMatrix ZZr(ZZs);                 // creates Rcpp matrix from SEXP
    int dimW = XXr.nrow()+ZZr.nrow(),dimX=XXr.nrow(),dimZ=ZZr.nrow();
    
	arma::mat coeff=zeros<mat>(dimW,dimW);
	

	coeff.submat(0,0,dimX-1,dimX-1)=XXr;
	//coeff.submat(0,dimX,dimX,dimW-1)=XZr;
	//coeff.submat(dimX,0,dimW-1,dimX)=ZXr;
	//coeff.submat(dimX,dimX,dimW-1,dimW-1)=ZZr;

	return Rcpp::wrap(coeff);
}*/


	
