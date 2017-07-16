#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector set_ZD(NumericMatrix Z,NumericVector D) {
  int  nanim = Z.nrow(),nSNP=Z.ncol();
  NumericMatrix ZD(nanim,nSNP);
  
  for (int i = 0; i < nanim; i++)
    for (int j = 0; j < nSNP; j++) 
      ZD(i,j)=Z(i,j)*D[j]; 
  return ZD;
}

// [[Rcpp::export]]
List get_pvalue(NumericVector U,NumericVector V,NumericMatrix G,NumericMatrix Z,NumericMatrix X){

mat z = as<mat>(Z);
mat g = as<mat>(G);
mat x = as<mat>(X);
vec u = as<vec>(U);
vec v = as<vec>(V);
int n=z.n_rows,m=z.n_cols;
mat res(2,m),Cuu(n,n),ginv(n,n),I(n,n),va(n,n),AiVAiZ(n,m);
double l=v(0)/v(1);
vec tmp(m);
res.fill(0);
I.eye();
ginv=inv(g);
Rprintf("Starting GWA and computing Cuu for the animal model \n");
Cuu=ginv*l+I-x*inv(x.t()*x)*x.t();
Cuu=inv(Cuu)*v(0);


Rprintf("Preparing to Compute diagonal of Caa \n");
va=g*v(1)-Cuu;

AiVAiZ=ginv*va*ginv*z;

tmp=z.t()*(ginv*u);

Rprintf("Computing diagonal of Caa \n");
for(int j=0;j<m;j++){
  res(0,j)=tmp(j);
  for(int i=0;i<n;i++){
    //<<i<<" "<<j<<endl;
    res(1,j)+=z(i,j)*AiVAiZ(i,j);
  }
  
}
Rprintf("Completed! \n");
return List::create(Named("gw") = res,
                    Named("AiVAiZ") = AiVAiZ);
}