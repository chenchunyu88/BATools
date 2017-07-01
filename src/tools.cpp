#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector set_ZD(NumericMatrix Z,NumericVector D) {
  int  nanim = Z.nrow(),nSNP=Z.ncol();
  NumericMatrix ZD(nanim,nSNP);
  
  for (int i = 0; i < nanim; i++)
    for (int j = 0; j < nSNP; j++) 
      ZD(i,j)=Z(i,j)*D[j]; 
  return ZD;
}