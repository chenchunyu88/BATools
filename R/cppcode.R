
library('Rcpp')
library('inline')
library("RcppArmadillo")

cppFunction('NumericVector set_ZD(NumericMatrix Z,NumericVector D) {
  int  nanim = Z.nrow(),nSNP=Z.ncol();
  NumericMatrix ZD(nanim,nSNP);

  for (int i = 0; i < nanim; i++)
    for (int j = 0; j < nSNP; j++) 
  	  ZD(i,j)=Z(i,j)*D[j]; 
  return ZD;
}')

cppFunction('NumericVector get_new_C10(NumericVector ZD,NumericMatrix AiVAiZD,NumericMatrix C10,double c) {
int len=C10.nrow()-1,nanim=AiVAiZD.nrow();
  NumericMatrix C10new(C10);
            
  for (int i = 0; i < len; i++){
     for (int j = 0; j < nanim; j++){
         C10new(len,i)-=ZD[j]*AiVAiZD(j,i); 
     }
   C10new(i,len)=C10new(len,i); 
  }
C10new(len,len)=c;
  return C10new;
}')

rcpp_inc <- '
using namespace Rcpp;
using namespace arma;
'

src <- '
mat z = as<mat>(Z);
vec w = as<vec>(a);
vec g=as<vec>(b);
int lenw=w.size(),p1=0,p2=0;
double nanim=z.n_rows;
vec varw(lenw),gw(nanim);



varw.fill(0);
gw.fill(0);

for(int k=0;k<lenw;k++){
	p2=p1+w(k)-1;
	if(p1==p2) gw=z(span::all,p1)*g(p1);
		
	else gw=z.cols(p1,p2)*g.subvec(p1,p2);
	
		
	varw(k)=accu(square(gw))/nanim-pow((accu(gw)/nanim),2.0);
	p1=p2+1;
}

return(wrap(varw));
'
fc <- cxxfunction(signature(a="numeric", Z="numeric", b="numeric"), src, plugin='RcppArmadillo', rcpp_inc)


calc.varw<-function(a,Z,b)
{
	as.numeric(fc(a,Z,b))
}
