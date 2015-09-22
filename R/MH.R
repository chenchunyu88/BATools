metropDf=function (def, scalea, G1, SNPeff2, priordf = 1, cdef = 0.5, max = 200,pi_math=3.141593)
{
          def0=def
          logdef0 = log(def0)
          logdef1 = rtnorm(n=1,mean=logdef0,sd=cdef,upper=log(max))   # proposal log df 
          def1 = exp(logdef1)
          lfullcd0 = G1*(lgamma((def0+1)/2)-lgamma(0.5*def0)+ 0.5*log(1/(pi_math*def0*scalea)))
          lfullcd1 = G1*(lgamma((def1+1)/2)-lgamma(0.5*def1)+ 0.5*log(1/(pi_math*def1*scalea)))
          # THE FOLLOWING MIGHT BE A LITTLE DANGEROUS WITH ZERO SNP EFFECTS..BUT IT APPEARS TO WORK
          lfullcd0 = lfullcd0 - sum((0.5*def0+0.5)*log(1+SNPeff2/(def0*scalea)));
          lfullcd1 = lfullcd1 - sum((0.5*def1+0.5)*log(1+SNPeff2/(def1*scalea)));
          lfullcd0 = lfullcd0 + log(priordf)- 2*log(priordf+def0)+logdef0   # prior + jacobian                                      
          lfullcd1 = lfullcd1 + log(priordf)- 2*log(priordf+def1)+logdef1   # prior + jacobian                          
          alpha1 = exp(lfullcd1-lfullcd0)
          alpha1 = min(alpha1,1)
          if (alpha1 < 1)
          { 
            u = runif(1)
            if (u<alpha1) {def = def1}
            else def = def0
          }
          else {def = def1;}
         
    	return(list(def=def,alpha1=alpha1))
}


metropScale=function (def, scalea, G1, SNPeff2, cscalea,alpha_scale, beta_scale, pi_math=3.141593)
{
         scalea0=scalea
         logscalea0 = log(scalea0)
         logscalea1 = rnorm(1,logscalea0,cscalea)   # proposal log df 
         scalea1 = exp(logscalea1)
         lfullcd0 = G1*(lgamma((def+1)/2)-lgamma(0.5*def)+ 0.5*log(1/(pi_math*scalea0*def)))
         lfullcd1 = G1*(lgamma((def+1)/2)-lgamma(0.5*def)+ 0.5*log(1/(pi_math*scalea1*def)))
         lfullcd0 = lfullcd0 - sum((0.5*def+0.5)*log(1+SNPeff2/(scalea0*def)));
         lfullcd1 = lfullcd1 - sum((0.5*def+0.5)*log(1+SNPeff2/(scalea1*def)));
         # SEE CHANGE BELOW!
         lfullcd0 = lfullcd0 + (alpha_scale-1)*log(scalea0) - beta_scale*scalea0 +logscalea0   # prior + jacobian                                      
         lfullcd1 = lfullcd1 + (alpha_scale-1)*log(scalea1) - beta_scale*scalea1 +logscalea1   # prior + jacobian                          
         alpha2 = exp(lfullcd1-lfullcd0)
         alpha2 = min(alpha2,1)
         if (alpha2 < 1)
          { 
            u = runif(1)
            if (u<alpha2) {scalea = scalea1}
            else scalea = scalea0
          }
          else {scalea = scalea1;}
    	
	return(list(scalea=scalea,alpha2=alpha2))
}
