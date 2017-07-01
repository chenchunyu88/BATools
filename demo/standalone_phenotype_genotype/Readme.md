##Demo code for using BATools directly with genotypes (Z), phenotypes (y) and map data:
The demo code performs whole genome regression as well as GWA at the same time. It implements windows based GWA if window size is sepcified. 

###Explanation of each file:
1. BA_MCMC.R			BayesA (t-distributed prior on marker effects) using MCMC
2. SSVS_MCMC.R			SSVS (stochastic search variable selection prior) using MCMC
3. BB_MCMC.R			BayesB (zero mixture of t-distributed prior) using MCMC
4. anteBA_MCMC.R 		anteBayesA (BayesA that account for association between adjacent SNPs) using MCMC
5. anteBB_MCMC.R 		anteBayesB (BayesB that account for association between adjacent SNPs) using MCMC
6. BA_EM.R				BayesA using EM algorithm
7. SSVS_EM.R			SSVS using EM algorithm
8. rrBLUP.R				rrBLUP (is equivalent to GBLUP when n is less than m) using REML
9. rrBLUP_MCMC.R		rrBLUP using MCMC