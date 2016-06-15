# BATools

This package implements Bayesian model using  Markov Chain Monte Carlo (MCMC) and expectation maximization (EM) algorithm for various Bayesian models.

The window based approach presented at ICGQ5 will be tested and uploaded in the next few day. Contact Chunyu Chen at chenchunyu88@gmail.com for more information. 

The package is under developement, but it is functional.

To install the package you can follow the 
steps (recommended: use Rstudio to run this): 

`install.packages("devtools")` 

`library(devtools)` 

`devtools::install_github("hadley/devtools")` 

`devtools::install_github("chenchunyu88/BATools",build_vignettes=T)`

Or, you can download it as an .zip file and use `R CMD INSTALL BATools` to install it in the command line.

A quick introduction: `vignette("BATools")`

Help can be accessed by typing: `help(BATools)`, `help(bafit)`

Example codes can be obtained by typing: `demo(vignette_demo)` and it's located in the demo folder

###Current working models:
- rrBLUP based on REML and MCMC
- BayesA based on EM and MCMC
- SSVS (BayesC) based on EM and MCMC
- BayesB based on MCMC

###Future commits:
1) ONGOING: 

- improve documentation to expand to other capabilities

- add optimized anteBayesA/B to BayesM function

- add optimized IW/CD-BayesA/B to the package

2) GWA capabilities based on EM-BayesA/SSVS

3) implement multi-core computation

4) write complete testing code



