# BATools

This package implements Bayesian model using  Markov Chain Monte Carlo (MCMC) and expectation maximization (EM) algorithm for various Bayesian models.

The package is under developement, but it is functional.

To install the package you can follow the 
steps (recommended: use Rstudio to run this): 

`install.packages("devtools")` 

`library(devtools)` 

`devtools::install_github("hadley/devtools")` 

`devtools::install_github("chenchunyu88/BATools",build_vignettes=T)`

A quick introduction: `vignette("BATools")`

###Current working models:
- 1) rrBLUP based on REML and MCMC
- 2) BayesA based on EM and MCMC
- 3) BayesC based on REML and MCMC
- 4) BayesB based on MCMC

###Future commits:
1) ONGOING: 

- improve documentation to expand to other capabilities

- add optimized anteBayesA/B to BayesM function

- add optimized IW/CD-BayesA/B to the package

2) GWA capabilities based on EM-BayesA/C

3) implement multi-core computation


