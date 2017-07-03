# BATools

This package implements Bayesian model using  Markov Chain Monte Carlo (MCMC) and Max a Posterior (MAP) algorithm for various Bayesian models for whole genome prediction (WGP) and genome-wide association (GWA).

The data in the demo is explained in https://github.com/chenchunyu88/batoolsdata/blob/master/MSUPRP.ipynb.

### Warning: 
The demo example is just showing how the packages is work, please refer to these papers appear in `help(BATools)` for previous study results how the prediction accuarcy and GWA can be affected by different factors. In order for MCMC to convergence, change number of iteration and burn-In to 100,000 or greater depending on the size of your data.  

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

Example codes can be found in demo/ folder. 





