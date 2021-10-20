######################################################################################################################
# EPICOH 2021 workshop: The Brothers Bayes’ Boisterous Bootcamp
# Program: bayes_shortcourse_ex1_indpriors.R
# Languages: R, JAGS
# Date: Thursday, October 21, 2021
# Project: Exercise 1
# Tasks: Fit a Bayesian linear regression model and interpret output
# Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
######################################################################################################################

# Some preliminaries to set up your R environment
OS = Sys.info()[["sysname"]]

if(OS=="Windows"){
  base = "F:/Downloads/BayesBootcamp"
}
if(OS %in% c("Darwin", "Linux")){
  base = path.expand("~/Downloads/BayesBootcamp")
}
setwd(base)
source("R/bayes_shortcourse_simdata.R")

# simulate data
set.seed(123)
data1 <- bayes.sim(m=250, beta1=.3, beta2=.3, beta3=.3, beta4=.3, beta5=0, std=1,
                   c12=.9, c13=.9, c14=.8, c15=.1, c23=.8, c24=.8, c25=.1, c34=.9, c35=.1, c45=.1)



#####################################################################################
##### Below is a Bayesian model with a distinct prior for each exposure.
##### These priors are set to be non-informative, as indicated by the large variance
##### This approach would be used if there were an informative prior for each
##### individual exposure of interest.
#####################################################################################

b.ind <- function() {
  for (n in 1:length(y)) {            ## FIRST, we specify a typical regression model.
    y[n] ~ dnorm(mu[n], y.tau)

    mu[n] <- b0 + b1*V1[n] + b2*V2[n] + b3*V3[n] + b4*V4[n] + b5*V5[n]
  }
  # NEXT come the priors.
  # NOTE: JAGS (and winBUGS) specify the inverse of the variance (precision matrix),
  # so below, 0.01 indicates a variance of 100.
  b1 ~ dnorm(0,0.01)
  b2 ~ dnorm(0,0.01)
  b3 ~ dnorm(0,0.01)
  b4 ~ dnorm(0,0.01)
  b5 ~ dnorm(0,0.01)

  b0 ~ dnorm(0,0.01)
  y.tau ~ dgamma(0.01, 0.01)
  y.sigma <- 1/sqrt(y.tau)
}


## Below we specify the parameters to monitor. These parameters will be summarized
## in the basic output and will be accessible via a dataset that we create to
## assess convergence

b.ind.parms <- c('b1','b2','b3','b4','b5')



set.seed(1) #MCMC generally is not deterministic, so use a seed for reproducibility (but not testing)
bayes.m1 <- R2jags::jags(data=data1,parameters.to.save=b.ind.parms,
                 n.chains=3, n.iter=5000, n.burnin=2500,n.thin=1,
                 model.file=b.ind)
print(bayes.m1)

## The code below transforms the output from the JAGS procedure to an MCMC object
## This allows us to assess convergence of parameters.

m1.output <- coda::as.mcmc(bayes.m1)

## Below is code to generate the primary convergence assessment tools discussed in lecture.

coda::autocorr.plot(m1.output[,'b1']) #this will generate a plot for each chain (here, 3)
coda::densplot(m1.output[,'b1'])
coda::traceplot(m1.output[,'b1'])


coda::HPDinterval(m1.output[,'b1'])
