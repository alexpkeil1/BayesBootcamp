######################################################################################################################
# EPICOH 2021 workshop: The Brothers Bayes’ Boisterous Bootcamp
# Program: bayes_shortcourse_ex2_hierarchical.R
# Languages: R, JAGS
# Date: Thursday, October 21, 2021
# Project: Exercise 2
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


#################################################################################
##### Below is the popular Bayesian Hierarchical model (a.k.a Semi or Empirical)
##### This model leverages the assumption that the effects of each exposure
##### may be partially exchangeable (i.e, share a distribution).
#################################################################################

b.shared <- function() {
  for (n in 1:length(y)) {
    y[n] ~ dnorm(mu[n], y.tau)

    mu[n] <- b0 + b[1]*V1[n] + b[2]*V2[n] + b[3]*V3[n] + b[4]*V4[n] + b[5]*V5[n]
  }
  # Below is the shared mean prior. It says that we think the effects of all exposures may be similar.
  # This prior does not preclude finding different effects for exposures when the data support them.

  for (j in 1:5) {                #here, 1:5 indicates 5 different exposures
    b[j] ~ dnorm(beta.m,beta.tau)
  }

  sigma.beta ~ dunif(0,100)
  beta.tau <- 1/(pow(sigma.beta,2))
  beta.m ~ dnorm(0,0.1)

  # This ends the shared mean prior

  # Below, we specify the priors for other variables not subject to the shared mean prior.

  b0 ~ dnorm(0,0.01)
  y.tau ~ dgamma(0.01, 0.01)
  y.sigma <- 1/sqrt(y.tau)
}


## Below we specify the parameters to monitor. These parameters will be summarized
## in the basic output and will be accessible via a dataset that we create to
## assess convergence

b.shared.parms <- c('b[1:5]', 'beta.m', 'sigma.beta') # note the (arbitrary) different specification of this from b1.ind.parms



set.seed(1) #MCMC generally is not deterministic, so use a seed for reproducibility (but not testing)
bayes.m2 <- R2jags::jags(data=data1,parameters.to.save=b.shared.parms,
                 n.chains=3, n.iter=5000, n.burnin=2500,n.thin=1,
                 model.file=b.shared)
print(bayes.m2)

## The code below transforms the output from the JAGS procedure to an MCMC object
## This allows us to assess convergence of parameters.

m2.output <- coda::as.mcmc(bayes.m2)

## Below is code to generate the primary convergence assessment tools discussed in lecture.

coda::autocorr.plot(m2.output[,'b[1]']) #this will generate a plot for each chain (here, 3)
coda::densplot(m2.output[,'b[1]'])
coda::traceplot(m2.output[,'b[1]'])


coda::HPDinterval(m2.output[,'b[1]'])
