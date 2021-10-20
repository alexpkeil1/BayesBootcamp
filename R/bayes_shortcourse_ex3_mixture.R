######################################################################################################################
# EPICOH 2021 workshop: The Brothers Bayes’ Boisterous Bootcamp
# Program: bayes_shortcourse_ex3_mixture.R
# Languages: R, JAGS
# Date: Thursday, October 21, 2021
# Project: Exercise 3
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



################################################################################
##### Below is code for a Bayesian mixture prior
##### This is a bit more complicated than the former models
##### However, I find it to be very applicable to environmental health research
################################################################################

b.mix <- function() {
  for (n in 1:length(y)) {
    y[n] ~ dnorm(mu[n], y.tau)

    mu[n] <- b0 + b[1]*V1[n] + b[2]*V2[n] + b[3]*V3[n] + b[4]*V4[n] + b[5]*V5[n]
  }
  # Below is the mixture prior. Here we allow 4 groups, one of which is a spiked null.
  # To increase/decrease the number of potential groups, change t[1]-t[4] to desired number.
  #   For example, if you want 5 possible distributions, that means addind a t[5].
  #   Also note that this means you need to adapt the stick breaking prior code below.
  # To increase/decrease the number of exposures in the model, change b[1]-b[5] to desired number.
  # Note: t[k[i]] should correspond to the number of exposures, which may be obvious from the code below.

  t[1]~dnorm(0,.01)
  t[2]~dnorm(0,.01)
  t[3]<-0
  t[4]~dnorm(0,.01)

  for (v in 1:5) {
    k[v]~dcat(p[])
  }

  b[1]<-t[k[1]]
  b[2]<-t[k[2]]
  b[3]<-t[k[3]]
  b[4]<-t[k[4]]
  b[5]<-t[k[5]]

  # Below is the machinery behind the mixture prior, known as the 'stick breaking prior' (described elsewhere)
  # The key points are as follow:
  # 1: The prior allows the parameters to explore multiple distributions.
  # 2: Mostly you should NOT change anything in this section of code.
  # 3: However, if you want to increase the number of distribution options,
  #    then you need to change the maximum number of sections, which here is 4.
  #    Also note that
  #    As an example, at the end of this script I comment out a section with 5 possible distributions.

  p[1]<-r[1]
  for (k in 1:3){
    r[k]~dbeta(1,1.5)
  }
  for (k in 2:3){
    p[k]<-r[k]*(1-r[k-1])*p[k-1]/r[k-1]
  }

  p[4]<-1-sum(p[1:3])

  # This ends the mixture prior

  # Here, we specify the priors for other variables not subject to the mixture distrbution.

  b0 ~ dnorm(0,0.01)
  y.tau ~ dgamma(0.01, 0.01)
  y.sigma <- 1/sqrt(y.tau)
}

## Below we specify the parameters to monitor. These parameters will be summarized
## in the basic output and will be accessible via a dataset that we create to
## assess convergence

b.mix.parms <- c('b[1:5]','k[1:5]')



set.seed(1) #MCMC generally is not deterministic, so use a seed for reproducibility (but not testing)
bayes.m3 <- R2jags::jags(data=data1,parameters.to.save=b.mix.parms,
                 n.chains=3, n.iter=5000, n.burnin=2500,n.thin=1,
                 model.file=b.mix)
print(bayes.m3)

## The code below transforms the output from the JAGS procedure to an MCMC object
## This allows us to assess convergence of parameters.

m3.output <- coda::as.mcmc(bayes.m3)

## Below is code to generate the primary convergence assessment tools discussed in lecture.

coda::autocorr.plot(m3.output[,'b[1]']) #this will generate a plot for each chain (here, 3)
coda::densplot(m3.output[,'b[1]'])
coda::traceplot(m3.output[,'b[1]'])


coda::HPDinterval(m3.output[,'b[1]'])
