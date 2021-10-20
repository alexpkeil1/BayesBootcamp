######################################################################################################################
# EPICOH 2021 workshop: The Brothers Bayes’ Boisterous Bootcamp
# Program: bayes_shortcourse_ex4_bma_gcomp.R
# Languages: R, JAGS
# Date: Thursday, October 21, 2021
# Project: Exercise 4
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
##### Below is code for a Bayesian model selection/averaging (BMA) with g-computation
#####
##### This model implements Bayesian model selection using a stochastic search
##### where individual model coefficients are probabilistically included in the
##### linear model according to the "posterior inclusion probability."
##### This approach may be useful with high dimensional exposures or where model
##### form is uncertain and you wish to average the fit over many possible model
##### forms. This approach is called Bayesian model averaging because it averages
##### the model fit across all possible combinations of coefficients in the model.
#####
##### This approach also implements a very simple form of g-computation where
##### the expected outcome under no intervention on exposures is compared with
##### the expected outcome under an intervention to set all exposures to exactly
##### their observed median values. The "meandiff" parameter is equal to the
##### expected reduction in the average outcome following the intervention to
##### change all exposures to the median.
################################################################################

b.bma <- function() {
  meandiff <- mean(muint) - mean(mu) # mean difference comparing "as exposed" with "exposed at the median"
  for (n in 1:length(y)) {
    y[n] ~ dnorm(mu[n], y.tau)
    muint[n] <- b0 + b[1]*V1int + b[2]*V2int + b[3]*V3int + b[4]*V4int + b[5]*V5int
    mu[n] <- b0 + b[1]*V1[n] + b[2]*V2[n] + b[3]*V3[n] + b[4]*V4[n] + b[5]*V5[n]
  }

    # BMA here is implemented by including a 1/0 indicator variable for every
    # coefficient, where delta = 1 indicates the exposure is EXCLUDED from the model
    # at a given iteration of the MCMC sampler. The
    for(j in 1:5){
      # selection probability (same prior across all exposures)
      delta[j] ~ dbern(pi)
      # "posterior inclusion probability": probability that an individual exposure is selected into the model
      pip[j] <- 1-delta[j]
      # beta coefficient prior to selection
      bpr[j] ~ dnorm(bpr.mu, bpr.tau)
      # actual model coefficient used in model
      b[j] <- bpr[j]*(1-delta[j])
    }
    # prior probability of exclusion (mean = 0.5)
    pi ~ dbeta(1,1)

  # note no model selection used on intercept
  b0 ~ dnorm(0,0.01)
  # other priors/hyperpriors
  y.tau ~ dgamma(0.01, 0.01)
  # "psuedo-priors" on regression coefficients ("pseudo" used because the actual coefficient is subject to selection)
  bpr.mu ~ dnorm(0,0.01)
  bpr.tau ~ dgamma(0.01, 0.01)
}

## Below we specify the parameters to monitor. These parameters will be summarized
## in the basic output and will be accessible via a dataset that we create to
## assess convergence

b.bma.parms <- c('meandiff', 'b[1:5]','pip[1:5]')

# modify the data to use "intervention" values of exposure
# jags also accepts R lists as data
data_list = as.list(data1)
data_list$V1int = median(data1$V1)
data_list$V2int = median(data1$V2)
data_list$V3int = median(data1$V3)
data_list$V4int = median(data1$V4)
data_list$V5int = median(data1$V5)


set.seed(1) #MCMC generally is not deterministic, so use a seed for reproducibility (but not testing)
bayes.m4 <- R2jags::jags(data=data_list,parameters.to.save=b.bma.parms,
                 n.chains=3, n.iter=5000, n.burnin=2500,n.thin=1,
                 model.file=b.bma)
print(bayes.m4)

## The code below transforms the output from the JAGS procedure to an MCMC object
## This allows us to assess convergence of parameters.

m4.output <- coda::as.mcmc(bayes.m4)

## Below is code to generate the primary convergence assessment tools discussed in lecture.

coda::autocorr.plot(m4.output[,'b[1]']) #this will generate a plot for each chain (here, 3)
coda::densplot(m4.output[,'b[1]'])
coda::densplot(m4.output[,'meandiff'])
coda::traceplot(m4.output[,'b[1]'])
coda::traceplot(m4.output[,'meandiff']) # a "mixture effect"


coda::HPDinterval(m4.output[,'b[1]'])
