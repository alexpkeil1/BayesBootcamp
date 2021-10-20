bayes.sim <- function(m=50, beta1=1, beta2=1, beta3=1, beta4=1, beta5=1, std=1,
                       c12=.9, c13=.1, c14=.1, c15=.1, c23=.1, c24=.1, c25=.1, c34=.1, c35=.1, c45=.1) {

  #set seed, establish matricies to store main results

  #means<-stddevs<-medians<-rhat<-LB<-UB<-NULL

  #create mean/sd and variance/covariance matrix for 5 normally distributed exposures
  mu <- c(0,0,0,0,0)
  stddev <- c(1,3,1,5,3)

  v <- matrix(c(1, c12, c13, c14, c15,
                c12, 1, c23, c24, c25,
                c13, c23, 1, c34, c35,
                c14, c24, c34, 1, c45,
                c15, c25, c35, c45, 1),
              ncol=5)

  covar <- stddev %*% t(stddev) * v

  #Now, create the exposure matrix
  data1 <- MASS::mvrnorm(m,mu,Sigma=covar, empirical=F)
  colnames(data1)<-c('V1','V2','V3','V4','V5')

  #simulate normally distributed outcome

  y <- rnorm(m,beta1*data1[,1]+beta2*data1[,2]+beta3*data1[,3]+beta4*data1[,4]+beta5*data1[,5],sd=std)
  #create dataset that JAGS can understand
  d <- as.data.frame(cbind(y, data1))

  return(d)

} # END HERE


