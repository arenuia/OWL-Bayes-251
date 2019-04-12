library(invgamma)

owlS01ST4 <- read.csv("rawOWLdata-S01ST4.csv", stringsAsFactors = FALSE)
owlS02ST1 <- read.csv("rawOWLdata-S02ST1.csv", stringsAsFactors = FALSE)

season1 <- owlS01ST4$UltHoldTime
season2 <- owlS02ST1$UltHoldTime
set.seed(411)
# Here we use the mean time to charge an ult as the parameter. Both mu and sigma^2
# are unknown so we must use the Normal-Normal-Inverse Gamma distribution. Our prior
# beliefs on the distributions of mu are that the mean time to charge an ult in 
# season 1 is 60 seconds with a standard deviation of 15 seconds. The mean time to
# charge an ult in season 2 is 50 seconds with a standard deviation of 12 seconds.
# These values were based on our personal experiences in the game. Our prior beliefs
# on the distributions of sigma^2 for both seasons are that a = b = 1. We picked 
# these values as we had no real prior information on this parameter and decided
# this would be a good place to start.



# Prior Values
n1 <- length(season1)
n2 <- length(season2)
m1 <- 60
m2 <- 50
v1 <- 15^2
v2 <- 12^2
a <- 1
b <- 1
# Prior Distributions
# Season 1: mu ~ Normal(m=60, v=225), sigma^2 ~ IG(a = 1, b = 1)
# Season 2: mu ~ Normal(m=50, v=144), sigma^2 ~ IG(a = 1, b = 1)


par(mfrow=c(1,2))

x1 <- seq(0, 120, by = .1)
prior1 <- dnorm(x1, m1, sqrt(v1))
plot(x1, prior1, type = "l", main = expression(mu %~% "N(60, 15)"), xlab = "θ", ylab = "π(θ)")

x2 <- seq(0, 15, by = .1)
prior2 <- dinvgamma(x2, a, b)
plot(x2, prior2, type = "l", main = expression(sigma ^ 2 %~% "IG(1, 1)"), xlab = "θ", ylab = "π(θ)")

par(mfrow=c(1,2))

x3 <- seq(10, 90, by = .1)
prior3 <- dnorm(x3, m2, sqrt(v2))
plot(x3, prior3, type = "l", main = expression(mu %~% "N(50, 12)"), xlab = "θ", ylab = "π(θ)")

x4 <- seq(0, 15, by = .1)
prior4 <- dinvgamma(x4, a, b)
plot(x4, prior4, type = "l", main = expression(sigma ^ 2 %~% "IG(1, 1)"), xlab = "θ", ylab = "π(θ)")


#################
#    Season 1   #
#################

# Populating Posterior Distribution
#----------------------------------------------------------------------#
nreps <- 10000
# create vectors to store draws
mus.1 <- rep(NA, nreps)
sigma2s.1 <- rep(NA, nreps)
# Set initial values
mus.1[1] <- mean(season1)
sigma2s.1[1] <- var(season1)
#----------------------------------------------------------------------#


# Sampling from complete conditionals
#----------------------------------------------------------------------#
for (i in 2:nreps){
  # update vstar and mstar
  temp.vstar <- 1/(n1/sigma2s.1[i-1] + 1/v1)
  temp.mstar <- temp.vstar*(sum(season1)/sigma2s.1[i-1] + m1/v1)
  # update mus.1
  mus.1[i] <- rnorm(1, temp.mstar, sqrt(temp.vstar))
  # update astar and bstar
  temp.astar <- astar <- a + n1/2
  temp.bstar <- b + .5*sum((season1 - mus.1[i])^2)
  # update sigma2s.1
  sigma2s.1[i] <- 1/rgamma(1, temp.astar, temp.bstar)
}
#----------------------------------------------------------------------#


# Trace plot of the mu and sigma2 draws
#----------------------------------------------------------------------#
# We found that there was no need to remove any burnout as the distribution was 
# fairly accurate from the first observation due to our initial values

plot(mus.1, type='l' , ylab=expression(mu) , 
     main=expression(paste("Trace plot for ",mu)))
plot(sigma2s.1, type='l', ylab=expression(sigma^2), 
     main=expression(paste("Trace plot for ",sigma^2)))
# Estimated Autocorrelations
acf(mus.1, main = expression(paste("Series ",mu, " Season 1")))
acf(sigma2s.1, main = expression(paste("Series ",sigma^2, " Season 1")))
#----------------------------------------------------------------------#


# Posterior Summary Statistics
#----------------------------------------------------------------------#
mean(mus.1) # est. posterior mean of mu
# 41.74092
mean(sigma2s.1) # est. posterior mean of sigma2
# 183.6061
sd(mus.1) # est. posterior std. dev of mu
# 1.20876
sd(sigma2s.1) # est. posterior std. dev of sigma2
# 23.61869
quantile(mus.1, c(.025, .975)) # Est. CI
# 2.5%      97.5%
# 39.39558  44.13925  
quantile(sigma2s.1, c(.025, .975)) # Est. CI
# 2.5%      97.5%
# 143.0798  235.5031
#----------------------------------------------------------------------#


# Plot Estimated posterior Densities
#----------------------------------------------------------------------#
plot(density(mus.1), 
     main = expression(paste("Density Plot for ", mu, " Season 1")))
plot(density(sigma2s.1), 
     main = expression(paste("Density Plot for ", sigma^2, " Season 1")))
#----------------------------------------------------------------------#




#################
#    Season 2   #
#################

# Populating Posterior Distribution
#----------------------------------------------------------------------#
nreps <- 10000
# create vectors to store draws
mus.2 <- rep(NA, nreps)
sigma2s.2 <- rep(NA, nreps)
# Set initial values
mus.2[1] <- mean(season2)
sigma2s.2[1] <- var(season2)
#----------------------------------------------------------------------#


# Sampling from complete conditionals
#----------------------------------------------------------------------#
for (i in 2:nreps){
  # update vstar and mstar
  temp.vstar <- 1/(n2/sigma2s.2[i-1] + 1/v2)
  temp.mstar <- temp.vstar*(sum(season2)/sigma2s.2[i-1] + m2/v2)
  # update mus.2
  mus.2[i] <- rnorm(1, temp.mstar, sqrt(temp.vstar))
  # update astar and bstar
  temp.astar <- astar <- a + n2/2
  temp.bstar <- b + .5*sum((season2 - mus.2[i])^2)
  # update sigma2s.2
  sigma2s.2[i] <- 1/rgamma(1, temp.astar, temp.bstar)
}
#----------------------------------------------------------------------#


# Trace plot of the mu and sigma2 draws
#----------------------------------------------------------------------#
# We found that there was no need to remove any burnout as the distribution was 
# fairly accurate from the first observation due to our initial values
plot(mus.2, type='l' , ylab=expression(mu) , 
     main=expression(paste("Trace plot for ",mu)))
plot(sigma2s.2, type='l', ylab=expression(sigma^2), 
     main=expression(paste("Trace plot for ",sigma^2)))
# Estimated Autocorrelations
acf(mus.2, main = expression(paste("Series ",mu, " Season 2")))
acf(sigma2s.2, main = expression(paste("Series ",sigma^2, " Season 2")))
#----------------------------------------------------------------------#


# Posterior Summary Statistics
#----------------------------------------------------------------------#
mean(mus.2) # est. posterior mean of mu
# 30.49075
mean(sigma2s.2) # est. posterior mean of sigma2
# 105.9311
sd(mus.2) # est. posterior std. dev of mu
# 0.8254836
sd(sigma2s.2) # est. posterior std. dev of sigma2
# 12.23851
quantile(mus.2, c(.025, .975)) # Est. CI
# 2.5%      97.5%
# 28.85884  32.08187 
quantile(sigma2s.2, c(.025, .975)) # Est. CI
# 2.5%      97.5%
# 84.61365  132.52618 
#----------------------------------------------------------------------#


# Plot Estimated posterior Densities
#----------------------------------------------------------------------#
plot(density(mus.2), 
     main = expression(paste("Density Plot for ", mu, " Season 2")))
plot(density(sigma2s.2), 
     main = expression(paste("Density Plot for ", sigma^2, " Season 2")))
#----------------------------------------------------------------------#



##################################
#    Comparing the Two Seasons   #
##################################

par(mfrow=c(1,1))
mu1_mu2 <- mus.1 - mus.2
# Create a histogram of the mean for the first season minus 
# the mean for the second season
hist(mu1_mu2, main = "1st Season Mean - 2nd Season Mean", 
     xlab = expression(mu[1]-mu[2]))
mean(mu1_mu2 > 0)
# 1
quantile(mu1_mu2, c(0.025, 0.975))
# 2.5%     97.5% 
# 8.408521 14.152023

# As this interval does not include 0, we can conclude that the mean ult hold 
# time for the first season is larger than the mean ult hold time for the second 
# season.



##############################
#   Prior versus Posterior   #
##############################

# Season 1
# Prior:      Normal(m1 = 60, v1= 15^2)
# Posterior:  Normal(m* = 41.74, v* = 183.61)
par(mfrow=c(1,1))
thetas <- seq(-750, 750, length=2000)
plot(thetas, dnorm(thetas, 41.74, 183.61), col="black", lwd=2.5, 
     main="Prior and Posterior Distributions Season 1", xlab=expression(theta), 
     lty=4, type="l", ylab=expression(paste(pi, "(", theta, ")", sep="")))
lines(thetas, dnorm(thetas, 60, 225), col="blue", lwd=2.5, lty=2)
legend("topright", c("Posterior", "Prior"), col=c("black","blue"), lwd=2.5, 
       lty=c(4,2))


# Season 2
# Prior:      Normal(m2 = 50, v2= 12^2)
# Posterior:  Normal(m* = 30.49, v* = 105.93)
thetas <- seq(-500, 500, length=2000)
plot(thetas, dnorm(thetas, 30.49, 105.93), col="black", lwd=2.5, 
     main="Prior and Posterior Distributions Season 2", xlab=expression(theta), 
     lty=4, type="l", ylab=expression(paste(pi, "(", theta, ")", sep="")))
lines(thetas, dnorm(thetas, 50, 144), col="blue", lwd=2.5, lty=2)
legend("topright", c("Posterior", "Prior"), col=c("black","blue"), lwd=2.5, 
       lty=c(4,2))
