library(nimble)
rm(list=ls())
set.seed(123)

dataset <- read.csv("simulated_data.csv")
NMarkers <- 20
y <- dataset[,1]
x <- as.matrix(dataset[,2:(NMarkers+1)])
N <- length(y)

# Set MCMC parameters
M.burnin <- 10000
M <- 5000
n.thin <- 5

myBUGScode <- nimbleCode({
  alpha ~ dnorm(log(10),1)                   # Intercept 
  tau ~ dgamma(1.0E-4,1.0E-4)
  for(j in 1:NMarkers)   {
    beta[j] ~ dnorm(0,1)                # Regression coefficient 
  }
  for(i in 1:N) {
    er[i] ~ dnorm(0, tau)
    mu[i] <- alpha + inprod(beta[], X[i,]) + er[i]
    lambda[i] <- exp(mu[i])
    y[i] ~ dpois(lambda[i])             # Likelihood 
  }
})

constants <- list(N = N, NMarkers=NMarkers)
dimensions = list(beta = NMarkers,
                  Ind =  NMarkers,
                  lambda =  N,
                  X = c(N,NMarkers),
                  er = N,
                  mu = N)

myModel <- nimbleModel(myBUGScode, 
                       constants = constants, 
                       dimensions = dimensions)

myModel$setData(list(y = dataset[,1], X = as.matrix(dataset[,-1])))
myModel$setInits(list(alpha=1, tau = 0.5, beta = matrix(0,nrow = 1, ncol = NMarkers), er = rnorm(N,0,1)))

myMCMC <- buildMCMC(myModel)
compiled <- compileNimble(myModel, myMCMC)

compiled$myMCMC$run(M.burnin + M*n.thin)
samples <- as.matrix(compiled$myMCMC$mvSamples)

thinnedsample <- samples[M.burnin + seq(from = 1,by = n.thin, 
                             to = M*n.thin),]
colMeans(thinnedsample)
plot(thinnedsample[,2])
plot(density(thinnedsample[,2]))

library(coda)
fit <- mcmc(thinnedsample)
s <- summary(fit, start=M.burnin+1)
barplot(s$statistics[grep('beta', rownames(s$statistics)),'Mean'],
        main = expression(E(beta~'|'~y)))

load('true.Rda') # populates 'theta' vector
library(tidyverse)
varnames = rownames(s$statistics)[grep('beta', rownames(s$statistics))]
Ebeta = s$statistics[varnames, 'Mean']
quantbeta = s$quantiles[varnames, c('2.5%','97.5%')]
posterior <- data.frame(variable=varnames, beta=as.factor(1:20), Ebeta, quantbeta, true.theta=theta)
ggplot(posterior, aes(x=beta, y=Ebeta, ymin=X2.5., ymax=X97.5.)) +
  geom_point() +
  geom_linerange() +
  geom_point(aes(y=true.theta), color='red') +
  ggtitle('Ordinary regression model', subtitle='Posterior estimate vs true value') +
  labs(caption="Source: Matt's simulated data and BUGS model") +
  geom_text(x=11, y=0.2, label='True value', color='red') +
  geom_text(x=6, y=0.4, label='Posterior 95% CI', color='black') +
  theme_bw(15)
ggsave('poisson_regression.pdf', width=8, height = 6)

save(thinnedsample, posterior, file='poisson-multiple.Rda')
