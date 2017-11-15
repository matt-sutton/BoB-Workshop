rm(list=ls())
set.seed(123)
library(nimble)

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
    betaT[j] ~ dnorm(0,1)                # Regression coefficient 
    Ind[j] ~ dbern(PInd)       # Indicator 
    beta [j] <- Ind[j]*betaT[j]
  }
  for(i in 1:N) {
    er[i] ~ dnorm(0, tau)
    mu[i] <- alpha + inprod(beta[], X[i,]) + er[i]
    lambda[i] <- exp(mu[i])
    y[i] ~ dpois(lambda[i])             # Likelihood 
  }
})

constants <- list(N = N, NMarkers=NMarkers, PInd = 0.2)
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
myModel$setInits(list(alpha=1, tau = 0.5, betaT = matrix(0,nrow = 1, ncol = NMarkers), er = rnorm(N,0,1)))

myMCMC <- buildMCMC(myModel)
compiled <- compileNimble(myModel, myMCMC)

compiled$myMCMC$run(M.burnin + M*n.thin)
samples <- as.matrix(compiled$myMCMC$mvSamples)

thinnedsample <- samples[M.burnin + seq(from = 1,by = n.thin, 
                             to = M*n.thin),]
colMeans(thinnedsample)
plot(thinnedsample[,2])
plot(density(thinnedsample[,2]))

cs <- colnames(thinnedsample)
beta.names <- cs[grep('beta', cs)]
ind.names  <- cs[grep('Ind', cs)]

thetas <- thinnedsample[,beta.names] * thinnedsample[,ind.names]

library(coda)
# fit <- mcmc(thinnedsample)
fit <- mcmc(thetas)
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
  ggtitle(label='Variable selection model', subtitle = "Posterior estimate vs true value") +
  labs(caption="Source: Matt's simulated data and BUGS model") +
  geom_text(x=3, y=0.2, label='True value', color='red') +
  geom_text(x=16, y=0.15, label='Posterior 95% CI', color='black') +
  xlab('theta') + ylab('Etheta') +
  theme_bw(15)
ggsave('poisson_regression_variable_selection.pdf', width=8, height = 6)
save(thetas, posterior, file='poisson-multiple_VS1.Rda')

#load('poisson-multiple_VS1.Rda')

# trace plots for variable inclusion
spikes <- thinnedsample[,ind.names] %>%
  as.tibble() %>%
  mutate(iter=1:nrow(thinnedsample)) %>%
  gather(key = 'variable', value='spike', -iter, factor_key=TRUE) %>%
  mutate(groupid = cumsum(c(1, abs(diff(spike))))) %>%
  group_by(groupid, variable, spike) %>%
  summarize(start=min(iter), end=max(iter))

ggplot(spikes, aes(x=start, y=variable, group=variable)) +
    geom_line() +
    geom_label(x=1500, y=5, label='Draws with zeroed variable') +
    ggtitle("'Spike' traceplot by variable") +
    ylab('Variable') + xlab('Iteration') +
    theme_bw(15)

# probability of inclusion
inc <- colMeans(thinnedsample[,ind.names])
inc.data <- data.frame(variable=as.factor(1:20), inc.prob = inc)
ggplot(inc.data, aes(x=variable, y=inc.prob)) +
  geom_bar(stat='identity') +
  ggtitle('Variable inclusion probability', subtitle="Spike-and-slab model (Matt's implementation") +
#  labs(caption="Source: Matt's simulated data and BUGS model") +
  xlab('Variable') + ylab('Inclusion probability') +
  theme_bw(15)
ggsave('inclusion_probability.pdf', width=8, height=6)
ggsave('inclusion_probability.png', width=6, height=6)
