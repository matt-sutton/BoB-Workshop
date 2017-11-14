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



