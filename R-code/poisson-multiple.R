N <- 200; NMarkers <- 20

myBUGScode <- nimbleCode({
  alpha ~ dnorm(log(10),1)                   # Intercept 
  tau ~ dgamma(1.0E-4,1.0E-4)
  for(j in 1:NMarkers)   {
    beta [j] ~ dnorm(0,1)                # Regression coefficient 
  }
  for(i in 1:N) {
    er ~ dnorm(0, tau)
    mu <- alpha + inprod(beta[], X[i,]) + er
    log(lambda[i]) <- mu
    y[i] ~ dpois(lambda[i])             # Likelihood 
  }
})

constants <- list(N = N, NMarkers=NMarkers)
dataset <- read.csv("simulated_data.csv")
dimensions = list(beta = NMarkers,
                  Ind =  NMarkers,
                  lambda =  N,
                  X = c(N,NMarkers))

myModel <- nimbleModel(myBUGScode, 
                       constants = constants, 
                       dimensions = dimensions)

myModel$setData(list(y = dataset[,1], X = as.matrix(dataset[,-1])))
myModel$setInits(list(alpha=2, tau = 1, beta = matrix(0,nrow = 1, ncol = NMarkers)))

myMCMC <- buildMCMC(myModel)
compiled <- compileNimble(myModel, myMCMC)

compiled$myMCMC$run(10000)
samples <- as.matrix(compiled$myMCMC$mvSamples)

plot(density(samples[-c(1:500),1]))

colMeans(samples[-c(1:500),])
