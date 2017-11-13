#install.packages("nimble")
library(nimble)

myBUGScode <- nimbleCode({
  mu ~ dnorm(0, sd = 100) ## uninformative prior
  sigma ~ dunif(0, 100)
  for(i in 1:10) y[i] ~ dnorm(mu, sd = sigma)
})

myModel <- nimbleModel(myBUGScode)

myData <- rnorm(10, mean = 2, sd = 5)
myModel$setData(list(y = myData))
myModel$setInits(list(mu = 0, sigma = 1))

myMCMC <- buildMCMC(myModel)
compiled <- compileNimble(myModel, myMCMC)

compiled$myMCMC$run(10000)

samples <- as.matrix(compiled$myMCMC$mvSamples)
plot(density(samples[,'mu']))
plot(density(samples[,'sigma']))
