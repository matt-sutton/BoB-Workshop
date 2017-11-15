library(tidyverse)

dat <- read.csv('simulated_data.csv')

X <- as.matrix(dat[,-1])
N <- nrow(X)
p <- ncol(X)
Sigma.hat <- 1/N * t(X) %*% X
image(Sigma.hat[p:1,], main='Data covariance matrix')

library(BoomSpikeSlab)
set.seed(123)
niter = 1e4
y <- dat[,'y']
X.design <- cbind(alpha=1, X)
prior <- SpikeSlabPrior(x=X.design, y=y,
                        expected.model.size = 0,
                        prior.inclusion.probabilities = rep(0.1, 21)
)
fit <- poisson.spike(y ~ ., data = dat, niter = niter,
                     prior = prior)
summary(fit)
plot(fit)

thinnedsample <- fit$beta

# probability of inclusion
inc <- colMeans(thinnedsample != 0)[-1]
inc.data <- data.frame(variable=as.factor(1:20), inc.prob = inc)
ggplot(inc.data, aes(x=variable, y=inc.prob)) +
  geom_bar(stat='identity') +
  ggtitle('Variable inclusion probability', subtitle='Spike-and-slab model (BoomSpikeSlab package)') +
  xlab('Variable') + ylab('Inclusion probability') +
  theme_bw(15)
ggsave('package_inclusion_probability.pdf', width=6, height=6)
ggsave('package_inclusion_probability.png', width=6, height=6)
