set.seed(1)
n <- 200
p <- 20

alpha <- log(10)
sig_e <- 0.75

x <- matrix(  rnorm(n = n*p, mean = 0, sd = 1), 
              nrow = n, ncol = p)

e <- rnorm(n = n, mean = 0, sd = sig_e)

(theta <- 0.3 + 0.3*(1:20/10.5 - 1))

nu <- alpha + x%*%theta + e
lambda <- exp(nu)
y <- rpois(n, lambda)

dataset <- data.frame(y=y, x=x)
write.csv(x = dataset, file = "simulated_data.csv", row.names = F)
