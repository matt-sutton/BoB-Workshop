
Bayes on the Beach Workshop 2017
================================

Workshop Aims:
--------------

-   Review some Bayesian Variable Selection methods (O’Hara and Sillanpaa suggested material)
-   Some coding in R (Visualization | Coding | Methods)
-   Analysis of a simple model with variable selection.

### Pages to visit

The Workshop is now complete. Thanks to all the participants: Sam, Grace, Anton, Saras, Anthony, Nick Golding, Leah, Andres for a great workshop!

-   [An example solution](Example-Solution/)
-   [Resources for the workshop](Slides-And-Resources/)
-   [Sour dough (Variable selection visualisation)](Visualisation/)
-   [Exploration and modelling](Modelling/)

The Task
--------

### Data Model:

<img src="datamodel.png" width="290px" style="display: block; margin: auto;" />

R code for simulation from data model (O’Hara and Sillanpaa):

``` r
set.seed(1)
## Data for the simulations see description under eq (2):

n <- 200; p <- 20
alpha <- log(10)
sig_e <- 0.75

x <- matrix(  rnorm(n = n*p, mean = 0, sd = 1), 
              nrow = n, ncol = p)

e <- rnorm(n = n, mean = 0, sd = sig_e)
theta <- 0.3 + 0.3*(1:20/10.5 - 1)

## Simulate for this dataset
mu <- alpha + x%*%theta + e
lambda <- exp(mu)
y <- rpois(n, lambda)

dataset <- data.frame(y=y, x=x)
#write.csv(x = dataset, file = "simulated_data.csv", row.names = F)
```

The regression parameters *θ*<sub>*j*</sub> have a mean of 0.3, rough range (0, 0.6). While not exact zeros, the aim is to investigate the inference (eg probability of inclusion) and sparse estimation in Bayesian framework.

``` r
#plot(1:length(theta), theta, xlab = "Regression parameter", ylab = "Theta")
```

### Modelling for variable seleciton

Recommended priors for the model (eq. 3, 4 of O’Hara and Sillanpaa)

<img src="priors.png" width="370px" style="display: block; margin: auto;" />

Each *θ*<sub>*j*</sub> has a spike-slab prior with spike at zero and slab elsewhere. Introduce an auxilary indicator to indicate presence of vat variable and auxilary *β*<sub>*j*</sub> to indicate effect size.
