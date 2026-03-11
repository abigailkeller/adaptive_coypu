library(nimble)
library(CARBayesdata)
library(sp)
library(spdep)

# read in data
data(GGHB.IZ)
data(respiratorydata)

respiratorydata_spatial <- merge(x = GGHB.IZ, y = respiratorydata, 
                                 by = "IZ", all.x = FALSE)
W.nb <- poly2nb(respiratorydata_spatial, 
                row.names =  rownames(respiratorydata_spatial))
## Determine neighborhood/adjacency information needed for neighborhood-based CAR model
nbInfo <- nb2WB(W.nb)

# A vector of indices indicating which regions are neighbors of which.
nbInfo$adj

nregions <- nrow(respiratorydata_spatial)

code <- nimbleCode({
  # priors
  beta ~ dnorm(0, sd = 100)
  sigma ~ dunif(0, 100)   # prior for variance components based on Gelman (2006)
  tau <- 1 / sigma^2
  # latent process
  s[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau, zero_mean = 0)
  # likelihood
  for(i in 1:N) {
    log(lambda[i]) <- log(expected[i]) + beta*x[i] + s[i]
    y[i] ~ dpois(lambda[i])
  }
})

x <- respiratorydata_spatial$incomedep
x <- x - mean(x)  # center for improved MCMC performance

set.seed(1)

constants <- list(N = nregions, L = length(nbInfo$adj), 
                  adj = nbInfo$adj, weights = nbInfo$weights, num = nbInfo$num,
                  x = x, expected = respiratorydata_spatial$expected)
data <- list(y = respiratorydata_spatial$observed)
inits <- list(beta = 0, sigma = 1, s = rnorm(nregions))

model <- nimbleModel(code, constants = constants, data = data, inits = inits)