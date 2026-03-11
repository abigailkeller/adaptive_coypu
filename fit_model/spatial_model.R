library(sf)
library(tidyverse)
library(janitor)
library(lubridate)
library(nimble)
#library(ubms)
library(MCMCvis)
library(sp)
library(spdep)

# import data
raw_dat <- readRDS("data/coypus.rds")
temperature <- readRDS("data/temperature.rds") %>% 
  mutate(avg = rowMeans(cbind(tdec, tjan, tfeb, tmar), na.rm = TRUE))
communes_shp <- st_read("data/shp/communes.shp")

# reorder communes
communes_shp <- communes_shp[match(temperature$commune, communes_shp$dpts_cl), ]

# create df of commune names to join with temp
commune_convert <- data.frame(
  temp_name = temperature$commune,
  dat_name = c("baillargues", "candillargues", "entre_vignes", 
               "la_grande_motte", "lansargues", "lunel", "lunel_viel",
               "marsillargues", "mauguio", "perols", 
               "saint_genies_des_mourgues", "saint_just",
               "saint_nazaire_de_pezan", "saint_vincent_de_barbeyrargues",
               "valergues")
)

# plot communes
ggplot(data = communes_shp) +
  geom_sf(aes(fill = dpts_cl)) +
  labs(fill = "Commune") +
  theme_minimal()

# get 2022 data
dat <- raw_dat$removal_coypus %>%
  # filter(Month %in% c(2, 3, 4, 5, 6)) %>%
  pivot_wider(names_from = Month, 
              values_from = n) %>%
  filter(Year == 2022)

# add temp and raw data to spatial data
communes_shp$temp <- temperature$avg
communes_shp$remove_2022 <- rowSums(dat[, 3:14])

# ggplot(data = communes_shp) +
#   geom_sf(aes(fill = remove_2022)) +
#   labs(fill = "Total\nremoved") +
#   theme_minimal()

# get neighborhood structure
nb <- poly2nb(communes_shp)
nbInfo <- nb2WB(nb)

# adjust island
## add to baillargues
nbInfo$num[1] <- 2
nbInfo$adj <- c(14, nbInfo$adj)
## add to st vincent de barbeyrargues
nbInfo$num[14] <- 1
nbInfo$adj <- c(nbInfo$adj[1:sum(nbInfo$num[1:13])], 1, 
                nbInfo$adj[(length(nbInfo$adj) - 
                              nbInfo$num[15] + 1):length(nbInfo$adj)])
## adjust weights
nbInfo$weights <- rep(1, length(nbInfo$adj))

# get data and constants ready
K <- nrow(dat) # nsites
J <- ncol(dat[, -c(1, 2)]) # nsurveys
constants <- list(K = K, J = J,
                  L = length(nbInfo$adj), 
                  adj = nbInfo$adj, weights = nbInfo$weights, num = nbInfo$num) 
data <- list(
  y = dat[, -c(1, 2)],
  n = apply(dat[, -c(1, 2)], 1, sum),
  temp = (temperature$avg - mean(temperature$avg)) / 
               sd(temperature$avg)
  )

# model code
# code
model_code <- nimbleCode({
  for (i in 1:K) { # loop over sites
    
    # prob of cells
    pi[i,1] <- p[i,1] # survey 1
    for(j in 2:J){ # loop over surveys > 1
      pi[i,j] <- prod(1 - p[i, 1:(j-1)]) * p[i,j]
    }
    pcap[i] <- sum(pi[i,1:J])
    
    # capture prob
    for(j in 1:J){ # loop over surveys
      cloglog(p[i,j]) <- alpha[j]
      pic[i,j] <- pi[i,j] / pcap[i]
    }
    
    # likelihood: observation
    y[i,1:J] ~ dmulti(pic[i,1:J], n[i]) # multinomial for each site
    n[i] ~ dbin(pcap[i], N[i]) # for each site
    # likelihood: process
    N[i] ~ dpois(lambda[i]) # for each site
    log(lambda[i]) <- beta[1] + beta[2] * temp[i] + s[i]

  } # end of loop over sites
  
  # spatial latent process
  s[1:K] ~ dcar_normal(adj[1:L], weights[1:L], num[1:K], tau, zero_mean = 0)
  sigma ~ dunif(0, 100)   # prior for variance components based on Gelman (2006)
  tau <- 1 / sigma^2
  
  for (j in 1:J){
    alpha[j] ~ dnorm(0, sd = 1.5)
  }
  for (k in 1:2){
    beta[k] ~ dnorm(0, sd = 1.5)
  }
  
  #---------------- check model goodness of fit
  
  # for(i in 1:K){ 
  #   ## replicate data set
  #   n.pred[i] ~ dbin(pcap[i], N[i])
  #   y.pred[i,1:J] ~ dmulti(pic[i,1:J], n.pred[i])
  #   
  #   ## Freeman–Tukey residuals: observation component
  #   for(k in 1:J){ 
  #     e1[i,k]      <- pic[i,k] * n[i] # expected counts (obs)
  #     resid1[i,k]  <- pow(sqrt(y[i,k]) - sqrt(e1[i,k]), 2)
  #     e1.pred[i,k] <- pic[i,k] * n.pred[i] # expected counts (rep)
  #     resid1.pred[i,k] <- pow(sqrt(y.pred[i,k]) - sqrt(e1.pred[i,k]), 2)
  #   }  
  #   ## Freeman–Tukey residuals: abundance component
  #   e2[i] <- pcap[i] * lambda[i] 
  #   resid2[i] <- pow(sqrt(n[i]) - sqrt(e2[i]), 2) 
  #   resid2.pred[i] <- pow(sqrt(n.pred[i]) - sqrt(e2[i]), 2) 
  # }  
  # 
  # # fit statistic
  # fit1.data <- sum(resid1[1:K, 1:J])      # observation part (data)
  # fit1.pred <- sum(resid1.pred[1:K, 1:J]) # observation part (replicates)
  # fit2.data <- sum(resid2[1:K])           # abundance part (data)
  # fit2.pred <- sum(resid2.pred[1:K])      # abundance part (replicates)
  
})

# get initial values
pic.init <- array(runif(K * J, 0.05, 0.15), c(K, J))
pic.init <- pic.init / apply(pic.init, 1, sum)
#apply(pic.init, 1, sum)
Nin <- apply(dat[, -c(1, 2)], 1, sum) + 10
inits <- function(){
  list(alpha = rep(0.1, J), 
       pic = pic.init, 
       beta = rep(0, 2), 
       N = Nin,
       sigma = 1, 
       s = rnorm(K))
  }

# create Nimble model
Rmodel <- nimbleModel(code = model_code, 
                      constants = constants, 
                      data = data, 
                      inits = inits())
Rmcmc <- compileNimble(Rmodel, showCompilerOutput = F)
ModSpec <- configureMCMC(Rmodel, onlySlice = TRUE) # slice sampling 
ModSpec$resetMonitors()
ModSpec$addMonitors(c("alpha", "N", "beta", "sigma", "s"#,
                      #"fit1.data","fit1.pred","fit2.data","fit2.pred"
                      )
                    )
Cmcmc <- buildMCMC(ModSpec)
Cmodel <- compileNimble(Cmcmc, project = Rmodel, resetFunctions = TRUE)

# run model
ni <- 50000 * 4
nb <- 5000 * 4
nc <- 2
nt <- 15
samp_spat <- runMCMC(Cmodel, niter = ni, nburnin = nb, 
                     nchains = nc, thin = nt, inits = inits,  
                     samplesAsCodaMCMC = TRUE)

# get summary of results
MCMCsummary(samp_spat)

# save samples
saveRDS(samp_spat, "samples/samples_spatial_2022.rds")

# check convergence
MCMCtrace(samp_spat, params = c("alpha", "beta"), pdf = FALSE)


