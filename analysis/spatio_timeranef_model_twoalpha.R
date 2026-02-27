library(sf)
library(tidyverse)
library(janitor)
library(lubridate)
library(nimble)
#library(ubms)
library(MCMCvis)
library(sp)
library(spdep)
library(parallel)
library(nimbleHMC)


# import data
raw_dat <- readRDS("data/coypus.rds")
temperature <- readRDS("data/temp/temperature.rds")
communes_shp <- st_read("data/communes_shp/communes.shp")

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

# convert temp data to matrix
temp_mat <- t(as.matrix(temperature[, 2:10]))

# plot communes
# ggplot(data = communes_shp) +
#   geom_sf(aes(fill = dpts_cl)) +
#   labs(fill = "Commune") +
#   theme_minimal()

# format data
dat <- raw_dat$removal_coypus %>%
  filter(Month %in% c(2, 3, 4, 5, 6),
         Year > 2015) %>%
  # mutate(Year = as_factor(Year)) %>%
  pivot_wider(names_from = Month, 
              values_from = n)
dat_st <- array(NA, dim = c(length(unique(dat$Year)),
                            length(unique(dat$commune)),
                            5))
year_vec <- unique(dat$Year)
for (i in seq_along(year_vec)) {
  dat_st[i, , ] <- as.matrix(dat[dat$Year == year_vec[i], 3:7])
}

# get total n
n <- matrix(NA, nrow = length(unique(dat$Year)), 
            ncol = length(unique(dat$commune)))
for (i in 1:length(year_vec)) {
  n[i, ] <- rowSums(dat_st[i, , ])
}

# add temp and raw data to spatial data
# communes_shp$temp <- temperature$avg
# communes_shp$remove_2022 <- rowSums(dat[dat$Year == 2022, 3:7])

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
K <- dim(dat_st)[2] # nsites
J <- dim(dat_st)[3] # nsurveys
M <- dim(dat_st)[1] # nyears
constants <- list(K = K, J = J, M = M,
                  L_s = length(nbInfo$adj), 
                  adj_s = nbInfo$adj, weights_s = nbInfo$weights, 
                  num_s = nbInfo$num) 
data <- list(
  y = dat_st,
  n = n,
  temp = (temp_mat - mean(temp_mat)) / 
               sd(temp_mat)
  )

# model code
# code
model_code <- nimbleCode({
  
  for (t in 1:M) { # loop over years
    
    # year ranef
    s_t[t] ~ dnorm(0, sd = 1.5)
    
    for (i in 1:K) { # loop over sites
      
      # likelihood: observation
      y[t, i, 1:J] ~ dmulti(pic[t, i, 1:J], n[t, i]) # multinomial for each site
      n[t, i] ~ dbin(pcap[t, i], N[t, i]) # for each site
      # likelihood: process
      N[t, i] ~ dpois(lambda[t, i]) # for each site
      log(lambda[t, i]) <- beta[1] + beta[2] * temp[t, i] + s_s[i] + s_t[t]
      
      # capture prob
      for (j in 1:J){
        sampled[t, i, j] ~ dbin(p_sample, 1)
        p[t, i, j] <- p_detect * sampled[t, i, j]
        pic[t, i, j] <- pi[t, i, j] / pcap[t, i]
      }
      
    } # end of loop over years
  } # end of loop over sites
  
  for (i in 1:K) { # loop over sites
    
    # prob of cells
    pi[1:M, i, 1] <- p[1:M, i, 1] # survey 1
    for (t in 1:M) { # loop over years
      for (j in 2:J) { # loop over surveys > 1
        pi[t, i, j] <- prod(1 - p[t, i, 1:(j - 1)]) * p[t, i, j] 
      }
      pcap[t, i] <- sum(pi[t, i, 1:J])
    }
    
  } # end of loop over sites
  
  # spatial latent process
  s_s[1:K] ~ dcar_normal(adj_s[1:L_s], weights_s[1:L_s], num_s[1:K], 
                         tau_s, zero_mean = 0)
  sigma_s ~ dunif(0, 100)   # prior for variance components based on Gelman (2006)
  tau_s <- 1 / sigma_s ^ 2
  
  p_detect ~ dunif(0, 1)
  p_sample ~ dunif(0, 1)
  
  for (k in 1:2){
    beta[k] ~ dnorm(0, sd = 1.5)
  }
  
  #---------------- check model goodness of fit
  
  # for(i in 1:K){ 
  #   ## replicate data set
  #   n.pred[i] ~ dbin(pcap[i], N[i])
  #   y.pred[i,1:J] ~ dmulti(pic[i,1:J], n.pred[i])
  #   
  #   ## Freeman???Tukey residuals: observation component
  #   for(k in 1:J){ 
  #     e1[i,k]      <- pic[i,k] * n[i] # expected counts (obs)
  #     resid1[i,k]  <- pow(sqrt(y[i,k]) - sqrt(e1[i,k]), 2)
  #     e1.pred[i,k] <- pic[i,k] * n.pred[i] # expected counts (rep)
  #     resid1.pred[i,k] <- pow(sqrt(y.pred[i,k]) - sqrt(e1.pred[i,k]), 2)
  #   }  
  #   ## Freeman???Tukey residuals: abundance component
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
pic.init <- array(NA, dim = c(M, K, J))
for (i in 1:M) {
  temp <- array(runif(K * J, 0.05, 0.15), c(K, J))
  pic.init[i, , ] <- temp / apply(temp, 1, sum)
}
#apply(pic.init, 1, sum)
Nin <- n + 10
inits <- function(){
  list(p_detect = runif(1, 0, 1), 
       p_sample = runif(1, 0, 1), 
       sampled = array(1, dim(dat_st)),
       pic = pic.init, 
       beta = rep(0, 2), 
       N = Nin,
       sigma_s = 1, 
       s_s = rnorm(K),
       s_t = rnorm(M))
  }

# run model
ni <- 2000000 * 25
nb <- 50000 * 6
nc <- 6
nt <- 10000
# ni <- 5000 * 4
# nb <- 1000 * 4
# nc <- 4
# nt <- 1

########################
# run MCMC in parallel #
########################

cl <- makeCluster(nc)

clusterExport(cl, c("model_code", "inits", "data", "constants",
                    "ni", "nt", "J", "K", "pic.init", "Nin", "M", "dat_st"))

out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  library(nimbleHMC)
  
  # create Nimble model
  Rmodel <- nimbleModel(code = model_code, 
                        constants = constants, 
                        data = data, 
                        inits = inits(),
                        buildDerivs = TRUE)
  
  # compile the model
  Rmcmc <- compileNimble(Rmodel, showCompilerOutput = F)
  
  # build the MCMC
  ModSpec <- configureMCMC(Rmodel, 
                           monitors = c("p_sample", "p_detect", "sampled",
                                        "N", "beta", "sigma_s",
                                        "s_s", "s_t")
  )
  
  # add HMC
  # ModSpec$removeSamplers(target = c("s_s", "s_t", "beta", "sigma_s"))
  # ModSpec$addSampler(target = c("s_s", "s_t", "beta", "sigma_s"),
  #                    type = "NUTS")
  
  Cmcmc <- buildMCMC(ModSpec)
  
  # compile MCMC and model
  Cmodel <- compileNimble(Cmcmc, project = Rmodel, resetFunctions = TRUE)
  
  # run MCMC
  Cmodel$run(niter = ni, thin = nt)
  
  samples <- as.mcmc(as.matrix(Cmodel$mvSamples))
  
  return(samples)
  
})

# discard burnin
sequence <- seq(nb / nt, ni / nt, 1)
out_sub <- list(out[[1]][sequence, ], out[[2]][sequence, ],
                out[[3]][sequence, ], out[[4]][sequence, ])

# save samples
saveRDS(out_sub, "samples/samples_spatio_timeranef_rw_twoalpha_longer.rds")

stopCluster(cl)



