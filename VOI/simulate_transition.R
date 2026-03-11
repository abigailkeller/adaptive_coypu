library(MASS)
library(sp)
library(spdep)
library(tidyverse)

# read samples
samp <- readRDS("samples/samples_spatio_timeranef_rw_twoalpha_noint.rds")

get_transition <- function(params, N, a_loc, a_adj, 
                           nsites, adj, num, weights) {
  
  
  
}

#################################################### 
# get spatial random effect neighborhood structure #
####################################################

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

N <- 15
adj <- nbInfo$adj
num <- nbInfo$num
weights <- nbInfo$weights

# Create an empty Adjacency Matrix (W)
W <- matrix(0, nrow = N, ncol = N)
pos <- 1
for(i in 1:N) {
  W[i, adj[pos:(pos + num[i] - 1)]] <- weights[pos:(pos + num[i] - 1)]
  pos <- pos + num[i]
}

# Create the Diagonal Matrix (D)
D <- diag(num)

simulate_car <- function(W, D, tau = 1, alpha = 0.99) {
  # Precision Matrix: Q = tau * (D - alpha * W)
  # alpha < 1 ensures the matrix is invertible (Proper CAR)
  Q <- tau * (D - alpha * W)
  
  # Covariance Matrix is the inverse of Precision
  Sigma <- solve(Q)
  
  # Simulate from Multivariate Normal
  return(mvrnorm(n = 1, mu = rep(0, nrow(W)), Sigma = Sigma))
}

tau_posterior <- 1 / c(samp[[1]][, "sigma_s"], samp[[2]][, "sigma_s"],
                       samp[[3]][, "sigma_s"], samp[[4]][, "sigma_s"]) ^ 2

# Run simulation
spatial_effects <- data.frame(
  min = simulate_car(W, D, tau = min(tau_posterior), alpha = 0.99),
  max = simulate_car(W, D, tau = max(tau_posterior), alpha = 0.99)
)

# add spatial effects to commune data
communes_shp <- cbind(communes_shp, spatial_effects)

ggplot(data = communes_shp) +
  geom_sf(aes(fill = max)) +
  labs(fill = "Commune") +
  ggtitle(paste0("tau = ", round(max(tau_posterior), 3))) +
  theme_minimal()

ggplot(data = communes_shp) +
  geom_sf(aes(fill = min)) +
  labs(fill = "Commune") +
  ggtitle(paste0("tau = ", round(min(tau_posterior), 3))) +
  theme_minimal()


