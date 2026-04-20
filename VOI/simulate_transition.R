library(MASS)
library(sp)
library(spdep)
library(tidyverse)

# transition function
get_transition <- function(params, N, neighbors, temp_data,
                           n_sites, n_visits, p_increase, site_increase) {
  
  ##########################
  # get annual temperature #
  ##########################
  
  # get landscape mean
  mean_temp <- rnorm(1, temp_data$landscape_mean, temp_data$landscape_sd)
  
  # get site-level z_score and temp
  z_score <- rep(NA, n_sites)
  for (i in 1:n_sites) {
    z_score[i] <- rnorm(1, temp_data$mean_z_score[i], temp_data$sd_z_score[i])
  }
  temp <- mean_temp + z_score * temp_data$sd
  
  
  ################################
  # get probability of detection #
  ################################
  
  p <- matrix(params["p_detect"], nrow = n_sites, ncol = n_visits)
  p[site_increase, ] <- params["p_detect"] * p_increase
  
  # get annual probability of capture
  pi <- matrix(NA, nrow = n_sites, ncol = n_visits)
  pcap <- rep(NA, n_sites)
  n <- rep(NA, n_sites)
  
  pi[, 1] <- p[, 1] # survey 1
  for (i in 1:n_sites) { # loop over sites
    
    # prob of cells
      for (j in 2:n_visits) { # loop over surveys > 1
        pi[i, j] <- prod(1 - p[i, 1:(j - 1)]) * p[i, j] 
      }
    pcap[i] <- sum(pi[i, ])
    
    # removal
    n[i] <- rbinom(1, round(N[i]), pcap[i])
  } 
  
  
  # remaining after removal
  remaining <- round(N) - n
  
  # get the mean of neighbors
  neighbors_rem <- rep(NA, n_sites)
  for (i in 1:n_sites) {
    neighbors_rem[i] <- mean(remaining[neighbors[[i]]])
  }
  
  # draw s_t
  s_t <- rnorm(1, params["st_mean"], sd = params["st_sd"])
  
  # advance to next year
  lambda_next <- (
    params["beta"] * temp + s_t + 
      params["rho"] * log(remaining + 1) +
      params["gamma"] * log(neighbors_rem + 1)
  )
  
  # get nextN
  nextN <- rep(NA, n_sites)
  for (i in 1:n_sites) {
    nextN[i] <- rpois(1, exp(lambda_next[i]))
  }
  
  return(nextN)
  
}

# read samples
samp <- do.call(rbind,
                readRDS("samples/samples_spillover_timeranef_rw_twoalpha_unif.rds"))

# read in simulation data
neighbors <- readRDS("data/simulation/neighbors.rds")
temp_data <- readRDS("data/simulation/temp_data.rds")

# get Nstart
Nstart <- colMeans(samp[, 1:15 * 9])

# get indices
index_min <- which.min(samp[, "gamma"])
index_max <- which.max(samp[, "gamma"])

# simulate 10 years of min index
N_sim_min <- matrix(NA, nrow = 11, ncol = 15)
N_sim_min[1, ] <- Nstart
for (t in 1:10) {
  N_sim_min[t + 1, ] <- get_transition(params = samp[index_min,], 
                                       N = N_sim_min[t, ], 
                                       neighbors = neighbors, 
                                       temp_data = temp_data, 
                                       n_sites = 15, n_visits = 8,
                                       p_increase = 2, site_increase = c(4, 13))
}

# simulate 10 years of max index
N_sim_max <- matrix(NA, nrow = 11, ncol = 15)
N_sim_max[1, ] <- Nstart
for (t in 1:10) {
  N_sim_max[t + 1, ] <- get_transition(params = samp[index_max,], 
                                       N = N_sim_max[t, ], 
                                       neighbors = neighbors, 
                                       temp_data = temp_data, 
                                       n_sites = 15, n_visits = 8,
                                       p_increase = 2, site_increase = c(4, 13))
}


########
# plot #
########

# read in commune shapefile
temperature <- readRDS("data/temp/temperature.rds")
communes_shp <- st_read("data/communes_shp/communes.shp")

# reorder communes
communes_shp <- communes_shp[match(temperature$commune, communes_shp$dpts_cl), ]

# add Nstart to commune data
communes_shp <- cbind(communes_shp, Nstart)

ggplot(data = communes_shp) +
  geom_sf(aes(fill = Nstart)) +
  labs(fill = "N") +
  theme_minimal()

ggplot(data = communes_shp) +
  geom_sf(aes(fill = dpts_cl)) +
  labs(fill = "Commune") +
  theme_minimal()

ggplot(data = communes_shp) +
  geom_sf(aes(fill = max)) +
  labs(fill = "Commune") +
  ggtitle(paste0("tau = ", round(max(tau_posterior), 3))) +
  theme_minimal()

ggplot(data = samp) +
  geom_point(aes(x = gamma, y = rho)) +
  theme_minimal()


