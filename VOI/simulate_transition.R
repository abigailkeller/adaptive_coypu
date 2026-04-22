library(MASS)
library(sp)
library(spdep)
library(tidyverse)
library(patchwork)

# transition function
get_transition <- function(params, N, dd,
                           n_sites, n_visits, p_increase, site_increase) {
  
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
  
  # draw s_t
  s_t <- rnorm(1, params["st_mean"], sd = params["st_sd"])
  
  # advance to next year
  lambda_next <- (
    s_t + log(params["rho"] * remaining * (1 - remaining / exp(dd)) + 1)
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
                readRDS("samples/samples_timeranef_rw_twoalpha_dd_nobeta.rds"))

# get Nstart indices
Nstart_ind <- 1:15 * 9

# get dd indices
dd_ind <- 136:150

# get indices
index_min <- which.min(samp[, "rho"])
index_max <- which.max(samp[, "rho"])


# simulate 10 years of min index
N_sim_min <- matrix(NA, nrow = 11, ncol = 15)
N_sim_min[1, ] <- samp[index_min, Nstart_ind]
for (t in 1:10) {
  N_sim_min[t + 1, ] <- get_transition(params = samp[index_min,], 
                                       N = N_sim_min[1, ], 
                                       dd = samp[index_min, dd_ind], 
                                       n_sites = 15, n_visits = 8,
                                       p_increase = 2, site_increase = c(4, 13))
}

# simulate 10 years of max index
N_sim_max <- matrix(NA, nrow = 11, ncol = 15)
N_sim_max[1, ] <- samp[index_max, Nstart_ind]
for (t in 1:10) {
  N_sim_max[t + 1, ] <- get_transition(params = samp[index_max,], 
                                       N = N_sim_max[1, ], 
                                       dd = samp[index_max, dd_ind], 
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
communes_shp$Nstart_min <- N_sim_min[1, ]
communes_shp$Nstart_max <- N_sim_max[1, ]

# add Nend to commune data
communes_shp$Nend_min <- N_sim_min[11, ]
communes_shp$Nend_max <- N_sim_max[11, ]

plot_Nstart_min <- ggplot(data = communes_shp) +
  geom_sf(aes(fill = Nstart_min)) +
  labs(fill = "N") +
  ggtitle("rho min, t = 1") +
  theme_minimal()
plot_Nend_min <- ggplot(data = communes_shp) +
  geom_sf(aes(fill = Nend_min)) +
  labs(fill = "N") +
  ggtitle("rho min, t = 11") +
  theme_minimal()

plot_Nstart_max <- ggplot(data = communes_shp) +
  geom_sf(aes(fill = Nstart_max)) +
  labs(fill = "N") +
  ggtitle("rho max, t = 1") +
  theme_minimal()
plot_Nend_max <- ggplot(data = communes_shp) +
  geom_sf(aes(fill = Nend_max)) +
  labs(fill = "N") +
  ggtitle("rho max, t = 11") +
  theme_minimal()

plot_Nstart_min + plot_Nend_min + plot_Nstart_max + plot_Nend_max +
  plot_layout(nrow = 2, guides = "collect") & 
  scale_fill_viridis_c(limits = c(0, 
                                  max(c(communes_shp$Nstart_min, 
                                        communes_shp$Nstart_max,
                                        communes_shp$Nend_min,
                                        communes_shp$Nend_max))))

##############################
# create transition matrices #
##############################

# state space
states <- seq(0, 300, 5)
actions <- combn(1:15, 2)

transition_min <- array(NA, nrow = )
get_transition(params, N, dd,
                           n_sites, n_visits, p_increase, site_increase)

