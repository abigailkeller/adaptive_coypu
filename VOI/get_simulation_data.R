##########################
# neighborhood structure #
##########################

# import data
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

# get list of neighbors
neighbors <- list()
start <- 1
for (i in 1:dim(commune_convert)[1]) {
  
  end <- start + nbInfo$num[i] - 1
  neighbors[[i]] <- nbInfo$adj[start:end]
  start <- start + nbInfo$num[i]
}

saveRDS(neighbors, "data/simulation/neighbors.rds")

########################
# get temperature data #
########################

temp_data <- list()

# get mean and sd landscape temperature
temp_data$landscape_mean <- mean(colMeans(temperature[, 2:10]))
temp_data$landscape_sd <- sd(colMeans(temperature[, 2:10]))

# get site-level z-score
z_scores <- apply(temperature[, 2:10], 2, function(x) (x - mean(x)) / sd(x))
temp_data$mean_z_score <- rowMeans(z_scores)
temp_data$sd_z_score <- apply(z_scores, 1, sd)
temp_data$sd <- mean(apply(temperature[, 2:10], 1, sd))

saveRDS(temp_data, "data/simulation/temp_data.rds")
