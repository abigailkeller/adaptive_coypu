library(tidyverse)
library(MCMCvis)
library(viridis)

# read in data
communes_shp <- st_read("data/shp/communes.shp")
temperature <- readRDS("data/temperature.rds")

# reorder communes
communes_shp <- communes_shp[match(temperature$commune, communes_shp$dpts_cl), ]


#####################
# plot all communes #
#####################

# plot communes
commune_names <- ggplot(data = communes_shp) +
  geom_sf(aes(fill = dpts_cl)) +
  labs(fill = "Commune") +
  theme_minimal()
ggsave("figures/commune_names.png", commune_names, height = 8, width = 8)

# plot 2022 N
samples_2022 <- MCMCsummary(readRDS("samples/samples_spatial_2022.rds"))
communes_shp$mean_2022 <- samples_2022[1:15, "mean"]
plot_2022 <- ggplot(data = communes_shp) +
  geom_sf(aes(fill = mean_2022)) +
  scale_fill_viridis() +
  ggtitle("2022 posterior") +
  labs(fill = "Mean N") +
  theme_minimal()
ggsave("figures/spatial_2022.png", plot_2022, height = 8, width = 8)

# plot multiyear N - full year
samples_my <- MCMCsummary(readRDS("samples/samples_multiyear.rds"))
df_my <- data.frame(
  N = samples_my[1:(10 * 15), "mean"],
  year = rep(2015:2024, times = 15)
) %>% 
  arrange(year) %>% 
  mutate(geometry = rep(communes_shp$geometry, 10))
plot_multiyear <- ggplot(data = st_as_sf(df_my)) +
  geom_sf(aes(fill = N)) +
  facet_wrap(~ year, ncol = 3) +
  scale_fill_viridis() +
  ggtitle("Months 1-12") +
  labs(fill = "Mean N") +
  theme_minimal()
ggsave("figures/multiyear_full.png", plot_multiyear, height = 8, width = 8)


# plot multiyear N - five months
samples_5m <- MCMCsummary(readRDS("samples/samples_multiyear_5mon.rds"))
df_5m <- data.frame(
  N = samples_5m[1:(10 * 15), "mean"],
  year = rep(2015:2024, times = 15)
) %>% 
  arrange(year) %>% 
  mutate(geometry = rep(communes_shp$geometry, 10))
plot_5m <- ggplot(data = st_as_sf(df_5m)) +
  geom_sf(aes(fill = N)) +
  facet_wrap(~ year, ncol = 3) +
  scale_fill_viridis() +
  ggtitle("Months 2-6") +
  labs(fill = "Mean N") +
  theme_minimal()
ggsave("figures/multiyear_5m.png", plot_5m, height = 8, width = 8)
