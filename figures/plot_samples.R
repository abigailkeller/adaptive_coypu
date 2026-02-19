library(tidyverse)
library(MCMCvis)

samp <- readRDS("samples/samples_spatio_timeranef_hmc_twoalpha.rds")

param1 <- "N[2, 4]"
param2 <- "p_detect"

ggplot() +
  geom_point(aes(x = samp[[2]][, param1],
                 y = samp[[2]][, param2])) +
  labs(x = param1, y = param2)

summary <- MCMCsummary(samp)

param <- "s_s[2]"

ggplot() +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[1]][, param]), color = "blue") +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[2]][, param]), color = "purple") +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[3]][, param]), color = "yellow") +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[4]][, param]), color = "green")

alpha <- summary[151:155, "mean"]
hist(alpha)
hist(1 - exp(-exp(alpha)))

 ggplot() +
  geom_density(aes(x = alpha), fill = "blue", alpha = 0.6) +
  geom_density(aes(x = rnorm(10000, 0, 1.5)), alpha = 0.6)

ggplot() +
  geom_density(aes(x = 1 - exp(-exp(alpha))), fill = "blue", alpha = 0.6) +
  geom_density(aes(x = 1 - exp(-exp(rnorm(10000, 0, 1.5)))), alpha = 0.6)

ggplot() +
  geom_density(aes(x = 1 - exp(-exp(alpha))), fill = "blue", alpha = 0.6) +
  geom_density(aes(x = 1 - exp(-exp(runif(10000, -5, 0)))), alpha = 0.6)
