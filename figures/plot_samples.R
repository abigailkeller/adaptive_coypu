library(tidyverse)
library(MCMCvis)

samp <- readRDS("samples/samples_spatio_timeranef_hmc_test.rds")

param1 <- "s_s[5]"
param2 <- "s_s[1]"

ggplot() +
  geom_point(aes(x = samp[[1]][, param1],
                 y = samp[[1]][, param2]))

summary <- MCMCsummary(samp)

param <- "N[4, 13]"

ggplot() +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[1]][, param]), color = "blue") +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[2]][, param]), color = "purple") +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[3]][, param]), color = "yellow") +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[4]][, param]), color = "green")

alpha <- out[151:1950, "mean"]
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
