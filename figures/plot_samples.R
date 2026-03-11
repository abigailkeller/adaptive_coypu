library(tidyverse)
library(MCMCvis)
library(viridis)

samp <- readRDS("samples/samples_spatio_timeranef_rw_twoalpha_noint.rds")

param1 <- "N[7, 2]"
param2 <- "N[7, 3]"
param3 <- "sigma_s"

ggplot() +
  geom_point(aes(x = c(samp[[1]][, param1], samp[[2]][, param1],
                       samp[[3]][, param1], samp[[4]][, param1]),
                 y = c(samp[[1]][, param2], samp[[2]][, param2], 
                       samp[[3]][, param2], samp[[4]][, param2]))) +
  labs(x = param1, y = param2)

ggplot() +
  geom_point(aes(x = c(samp[[1]][, param1], samp[[2]][, param1],
                       samp[[3]][, param1], samp[[4]][, param1]),
                 y = c(samp[[1]][, param2], samp[[2]][, param2], 
                       samp[[3]][, param2], samp[[4]][, param2]),
                 color = c(samp[[1]][, param3], samp[[2]][, param3], 
                           samp[[3]][, param3], samp[[4]][, param3]))) +
  labs(x = param1, y = param2, color = param3) +
  scale_color_viridis()

summary <- MCMCsummary(samp)

param <- "sigma_s"

ggplot() +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[1]][, param]), color = "blue") +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[2]][, param]), color = "purple") +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[3]][, param]), color = "yellow") +
  geom_line(aes(x = 1:length(samp[[1]][, param]),
                y = samp[[4]][, param]), color = "green") +
  labs(x = "iteration", y = param)

ggplot() +
  geom_line(aes(x = 1:length(samp2[[1]][, param]),
                y = samp2[[1]][, param]), color = "blue") +
  geom_line(aes(x = 1:length(samp2[[1]][, param]),
                y = samp2[[2]][, param]), color = "purple") +
  geom_line(aes(x = 1:length(samp2[[1]][, param]),
                y = samp2[[3]][, param]), color = "yellow") +
  geom_line(aes(x = 1:length(samp2[[1]][, param]),
                y = samp2[[4]][, param]), color = "green") +
  labs(x = "iteration", y = param)

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

Ntime <- data.frame(
  N = summary[1:135, "mean"],
  t = rep(1:9, 15),
  site = rep(1:15, each = 9)
)

ggplot() +
  geom_line(data = Ntime,
            aes(x = t, y = N, color = as.factor(site))) +
  labs(x = "t", y = "N", color = "site")

