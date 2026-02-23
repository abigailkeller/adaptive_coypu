`samples_spatio_timeranef_hmc_twoalpha_temp.rds`

- log(lambda[t, i]) <- beta[1] + beta[2] * temp[t, i] + s_s[i] + s_t[t]
- p[t, i, j] <- p_detect * sampled[t, i, j]
- p_detect ~ dunif(0, 1)
- p_sample ~ dunif(0, 1)

`samples_spatio_timeranef_hmc_twoalpha_prior.rds`

- p_detect ~ dlnorm(meanlog = log(0.02), sdlog = 0.3)

`samples_spatio_timeranef_hmc_twoalpha.rds`

- p[t, i, j] <- p_detect * sampled[t, i, j]
- p_detect ~ dunif(0, 1)
- p_sample ~ dunif(0, 1)

`samples_spatio_timeranef_hmc_twoalpha_transform.rds`

- p[t, i, j] <- ilogit(p_detect) * sampled[t, i, j]
- p_detect ~ dunif(-10, 10)
- p_sample ~ dunif(0, 1)

`samples_spatio_timeranef_hmc_twoalpha_mu.rds`

- N[t, i] ~ dpois(mu[t, i] / p_detect) 
- log(mu[t, i]) <- beta[1] + beta[2] * temp[t, i] + s_s[i] + s_t[t]

`samples_spatio_timeranef_hmc_onealpha.rds`

- alpha ~ dunif(-10, 10)

`samples_spatio_timeranef_hmc_mean-1.rds`
 
 - alpha[j] ~ dnorm(-1, sd = 1)
 
 `samples_spatio_timeranef_hmc_mean0.rds`
 
 - alpha[j] ~ dnorm(0, sd = 1.5)
 
 `samples_spatio_timeranef_hmc_betadist.rds`

 - alpha[j] ~ dbeta(1, 1)
 
 `samples_spatio_timeranef_hmc_unifdist.rds`
 
  - alpha[j] ~ dunif(0, 1)