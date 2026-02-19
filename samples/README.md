`samples_spatio_timeranef_hmc_twoalpha.rds`

- p[t, i, j] <- p_detect * sampled[t, i, j]
- p_detect ~ dunif(0, 1)
- p_sample ~ dunif(0, 1)

`samples_spatio_timeranef_hmc_twoalpha_transform.rds`

- p[t, i, j] <- ilogit(p_detect) * sampled[t, i, j]
- p_detect ~ dunif(-10, 10)
- p_sample ~ dunif(0, 1)

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