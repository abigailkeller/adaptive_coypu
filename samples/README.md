`samples_spatio_timeranef_rw_twoalpha_longer.rds`

- like `samples_spatio_timeranef_hmc_twoalpha_long.rds` but just sampled even longer
- ni = 2000000 * 25
- nb = 50000 * 6
- nc = 5
- nt = 5000

`samples_spatio_timeranef_rw_twoalpha_long.rds`

- like `samples_spatio_timeranef_hmc_twoalpha.rds` but just sampled longer
- ni = 2000000 * 12
- nb = 50000 * 6
- nc = 4
- nt = 5000

`samples_spatio_timeranef_rw_twoalpha.rds`

- like `samples_spatio_timeranef_hmc_twoalpha_temp.rds` but just sampled longer with RW
- ni = 2000000 * 4
- nb = 50000 * 4
- nc = 4
- nt = 1000

`samples_spatio_timeranef_rw&amp;hmc_twoalpha.rds`

- HMC for ranefs, slice sampler for N, and RW for everything else

`samples_spatio_timeranef_hmc_twoalpha_contN.rds`

- made N a continuous variable, used HMC for everything

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