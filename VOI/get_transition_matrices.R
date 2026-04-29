library(MASS)
library(sp)
library(spdep)
library(tidyverse)
library(patchwork)
library(MCMCvis)

get_transition_simple <- function(params, dd, N_seq, seq_st,
                                  n_visits, p_increase, a) {
  
  out <- matrix(NA, nrow = length(N_seq), ncol = length(N_seq))
  
  # get probability of detection
  p <- pi <- if (a) {
    as.vector(rep(params["p_detect"] * p_increase, n_visits))
  } else {
    as.vector(rep(params["p_detect"], n_visits))
  }
  pi[2:length(pi)] <- cumprod(1 - p[-length(pi)]) * p[-1]   
  
  pcap <- sum(pi)
  
  # draw s_t
  s_t <- dnorm(seq_st, params["st_mean"], sd = params["st_sd"])
  
  # loop through N
  for (j in seq_along(N_seq)) {
    
    survived <- N_seq[j] * (1 - pcap)
    
    lambda_next <- (
      exp(seq_st) + params["rho"] * survived * (1 - survived / exp(dd))
    )
    
    prob <- rep(0, length(N_seq))
    
    for (i in seq_along(lambda_next)) {
      
      index <- which.min(abs(lambda_next[i] - N_seq))
      prob[index] <- prob[index] + s_t[i]
      
    }
    
    out[, j] <- prob / sum(prob)
    
  }
  
  return(out)
  
}

# read samples
sample_file <- "samples/samples_timeranef_rw_twoalpha_dd_nobeta.rds"
samp <- do.call(rbind, readRDS(sample_file))

# get summary of samples
summary <- MCMCsummary(readRDS(sample_file))

# get Nstart indices
Nstart_ind <- 1:15 * 9

# get dd indices
dd_ind <- 136:150

# get indices
index_min <- which.min(samp[, "rho"])
index_max <- which.max(samp[, "rho"])


##############################
# create transition matrices #
##############################

# state and action spaces

# create states based on carrying capacity of each site
n_states <- 8
n_sites <- 15
n_choice <- 2

states <- matrix(NA, nrow = n_states, ncol = n_sites)
for (i in 1:n_sites) {
  
  states[, i] <- seq(from = 0, to = exp(summary[dd_ind[i], "mean"]), 
                     length.out = n_states)
  
}

action_set <- combn(n_sites, n_choice)
actions <- vapply(seq_len(ncol(action_set)), 
                  function(i) tabulate(action_set[, i], nbins = n_sites), 
                  integer(n_sites))
# do nothing action
actions <- cbind(actions, rep(0, n_sites))


# sequence of possible s_t
seq_st <- seq(-2, 2, 0.1)


transition_min <- array(NA, dim = c(n_states, 
                                    n_states * n_sites, 
                                    length(actions)))

for (i in seq_along(1:n_sites)) {
  
  col_start <- (i - 1) * n_states + 1
  col_end <- i * n_states
  
  for (j in seq_along(1:ncol(actions))) {
    
    transition_min[, col_start:col_end, j] <- get_transition_simple(
      params = samp[index_min, ], dd = samp[index_min, dd_ind[i]], 
      N_seq = states[, i], seq_st, n_visits = 8, p_increase = 2, 
      a = actions[i, j]
      )
  }
}


transition_max <- array(NA, dim = c(n_states, 
                                    n_states * n_sites, 
                                    length(actions)))

for (i in seq_along(1:n_sites)) {
  
  col_start <- (i - 1) * n_states + 1
  col_end <- i * n_states
  
  for (j in seq_along(1:ncol(actions))) {
    
    transition_max[, col_start:col_end, j] <- get_transition_simple(
      params = samp[index_max, ], dd = samp[index_max, dd_ind[i]], 
      N_seq = states[, i], seq_st, n_visits = 8, p_increase = 2, 
      a = actions[i, j]
    )
  }
}

# save transition matrices
saveRDS(transition_min, "transition_matrices/t_min.rds")
saveRDS(transition_max, "transition_matrices/t_max.rds")
