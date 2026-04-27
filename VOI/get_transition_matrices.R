library(MASS)
library(sp)
library(spdep)
library(tidyverse)
library(patchwork)

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
samp <- do.call(rbind,
                readRDS("samples/samples_timeranef_rw_twoalpha_dd_nobeta.rds"))

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
N_seq <- seq(0, 300, 20)
states <- expand.grid(N_seq, N_seq, N_seq)
# actions <- combn(1:15, 2)
actions <- combn(1:3, 2)

n_sites <- 15

seq_st <- seq(-2, 2, 0.1)

# Initialize matrix
transition_min <- array(0, dim = c(nrow(states), nrow(states), ncol(actions)))

for (s in 1:n_states) {
  for (s_next in 1:n_states) {
    # Current components
    curr_x <- states$x[s]
    curr_y <- states$y[s]
    
    # Next components
    next_x <- states$x[s_next]
    next_y <- states$y[s_next]
    
    # Probabilities
    # Prob_y only looks at y, ignores x
    prob_x <- get_prob_x(next_x, curr_x, curr_y, action)
    prob_y <- get_prob_y(next_y, curr_y, action) 
    
    T_matrix[s, s_next] <- prob_x * prob_y
  }
}





transition_min <- array(NA, dim = c(length(states) * n_sites, 
                                    length(states) * n_sites, 
                                    length(actions)))

for (i in seq_along(1:n_sites)) {
  for (j in seq_along(1:ncol(actions))) {
    get_transition_simple(samp[index_min,], N, dd,
                          n_sites, n_visits, p_increase, site_increase)
  }
}
test <- get_transition_simple(params = samp[index_min,], 
                              dd = samp[index_min, dd_ind[1]], 
                              N_seq = states, seq_st,
                              n_visits = 8, p_increase = 2, a = 1)

(samp[index_min,], N, dd,
               n_sites, n_visits, p_increase, site_increase)

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

# Define state spaces
X <- 1:5
Y <- 1:3
states <- expand.grid(y = Y, x = X) # Y varies fastest to match matrix row-major
n_states <- nrow(states)

# Initialize matrix
T_matrix <- matrix(0, nrow = n_states, ncol = n_states)

for (s in 1:n_states) {
  for (s_next in 1:n_states) {
    # Current components
    curr_x <- states$x[s]
    curr_y <- states$y[s]
    
    # Next components
    next_x <- states$x[s_next]
    next_y <- states$y[s_next]
    
    # Probabilities
    # Prob_y only looks at y, ignores x
    prob_x <- get_prob_x(next_x, curr_x, curr_y, action)
    prob_y <- get_prob_y(next_y, curr_y, action) 
    
    T_matrix[s, s_next] <- prob_x * prob_y
  }
}
