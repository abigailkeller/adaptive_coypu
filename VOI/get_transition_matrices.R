library(MASS)
library(sp)
library(spdep)
library(tidyverse)
library(patchwork)

# transition function
get_transition <- function(params, N, dd, N_seq,
                           n_sites, n_visits, p_increase, site_increase) {
  
  ################################
  # get probability of detection #
  ################################
  
  p <- matrix(params["p_detect"], nrow = n_sites, ncol = n_visits)
  p[site_increase, ] <- params["p_detect"] * p_increase
  
  # get annual probability of capture
  pi <- matrix(NA, nrow = n_sites, ncol = n_visits)
  pcap <- rep(NA, n_sites)
  n <- rep(NA, n_sites)
  
  pi[, 1] <- p[, 1] # survey 1
  for (i in 1:n_sites) { # loop over sites
    
    # prob of cells
      for (j in 2:n_visits) { # loop over surveys > 1
        pi[i, j] <- prod(1 - p[i, 1:(j - 1)]) * p[i, j] 
      }
    pcap[i] <- sum(pi[i, ])
    
    prob_S <- dbinom(N_seq, size = round(N[i]), prob = 1 - pcap[i])
    
    growth_term <- params["rho"] * N_seq * (1 - N_seq / exp(dd[i])) + 1
    
    
    
    
    # removal
    n[i] <- rbinom(1, round(N[i]), pcap[i])
  } 
  
  
  # remaining after removal
  remaining <- round(N) - n
  
  # draw s_t
  s_t <- rnorm(1, params["st_mean"], sd = params["st_sd"])
  
  # advance to next year
  lambda_next <- (
    s_t + log(params["rho"] * remaining * (1 - remaining / exp(dd)) + 1)
  )
  
  # get nextN
  nextN <- rep(NA, n_sites)
  for (i in 1:n_sites) {
    nextN[i] <- rpois(1, exp(lambda_next[i]))
  }
  
  return(nextN)
  
}

# transition function for one site
seq_st <- seq(-2, 2, 0.1)


get_transition_onesite <- function(params, N, dd, N_seq,
                                   n_visits, p_increase, a) {
  
  ################################
  # get probability of detection #
  ################################
  
  p <- pi <- if (a) {
    as.vector(rep(params["p_detect"] * p_increase, n_visits))
  } else {
    as.vector(rep(params["p_detect"], n_visits))
  }
  pi[2:length(pi)] <- cumprod(1 - p[-length(pi)]) * p[-1]   
  
  pcap <- sum(pi)
  
  prob_S <- dbinom(N_seq, size = round(N), prob = 1 - pcap)
  
  growth_term <- params["rho"] * N_seq * (1 - N_seq / exp(dd)) + 1
  
  # draw s_t
  s_t <- dnorm(seq_st, params["st_mean"], sd = params["st_sd"])
  
  # advance to next year
  for (i in seq_along(N_seq)) {
    lambda_next <- (
      log(exp(seq_st) + params["rho"] * N_seq[i] * (1 - N_seq[i] / exp(dd)))
    )
  }
  
  for(c in seq_along(Sa)){
    for(d in seq_along(f)){
      for(e in seq_along(Sl)){
        for(g in seq_along(alpha_seq)){
          for(h in seq_along(removal_seq_sub)){
            
            out <- growth(Sa[c],f[d],Sl[e],ratio,states[a],K,
                          alpha_seq[g],removal_seq_sub[h])
            index <- which.min(abs(states-out))
            sprime[index] <- (
              sprime[index] +
                Sa_prob[c]*f_prob[d]*Sl_prob[e]*alpha_prob[g]*
                removal_prob[h]
            )
            
          }
        }
      }
    }
  }
  # print(Sys.time()-time1)
  # add to transition matrix
  mat[,b] <- sprime/sum(sprime)
  
  # get nextN
  nextN <- rep(NA, n_sites)
  for (i in 1:n_sites) {
    nextN[i] <- rpois(1, exp(lambda_next[i]))
  }
  
  return(nextN)
  
}

get_transition_simple <- function(params, N, dd, N_seq,
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
  for (j in seq_alon(N_seq)) {
    
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
  
  return(nextN)
  
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


##############################
# create transition matrices #
##############################

# state space
states <- seq(0, 300, 10)
actions <- combn(1:15, 2)

transition_min <- array(NA, dim = c(length(states) * n_sites, 
                                    length(states) * n_sites, 
                                    length(actions)))
get_transition(params, N, dd,
                           n_sites, n_visits, p_increase, site_increase)

