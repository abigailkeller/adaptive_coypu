# read samples
sample_file <- "samples/samples_timeranef_rw_twoalpha_dd_nobeta.rds"
samp <- do.call(rbind, readRDS(sample_file))

# get summary of samples
summary <- MCMCsummary(readRDS(sample_file))

# get Nstart indices
Nstart_ind <- 1:15 * 9

# get dd indices
dd_ind <- 136:150

# create states based on carrying capacity of each site
n_states <- 4
n_sites <- 5
n_choice <- 2

# get states and actions
action_set <- combn(n_sites, n_choice)
actions <- vapply(seq_len(ncol(action_set)), 
                  function(i) tabulate(action_set[, i], nbins = n_sites), 
                  integer(n_sites))
# do nothing action
actions <- cbind(actions, rep(0, n_sites))

states <- matrix(NA, nrow = n_states, ncol = n_sites)
for (i in 1:n_sites) {
  
  states[, i] <- seq(from = 0, to = exp(summary[dd_ind[i], "mean"]), 
                     length.out = n_states)
  
}

n_states_all <- n_states ^ n_sites

state_indices <- expand.grid(rep(list(c(1, 2, 3, 4, 5)), 5))

utility <- matrix(NA, nrow = n_states_all, ncol = col(actions))
for (i in 1:n_states_all) {
  
  
  
}

