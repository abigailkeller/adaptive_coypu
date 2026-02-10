
coypus <- readRDS("counting-rodents/analyses/data/coypus.rds")

net <- read.csv("netherlands_data/occurrence.csv")


# code
code <- nimbleCode({
  for (i in 1:K) { # loop over sites
    
    # prob of cells
    pi[i,1] <- p[i,1] # survey 1
    for(j in 2:J){ # loop over surveys > 1
      pi[i,j] <- prod(1 - p[i, 1:(j-1)]) * p[i,j]
    }
    pcap[i] <- sum(pi[i,1:J])
    
    # capture prob
    for(j in 1:J){ # loop over surveys
      cloglog(p[i,j]) <- alpha[i, j]
      pic[i,j] <- pi[i,j] / pcap[i]
      alpha[i, j] ~ dnorm(0, sd = 1.5)
    }
    
    # likelihood: observation
    y[i,1:J] ~ dmulti(pic[i,1:J], n[i]) # multinomial for each site
    n[i] ~ dbin(pcap[i], N[i]) # for each site
    # likelihood: process
    N[i] ~ dpois(lambda[i]) # for each site
    log(lambda[i]) <- beta[1] + beta[2] * temp[i]
  } # end of loop over sites
  
  for (k in 1:2){
    beta[k] ~ dnorm(0, sd = 1.5)
  }
  
  #---------------- check model goodness of fit
  
  for(i in 1:K){ 
    ## replicate data set
    n.pred[i] ~ dbin(pcap[i], N[i])
    y.pred[i,1:J] ~ dmulti(pic[i,1:J], n.pred[i])
    
    ## Freeman–Tukey residuals: observation component
    for(k in 1:J){ 
      e1[i,k]      <- pic[i,k] * n[i] # expected counts (obs)
      resid1[i,k]  <- pow(sqrt(y[i,k]) - sqrt(e1[i,k]), 2)
      e1.pred[i,k] <- pic[i,k] * n.pred[i] # expected counts (rep)
      resid1.pred[i,k] <- pow(sqrt(y.pred[i,k]) - sqrt(e1.pred[i,k]), 2)
    }  
    ## Freeman–Tukey residuals: abundance component
    e2[i] <- pcap[i] * lambda[i] 
    resid2[i] <- pow(sqrt(n[i]) - sqrt(e2[i]), 2) 
    resid2.pred[i] <- pow(sqrt(n.pred[i]) - sqrt(e2[i]), 2) 
  }  
  
  # fit statistic
  fit1.data <- sum(resid1[1:K, 1:J])      # observation part (data)
  fit1.pred <- sum(resid1.pred[1:K, 1:J]) # observation part (replicates)
  fit2.data <- sum(resid2[1:K])           # abundance part (data)
  fit2.pred <- sum(resid2.pred[1:K])      # abundance part (replicates)
  
})