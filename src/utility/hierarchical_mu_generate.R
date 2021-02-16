#### Function for hierarchically creating participant means ####
hierarchical.mu <- function(popmus, popsds, npp){
  ######################################################################################
  # popmus = K x M matrix of population-level means of a feature (K is nr. of data features; M is nr. of hidden states)
  # popsds = K x M matrix of std.devs of participant means drawn from population distribution
  # npp    = nr. of simulated participants ergo unique mu-matrices generated
  # ++++
  # mus = K x M x npp array of participant means
  ######################################################################################
  mus <- array(dim=c(nrow(popmus),ncol(popmus),npp))
  for(pp in 1:npp){
    mus[,,pp] <- sapply(1:length(popmus), function(mu){rnorm(1, mean=popmus[mu],sd=popsds[mu])})
  }
  return(mus)
}