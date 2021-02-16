#### Function for generating tripartite PHMM data ####
PHMM.generate <- function(nt, pIdx, nk, nm, TPM, mus, initState=FALSE, delta=FALSE){
  ######################################################################################
  # nt    = nr. of trials in total
  # pIdx  = vector of trials that are probed (time indices)
  # nk    = nr. of features (All Gaussian distributed)
  # nm    = nr. of hidden states
  # TPM   = M x M matrix of transition probabilities
  # mus   = KxM matrix of means for the features in the different hidden states
  #   - Note that standard deviation will always be 1. Mu thus defines the 'effect size'
  # ++++
  # X    = T x K matrix of generated emission data
  # P    = vector of probe values at pIdx
  # pIdx = identical to input pIdx
  # hMC  = vector sequence of hidden states (hidden Markov Chain)
  ######################################################################################
  
  # Initialize hidden markov chain
  hMC <- rep(0, nt)
  
  # Set the first state of the chain
  if(initState){
    hMC[1] <- initState
  }else if(!initState & delta){
    hMC[1] <- sample(x=1:nm, size=1, prob=delta)
  }else{
    stop('Define either the initial hidden state (initState) or a vector defining the starting state probabilities (delta).')
  }
  
  # Generate the rest of the chain
  for(t in 2:nt){
    hMC[t] <- sample(x=1:nm, size=1, prob=TPM[hMC[t-1],])
  }
  
  # Initialize matrix of emission data
  X <- matrix(nrow=nt, ncol=nk)
  
  # Fill matrix of emission data
  for(t in 1:nt){
    X[t,] <- sapply(mus[,hMC[t]],function(mu){rnorm(1,mu,sd=1)})
  }
  
  # Get probe values with indices
  P <- hMC[pIdx]
  
  return(list(X=X, P=P, pIdx=pIdx, hMC=hMC))
}