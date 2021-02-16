#### compute a diagonal matrix of emission probabilities under hidden state parameters ####
P.Gaus <- function(x, mus, sds){
  probs <- c()
  for(i in 1:ncol(mus)){
    probs[i] <- mvtnorm::dmvnorm(x, mean=mus[,i], sigma=diag(sds[,i]))
  }
  return(diag(probs))
}