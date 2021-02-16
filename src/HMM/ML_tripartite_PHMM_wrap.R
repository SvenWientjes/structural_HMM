#### Function for recovering the input parameter vector from an optim call and parsing it to the multilevel tripartite HMM likelihood function ####
ML.tripartite.PHMM.wrap <- function(params, X.list, G, negative=T, reversed=F){
  ######################################################################################
  # params = vector of unconstrained (?) parameter values that will be transformed to tripartite.PHMM parameters
  # X.list = list of participant T x K matrices of emission data
  # P.list = list of participant vectors of probed state values
  #   - 0 is on-task, 1 is median, 2 is off-task
  # pIdx.list = list of participant vectors of trial indices matching the probes from P
  # G    = function object that returns emission probabilities as M x M diagonal matrix
  # negative = logical argument if log-likelihood should be return negative. Optim requires this, so default = T 
  # ++++
  # l = computed log-likelihood of the data under the given parameters. Negative LL if argument negative = T
  ######################################################################################
  # Transform Transition Probability Matrix (should this become transposed??)
  TPM <- matrix(0, nrow=4, ncol=4)
  TPM[1,1] <- exp(params[1])/sum(1,exp(params[1]))
  TPM[1,2] <- 1/sum(1,exp(params[1]))
  TPM[2,1] <- exp(params[2])/sum(1,exp(params[2]),exp(params[3]))
  TPM[2,2] <- 1/sum(1,exp(params[2]),exp(params[3]))
  TPM[2,4] <- exp(params[3])/sum(1,exp(params[2]),exp(params[3]))
  TPM[3,1] <- exp(params[4])/sum(1,exp(params[4]),exp(params[5]))
  TPM[3,3] <- 1/sum(1,exp(params[4]),exp(params[5]))
  TPM[3,4] <- exp(params[5])/sum(1,exp(params[4]),exp(params[5]))
  TPM[4,4] <- exp(params[6])/sum(1,exp(params[6]))
  TPM[4,3] <- 1/sum(1,exp(params[6]))
  
  # Recover popmus
  popmus <- matrix(params[7:26], nrow=5, ncol=4)
  # Recover popmuvar
  popmuvar <- matrix(params[27:46], nrow=5, ncol=4)
  # Recover popsdvar
  popsdvar <- matrix(params[47:66], nrow=5, ncol=4)
  # Recover mus
  mus <- array(params[67:(67+5*4*length(X.list)-1)], dim=c(5,4,length(X.list)))
  # Recover sds
  sds <- array(params[(67+5*4*length(X.list)):((67+5*4*length(X.list))+5*4*length(X.list)-1)], dim=c(5,4,length(X.list)))
  ML.tripartite.PHMM(X.list, P.list, pIdx.list, pr.reg, TPM, popmus, popmuvar, popsdvar, mus, sds, P.Gaus, negative, reversed)
}
