##############################################################################################
################# New script and folder for pHMM type analyses and tests #####################
##############################################################################################
funcSrc.path <- "/src/"
funcSrc.list <- list.files(path=paste('.',funcSrc.path,sep=''), recursive=T)
sapply(funcSrc.list, function(i){source(paste('.',funcSrc.path,'/',i,sep=''))})
#######################################################################################################
nk = 5 #nr. of features
nm = 4 #nr. of states
np = 27 #nr. of participants
#######################################################################################################
# Generate TPM (same all over population)
TPM <- rbind(c(0.95, 0.05, 0,   0),
             c(0.1, 0.85, 0,   0.05),
             c(0.05, 0,   0.85, 0.1),
             c(0  , 0,   0.05, 0.95))
# Generate population mus (participant mus will be drawn from this)
popmus <- rbind(c(-0.5, 0.5,  0.5, -0.5), #BV
                c(1,    0.2, -0.2, -1),   #ApEn
                c(0,   -0.2, -0.2,  0),   #BV x ApEn
                c(-0.4, 0.4,  0.4, -0.4), #Pupil baseline
                c(0,   -0.2, -0.2,  0))   #Generic goal focus measure
# Generate standard deviation of population mus (participant mus will be drawn from this)
popmuvar <- rbind(c(0.2, 0.1, 0.1, 0.1), #BV
                  c(0.1, 0.1, 0.1, 0.1), #ApEn
                  c(0.1, 0.1, 0.1, 0.1), #BV x ApEn
                  c(0.2, 0.1, 0.1, 0.1), #Pupil baseline
                  c(0.1, 0.1, 0.2, 0.1)) #Generic goal focus measure
# Generate half-cauchy deviation for standard deviation of emissions ## make a bit smaller
popsdvar <- rbind(c(2, 2, 2, 2), #BV
                  c(2, 2, 2, 2), #ApEn
                  c(2, 2, 2, 2), #BV x ApEn
                  c(2, 2, 2, 2), #Pupil baseline
                  c(2, 2, 2, 2)) #Generic goal focus measure
# Generate participant standard deviations for emissions
sds      <- rbind(c(1.1, 1.1, 1.1, 1.1), #BV
                  c(1.1, 1.1, 1.1, 1.1), #ApEn
                  c(1.1, 1.1, 1.1, 1.1), #BV x ApEn
                  c(1.1, 1.1, 1.1, 1.1), #Pupil baseline
                  c(1.1, 1.1, 1.1, 1.1)) #Generic goal focus measure
sds <- array(rep(sds, np), dim=c(nk, nm, np))
# Draw participant mus from N(popmus, popmuvar)
mus <- hierarchical.mu(popmus,popmuvar, np)
#######################################################################################################
# Generate a data set with tripartite PHMM characteristics
data.list <- lapply(1:np, function(pp){PHMM.generate(nt=1200, pIdx=seq(80,1200,80), nk=nk, nm=nm, TPM=TPM, mus=mus[,,pp], initState=1)})
X.list <- lapply(1:length(data.list), function(p){data.list[[p]]$X})

# Get parameters into a vector for optimization
params <- c(log(TPM[1]/(1-TPM[1])), log(TPM[2]/(1-TPM[2]-TPM[14])), log(TPM[3]/(1-TPM[3]-TPM[15])),
            log(TPM[14]/(1-TPM[2]-TPM[14])), log(TPM[15]/(1-TPM[3]-TPM[15])), log(TPM[16]/(1-TPM[16])),
            popmus, popmuvar, popsdvar, mus, sds)

# Test likelihood computations
ML.tripartite.PHMM.wrap(params, X.list, G=P.Gaus)



