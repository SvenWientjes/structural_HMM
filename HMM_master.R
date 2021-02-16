##############################################################################################
################# New script and folder for pHMM type analyses and tests #####################
##############################################################################################
funcSrc.path <- "/src/utility/"
funcSrc.list <- list.files(path=paste('.',funcSrc.path,sep=''))
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



