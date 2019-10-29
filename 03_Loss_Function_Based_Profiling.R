
###
rm(list = ls());
###
source("00_Functions.R")
dyn.load("PEL_biv.so")

load("theta1.RData")
load("theta2.RData")

############
### 1. Classifying hospitals according to whether they are ranked
###    in the top 10% for 90-day readmission and in the top 10% for 90-day mortality

J <- 264
nS <- 1000
nT <- 27 # the number of hospital in the top 10%

### Based on the minimizer of the approximate Bayes risks
##
L01.R <- L_01(theta1, nT)
L01.D <- L_01(theta2, nT)

# save(L01.R, file = "L01.R.RData")
# save(L01.D, file = "L01.D.RData")

# load(file = "L01.R.RData")
# load(file = "L01.D.RData")

## 
L01 <- rep(NA, J)
L01[L01.R == 0 & L01.D == 0] <- "No/No"
L01[L01.R == 0 & L01.D == 1] <- "No/Yes"
L01[L01.R == 1 & L01.D == 0] <- "Yes/No"
L01[L01.R == 1 & L01.D == 1] <- "Yes/Yes"

### Based on the posterior medians
PM.R <- PM.D <- rep(0, J)
PM.R[which(rank(apply(theta1, 2, median)) <= nT)] <- 1
PM.D[which(rank(apply(theta2, 2, median)) <= nT)] <- 1

PM <- rep(NA, J)
PM[PM.R == 0 & PM.D == 0] <- "No/No"
PM[PM.R == 0 & PM.D == 1] <- "No/Yes"
PM[PM.R == 1 & PM.D == 0] <- "Yes/No"
PM[PM.R == 1 & PM.D == 1] <- "Yes/Yes"

table(PM, L01)

############
### 2. Classifying hospitals according to whether they are found to have higher- or lower-
###    than expected 90-day readmission and 90-day mortality.
BC <- L_bc(thata1, theta2)

### Based on the minimizer of the approximate Bayes risks
L.BC <- BC$val
L.BC[L.BC == 1] <- "Higher/Higher"
L.BC[L.BC == 2] <- "Lower/Higher"
L.BC[L.BC == 3] <- "Lower/Lower"
L.BC[L.BC == 4] <- "Higher/Lower"

### Based on the posterior medians
PM.BC <- BC$phi.med
PM.BC[PM.BC == 1] <- "Higher/Higher"
PM.BC[PM.BC == 2] <- "Lower/Higher"
PM.BC[PM.BC == 3] <- "Lower/Lower"
PM.BC[PM.BC == 4] <- "Higher/Lower"

table(PM.BC, L.BC)
