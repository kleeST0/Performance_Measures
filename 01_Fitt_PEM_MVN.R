###
rm(list = ls());

###
set.seed(1234)
###
library(SemiCompRisks)
load("data.RData")

###
id=data$cluster
form <- Formula(time1 + event1 | time2 + event2 ~ deyo2_+race_noWhite+adm+age_std+sex_female+los_std+disc_homecare+disc_hospice+disc_snf_icf+disc_other | deyo2_+race_noWhite+adm+age_std+sex_female+los_std+disc_homecare+disc_hospice+disc_snf_icf+disc_other | deyo2_+race_noWhite+adm+age_std+sex_female+los_std+disc_homecare+disc_hospice+disc_snf_icf+disc_other)


## Subject-specific frailty variance component
##  - prior parameters for 1/theta
##
theta.ab <- c(0.7, 0.7)

## PEM baseline hazard function
##
PEM.ab1 <- c(0.7, 0.7) # prior parameters for 1/sigma_1^2
PEM.ab2 <- c(0.7, 0.7) # prior parameters for 1/sigma_2^2
PEM.ab3 <- c(0.7, 0.7) # prior parameters for 1/sigma_3^2
##
PEM.alpha1 <- 10 # prior parameters for K1
PEM.alpha2 <- 10 # prior parameters for K2
PEM.alpha3 <- 10 # prior parameters for K3

## MVN cluster-specific random effects
##
Psi_v <- diag(1, 3)
rho_v <- 5

##
hyperParams <- list(theta=theta.ab, PEM=list(PEM.ab1=PEM.ab1, PEM.ab2=PEM.ab2, PEM.ab3=PEM.ab3, PEM.alpha1=PEM.alpha1, PEM.alpha2=PEM.alpha2, PEM.alpha3=PEM.alpha3), MVN=list(Psi_v=Psi_v, rho_v=rho_v))

###################
## MCMC SETTINGS ##
###################

## Setting for the overall run
##
numReps    <- 2000000
thin       <- 1000
burninPerc <- 0.5

## Settings for storage
##
nGam_save <- dim(data)[1]
storeV    <- rep(TRUE, 3)

## Tuning parameters for specific updates
##
##  - those common to all models
mhProp_theta_var  <- 0.05
mhProp_Vg_var     <- c(0.05, 0.05, 0.05)
##
## - those specific to the PEM specification of the baseline hazard functions
Cg        <- c(0.2, 0.2, 0.2)
delPertg  <- c(0.5, 0.5, 0.5)
rj.scheme <- 1
Kg_max    <- c(50, 50, 50)
sg_max    <- c(max(data$time1[data$event1 == 1]),
max(data$time2[data$event1 == 0 & data$event2 == 1]),
max(data$time2[data$event1 == 1 & data$event2 == 1]))

time_lambda1 <- seq(1, sg_max[1], 1)
time_lambda2 <- seq(1, sg_max[2], 1)
time_lambda3 <- seq(1, sg_max[3], 1)               

##
mcmc.PEM <- list(run=list(numReps=numReps, thin=thin, burninPerc=burninPerc),
storage=list(nGam_save=nGam_save, storeV=storeV),
tuning=list(mhProp_theta_var=mhProp_theta_var,
mhProp_Vg_var=mhProp_Vg_var, Cg=Cg, delPertg=delPertg,
rj.scheme=rj.scheme, Kg_max=Kg_max,
time_lambda1=time_lambda1, time_lambda2=time_lambda2,
time_lambda3=time_lambda3))

#############
## PEM-MVN ##
#############

##
myModel <- c("semi-Markov", "PEM", "MVN")
myPath  <- "Output/05-Results-PEM_MVN/"

startValues      <- initiate.startValues_HReg(form, data, model=myModel, id, nChain=2)

##
fit_PEM_MVN <- BayesID_HReg(form, data, id, model=myModel,
hyperParams, startValues, mcmc.PEM, path=myPath)

fit_PEM_MVN
summ.fit_PEM_MVN <- summary(fit_PEM_MVN); names(summ.fit_PEM_MVN)
summ.fit_PEM_MVN
pred_PEM_MVN <- predict(fit_PEM_MVN)
plot(pred_PEM_MVN, plot.est="Haz")
plot(pred_PEM_MVN, plot.est="Surv")








