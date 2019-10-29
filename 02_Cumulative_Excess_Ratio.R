
###
rm(list = ls());
###
library(SemiCompRisks)
library(gaussquad)

dyn.load("RASCRF_PM_SM.so")
source("00_Functions.R")

### The number of nodes for GH
K <- 3
L <- 3
M <- 3

rule3 <- hermite.h.quadrature.rules(K)
rule4 <- hermite.h.quadrature.rules(L)
rule5 <- hermite.h.quadrature.rules(M)

x.k <- rule3[[K]]$x
w.k <- rule3[[K]]$w
x.l <- rule4[[L]]$x
w.l <- rule4[[L]]$w
x.m <- rule5[[M]]$x
w.m <- rule5[[M]]$w

### Loading data
load("data.RData")

### Output from PEM-MVN model fit
load("Output_Fit_PEM_MVN/fit.RData")
load("Output_Fit_PEM_MVN/gammaPch1.RData") ## make sure to unzip the "gam" zip files first!
load("Output_Fit_PEM_MVN/V1Pch1.RData")
load("Output_Fit_PEM_MVN/V2Pch1.RData")
load("Output_Fit_PEM_MVN/V3Pch1.RData")

form <- fit$setup$Formula
X <- X1 <- X2 <- X3 <- model.frame(formula(form, lhs=0, rhs=1), data=data)
cluster <- data$cluster

#########
### Cumulative excess readmission ratio (theta1)

## adjusted cumulative readmission rate
muA1 <- CER.A1(fit, gamma.p, V1.p, V2.p, X1, X2, cluster)

## standardized adjusted cumulative readmission rate
muS1 <- CER.S1(fit, gamma.p, X1, X2, cluster, K, L, x.k, w.k, x.l, w.l)

theta1 <- muA1/muS1
save(theta1, file = "theta1.RData")


#########
### Cumulative excess mortality ratio (theta2)

## adjusted cumulative mortality rate
muA2 <- CER.A2(fit, gamma.p, V1.p, V2.p, V3.p, X1, X2, X3, cluster)

## standardized adjusted cumulative mortality rate
muS2 <- CER.S2(fit, gamma.p, X1, X2, X3, cluster, K, L, M, x.k, w.k, x.l, w.l, x.m, w.m)

theta2 <- muA2/muS2
save(theta2, file = "theta2.RData")


