
####
####
CER.A1 <- function(fit, gamma.p, V1.p, V2.p, X1, X2, cluster)
{
    n <- dim(X1)[1]
    p1 <- dim(X1)[2]
    p2 <- dim(X2)[2]
    J <- dim(V1.p)[2]
    
    nj <- rep(NA, J); for(i in 1:J) nj[i]  <- length(which(cluster == i))
    
    nS <- dim(V1.p)[1]
    s1 <- fit$chain1$time_lambda1
    s2 <- fit$chain1$time_lambda2
    s3 <- fit$chain1$time_lambda3
    s1_diff <- diff(c(0, s1))
    s2_diff <- diff(c(0, s2))
    s3_diff <- diff(c(0, s3))
    K1 <- length(s1)
    K2 <- length(s2)
    K3 <- length(s3)
    
    muA1 <- matrix(NA, nS, J)
    
    for(MM in 1:nS)
    {
        if(MM %% 10 == 0) cat("MM=", MM, "/", nS, fill =T); Sys.time()
        
        beta1 <- fit$chain1$beta1.p[MM,]
        beta2 <- fit$chain1$beta2.p[MM,]
        lam1        <- fit$chain1$lambda1.fin[MM,]
        lam2        <- fit$chain1$lambda2.fin[MM,]
        gam     <- gamma.p[MM, ]
        V1      <- V1.p[MM, ]
        V2      <- V2.p[MM, ]
        
        mcmc <- .C("RASCRF_PM_SM_R_muA",
        n        = as.integer(n),
        p1        = as.integer(p1),
        p2        = as.integer(p2),
        J        = as.integer(J),
        nj        = as.double(nj),
        X1mat         = as.double(as.matrix(X1)),
        X2mat         = as.double(as.matrix(X2)),
        clust = as.double(cluster),
        be1 = as.double(beta1),
        be2 = as.double(beta2),
        lam1_vec = as.double(lam1),
        lam2_vec = as.double(lam2),
        s1_vec = as.double(s1),
        s2_vec = as.double(s2),
        K1 = as.integer(K1),
        K2 = as.integer(K2),
        gam = as.double(gam),
        V1vec = as.double(V1),
        V2vec = as.double(V2),
        muA    = as.double(rep(0, J)))
        
        muA1[MM, ] <- mcmc$muA
    }
    return(muA1)
}


####
####
CER.S1 <- function(fit, gamma.p, X1, X2, cluster, K, L, x.k, w.k, x.l, w.l)
{
    n <- dim(X1)[1]
    p1 <- dim(X1)[2]
    p2 <- dim(X2)[2]
    J <- dim(V1.p)[2]
    
    nj <- rep(NA, J); for(i in 1:J) nj[i]  <- length(which(cluster == i))
    
    nS <- dim(V1.p)[1]
    s1 <- fit$chain1$time_lambda1
    s2 <- fit$chain1$time_lambda2
    s3 <- fit$chain1$time_lambda3
    s1_diff <- diff(c(0, s1))
    s2_diff <- diff(c(0, s2))
    s3_diff <- diff(c(0, s3))
    K1 <- length(s1)
    K2 <- length(s2)
    K3 <- length(s3)
    
    ###
    muS1 <- matrix(NA, nS, J)
    
    for(MM in 1:nS)
    {
        if(MM %% 10 == 0) cat("MM=", MM, "/", nS, fill =T); Sys.time()
        
        beta1 <- fit$chain1$beta1.p[MM,]
        beta2 <- fit$chain1$beta2.p[MM,]
        lam1        <- fit$chain1$lambda1.fin[MM,]
        lam2        <- fit$chain1$lambda2.fin[MM,]
        gam     <- gamma.p[MM, ]
        Sigma_V <- fit$chain1$Sigma_V.p[,,MM]
        
        mcmc <- .C("RASCRF_PM_SM_R_muS",
        n        = as.integer(n),
        p1        = as.integer(p1),
        p2        = as.integer(p2),
        J        = as.integer(J),
        nj        = as.double(nj),
        X1mat         = as.double(as.matrix(X1)),
        X2mat         = as.double(as.matrix(X2)),
        clust = as.double(cluster),
        be1 = as.double(beta1),
        be2 = as.double(beta2),
        lam1_vec = as.double(lam1),
        lam2_vec = as.double(lam2),
        s1_vec = as.double(s1),
        s2_vec = as.double(s2),
        K1 = as.integer(K1),
        K2 = as.integer(K2),
        Sig_V = as.double(Sigma_V),
        gam = as.double(gam),
        K   = as.integer(K),
        L   = as.integer(L),
        x_k_vec = as.double(x.k),
        w_k_vec = as.double(w.k),
        x_l_vec = as.double(x.l),
        w_l_vec = as.double(w.l),
        muS    = as.double(rep(0, J)))
        
        muS1[MM, ] <- mcmc$muS
    }
    return(muS1)
}


####
####
CER.A2 <- function(fit, gamma.p, V1.p, V2.p, V3.p, X1, X2, X3, cluster)
{
    n <- dim(X1)[1]
    p1 <- dim(X1)[2]
    p2 <- dim(X2)[2]
    p3 <- dim(X3)[2]
    J <- dim(V1.p)[2]
    
    nj <- rep(NA, J); for(i in 1:J) nj[i]  <- length(which(cluster == i))
    
    nS <- dim(V1.p)[1]
    s1 <- fit$chain1$time_lambda1
    s2 <- fit$chain1$time_lambda2
    s3 <- fit$chain1$time_lambda3
    s1_diff <- diff(c(0, s1))
    s2_diff <- diff(c(0, s2))
    s3_diff <- diff(c(0, s3))
    K1 <- length(s1)
    K2 <- length(s2)
    K3 <- length(s3)
    
    muA2 <- matrix(NA, nS, J)
    
    for(MM in 1:nS)
    {
        if(MM %% 10 == 0) cat("MM=", MM, "/", nS, fill =T); Sys.time()
        
        beta1 <- fit$chain1$beta1.p[MM,]
        beta2 <- fit$chain1$beta2.p[MM,]
        beta3 <- fit$chain1$beta3.p[MM,]
        lam1        <- fit$chain1$lambda1.fin[MM,]
        lam2        <- fit$chain1$lambda2.fin[MM,]
        lam3        <- fit$chain1$lambda3.fin[MM,]
        gam     <- gamma.p[MM, ]
        V1      <- V1.p[MM, ]
        V2      <- V2.p[MM, ]
        V3      <- V3.p[MM, ]
        
        mcmc <- .C("RASCRF_PM_SM_D_muA",
        n        = as.integer(n),
        p1        = as.integer(p1),
        p2        = as.integer(p2),
        p3        = as.integer(p3),
        J        = as.integer(J),
        nj        = as.double(nj),
        X1mat         = as.double(as.matrix(X1)),
        X2mat         = as.double(as.matrix(X2)),
        X3mat         = as.double(as.matrix(X3)),
        clust = as.double(cluster),
        be1 = as.double(beta1),
        be2 = as.double(beta2),
        be3 = as.double(beta3),
        lam1_vec = as.double(lam1),
        lam2_vec = as.double(lam2),
        lam3_vec = as.double(lam3),
        s1_vec = as.double(s1),
        s2_vec = as.double(s2),
        s3_vec = as.double(s3),
        K1 = as.integer(K1),
        K2 = as.integer(K2),
        K3 = as.integer(K3),
        gam = as.double(gam),
        V1vec = as.double(V1),
        V2vec = as.double(V2),
        V3vec = as.double(V3),
        muA    = as.double(rep(0, J)))
        
        muA2[MM, ] <- mcmc$muA
    }
    return(muA2)
}


####
####
CER.S2 <- function(fit, gamma.p, X1, X2, X3, cluster, K, L, M, x.k, w.k, x.l, w.l, x.m, w.m)
{
    n <- dim(X1)[1]
    p1 <- dim(X1)[2]
    p2 <- dim(X2)[2]
    p3 <- dim(X3)[2]
    J <- dim(V1.p)[2]
    
    nj <- rep(NA, J); for(i in 1:J) nj[i]  <- length(which(cluster == i))
    
    nS <- dim(V1.p)[1]
    s1 <- fit$chain1$time_lambda1
    s2 <- fit$chain1$time_lambda2
    s3 <- fit$chain1$time_lambda3
    s1_diff <- diff(c(0, s1))
    s2_diff <- diff(c(0, s2))
    s3_diff <- diff(c(0, s3))
    K1 <- length(s1)
    K2 <- length(s2)
    K3 <- length(s3)
    
    muS2 <- matrix(NA, nS, J)
    
    for(MM in 1:nS)
    {
        if(MM %% 10 == 0) cat("MM=", MM, "/", nS, fill =T); Sys.time()
        
        beta1 <- fit$chain1$beta1.p[MM,]
        beta2 <- fit$chain1$beta2.p[MM,]
        beta3 <- fit$chain1$beta3.p[MM,]
        lam1        <- fit$chain1$lambda1.fin[MM,]
        lam2        <- fit$chain1$lambda2.fin[MM,]
        lam3        <- fit$chain1$lambda3.fin[MM,]
        gam     <- gamma.p[MM, ]
        Sigma_V <- fit$chain1$Sigma_V.p[,,MM]
        
        mcmc <- .C("RASCRF_PM_SM_D_muS",
        n        = as.integer(n),
        p1        = as.integer(p1),
        p2        = as.integer(p2),
        p3        = as.integer(p3),
        J        = as.integer(J),
        nj        = as.double(nj),
        X1mat         = as.double(as.matrix(X1)),
        X2mat         = as.double(as.matrix(X2)),
        X3mat         = as.double(as.matrix(X3)),
        clust = as.double(cluster),
        be1 = as.double(beta1),
        be2 = as.double(beta2),
        be3 = as.double(beta3),
        lam1_vec = as.double(lam1),
        lam2_vec = as.double(lam2),
        lam3_vec = as.double(lam3),
        s1_vec = as.double(s1),
        s2_vec = as.double(s2),
        s3_vec = as.double(s3),
        K1 = as.integer(K1),
        K2 = as.integer(K2),
        K3 = as.integer(K3),
        Sig_V = as.double(Sigma_V),
        gam = as.double(gam),
        K   = as.integer(K),
        L   = as.integer(L),
        M   = as.integer(M),
        x_k_vec = as.double(x.k),
        w_k_vec = as.double(w.k),
        x_l_vec = as.double(x.l),
        w_l_vec = as.double(w.l),
        x_m_vec = as.double(x.m),
        w_m_vec = as.double(w.m),
        muS    = as.double(rep(0, J)))
        
        muS2[MM, ] <- mcmc$muS
    }
    return(muS2)
}


####
####
## Note: change the values for "split.number" and "nSplits"
## only when spliting the computation to multiple jobs
## (e.g. using cluster)
L_01 <- function(theta, nT, split.number=1, nSplits=1000)
{
    J <- dim(theta)[2]
    nS <- dim(theta)[1]
    
    theta.med <- apply(theta, 2, median)
    
    thetaR <- theta
    for(i in 1:dim(theta)[1])
    {
        thetaR[i,] <- rank(theta[i,])
    }
    
    ## a tight constraint only for illustration purpose
    excl <- order(theta.med, decreasing=T)[1:214]
    incl <- which(rank(theta.med) <= 15)
    
    Ind <- seq_len(dim(theta)[2])[-sort(c(incl, excl))]
    theta.sub <- theta[, -sort(c(incl, excl))]
    
    KK <- dim(theta.sub)[2]
    m <- nT-length(incl)
    nC <- choose(KK, m)
    x.vec <- seq_len(KK)
    
    nIter <- choose(KK, m) %/% (nSplits-1)
    nSplits.last <- choose(KK, m) %% (nSplits-1)
    
    param <- list()
    start.ind = 1
    
    param$e <- 0
    param$h <- m
    param$a <- x.vec[1:m]
    
    if(split.number == 1)
    {
        start.ind = 1
        end.ind = nIter
    }else if(split.number > 1 & split.number < nSplits)
    {
        start.ind = nIter * (split.number-1) + 1
        end.ind =  start.ind + nIter - 1
        param <- findConf(x.vec, m, 1, start.ind-1, param)
    }else if(split.number == nSplits)
    {
        start.ind = nIter * (split.number-1) + 1
        end.ind =  start.ind + nSplits.last - 1
        param <- findConf(x.vec, m, 1, start.ind-1, param)
    }
    
    results <- ExpPR_L01_uni(theta, x.vec, m, incl, excl, Ind, start.ind, end.ind, param)
    
    val <- rep(0, J)
    val[results$Rhat[1,]] <- 1
    
    return(val)
}


####
####
ExpPR_L01_uni <- function(theta.post, x.vec, m, incl, excl, Ind, start.ind, end.ind, param)
{
    conf = rep(0, m)
    
    if(start.ind == 1)
    {
        a = seq_len(m)
        e = 0
        h = m
    }else
    {
        a = param$a
        e = param$e
        h = param$h
    }
    
    ii = start.ind
    jj = end.ind
    
    KK = length(x.vec)
    
    nS = dim(theta.post)[1]
    
    n_incl <- length(incl)
    n_excl <- length(excl)
    
    sel.true <- matrix(NA, nS, m + n_incl)
    for(j in 1:nS)
    {
        sel.true[j,] <- sort(which(rank(theta.post[j,]) <= m + n_incl))
    }
    
    val <- .C("EPR_L01_uni",
    sel_true = as.double(as.matrix(sel.true)),
    x_vec = as.double(x.vec),
    nS = as.integer(nS),
    N = as.integer(KK),
    m = as.integer(m),
    n_incl = as.integer(n_incl),
    n_excl = as.integer(n_excl),
    incl = as.double(incl),
    Ind_vec = as.double(Ind),
    ii = as.double(ii),
    jj = as.double(jj),
    conf_vec = as.double(conf),
    a_vec = as.double(a),
    e = as.integer(e),
    h = as.integer(h),
    nMin = as.integer(1),
    minEPR = as.double(0),
    inx_mat = as.double(rep(0, 100*(m+n_incl))))
    
    Rhat <- matrix(val$inx_mat, nrow = 100, byrow = TRUE)
    
    for(i in 1:100)
    {
        if(Rhat[i,1] != 0)
        {
            Rhat[i,] <- sort(Rhat[i,])
        }
    }
    
    param <- list()
    param$e <- val$e
    param$h <- val$h
    param$a <- val$a_vec
    
    list(param=param, nMin = val$nMin, Rhat=Rhat, minEPR = val$minEPR)
}


####
####
L_01.LO <- function(theta, m)
{
    J <- dim(theta)[2]
    nS <- dim(theta)[1]
    
    L_01 <- rep(0, J)
    
    theta.org <- theta
    
    sel.true <- matrix(NA, nS, m)
    for(j in 1:nS) sel.true[j,] <- which(rank(theta.org[j,]) <= m)
    
    temp <- matrix(0, nS, J)
    for(i in 1:nS) temp[i,] <- as.numeric(theta.org[i,] < sort(theta.org[i,], decreasing = F)[m])
    tempS <- apply(temp, 2, sum)
    
    ind <- sort(tempS)[J:1][m]
    R.prop <- c(1:J)[tempS >= ind]
    
    summ = 0
    for(k in 1:nS)
    {
        Rtr <- sel.true[k,]
        summ = summ + 2*(m-length(intersect(Rtr, R.prop)))/(J)
    }
    
    L_01.PEL <- summ/nS
    L_01[R.prop] <- 1
    
    list(PEL = L_01.PEL, val = L_01)
}


####
####
L_bc <- function(thata1, theta2)
{
    J <- dim(theta1)[2]
    nS <- dim(theta1)[1]
    
    theta1.med <- apply(theta1, 2, median)
    theta2.med <- apply(theta2, 2, median)
    
    BC <- rep(0, J)
    
    phi.true.mat <- matrix(NA, nS, J)
    for(i in 1:nS)
    {
        for(j in 1:J)
        {
            phi.true.mat[i,j] <- assign.phi(theta1[i, j], theta2[i, j])
        }
    }
    
    Cal.case <- Gen_CaseMat(10, theta1.med, theta2.med) ## no constraint
    case.mat <- Cal.case$case.mat
    nCase.vec <- Cal.case$nCase.vec
    
    phi.med <- rep(NA, J)
    for(j in 1:J)
    {
        phi.med[j] <- assign.phi(theta1.med[j], theta2.med[j])
    }
    
    phi.ini <- rep(1, J)
    
    results <- PEL_biv(phi.post=phi.true.mat, case.mat, nCase.vec, phi.ini=phi.ini)
    
    list(val = results$Phi.hat[1,], phi.med = phi.med)
}


####
####
Gen_CaseMat <- function(r, theta1, theta2)
{
    J <- length(theta1)
    
    case4 <- which((theta1-1)^2+(theta2-1)^2 <= r^2)
    case3 <- setdiff(which(theta1<=1+r & theta1>1-r & theta2<=1+r & theta2>1-r), case4)
    case2 <- setdiff(union(which(theta1<=1+r & theta1>1-r), which(theta2<=1+r & theta2>1-r)), union(case3, case4))
    case1 <- setdiff(1:J, union(union(case2, case3), case4))
    
    case3_123 <- intersect(which(theta1 <= 1 & theta2 > 1), case3)
    case3_234 <- intersect(which(theta1 <= 1 & theta2 <= 1), case3)
    case3_341 <- intersect(which(theta1 > 1 & theta2 <= 1), case3)
    case3_412 <- intersect(which(theta1 > 1 & theta2 > 1), case3)
    
    case2_12 <- intersect(which(theta2 > 1+r), case2)
    case2_23 <- intersect(which(theta1 <= 1-r), case2)
    case2_34 <- intersect(which(theta2 <= 1-r), case2)
    case2_41 <- intersect(which(theta1 > 1+r), case2)
    
    case1_1 <- intersect(which(theta1 > 1 & theta2 > 1), case1)
    case1_2 <- intersect(which(theta1 <= 1 & theta2 > 1), case1)
    case1_3 <- intersect(which(theta1 <= 1 & theta2 <= 1), case1)
    case1_4 <- intersect(which(theta1 > 1 & theta2 <= 1), case1)
    
    nCase.vec <- rep(NA, J)
    
    nCase.vec[case1] <- 1
    nCase.vec[case2] <- 2
    nCase.vec[case3] <- 3
    nCase.vec[case4] <- 4
    
    case.mat <- matrix(0, 4, J)
    for(i in case4) case.mat[,i] <-  1:4
    for(i in case3_123) case.mat[1:3,i] <-  1:3
    for(i in case3_234) case.mat[1:3,i] <-  c(2,3,4)
    for(i in case3_341) case.mat[1:3,i] <-  c(1,3,4)
    for(i in case3_412) case.mat[1:3,i] <-  c(1,2,4)
    for(i in case2_12) case.mat[1:2,i] <-  1:2
    for(i in case2_23) case.mat[1:2,i] <-  c(2,3)
    for(i in case2_34) case.mat[1:2,i] <-  c(3,4)
    for(i in case2_41) case.mat[1:2,i] <-  c(1,4)
    for(i in case1_1) case.mat[1,i] <-  1
    for(i in case1_2) case.mat[1,i] <-  2
    for(i in case1_3) case.mat[1,i] <-  3
    for(i in case1_4) case.mat[1,i] <-  4
    
    list(case.mat=case.mat, nCase.vec=nCase.vec, case1=case1, case2=case2, case3=case3, case4=case4)
}


####
####
assign.phi <- function(theta.X, theta.Y)
{
    if(theta.X > 1 & theta.Y > 1) phi = 1
    else if(theta.X <= 1 & theta.Y > 1) phi = 2
    else if(theta.X <= 1 & theta.Y <= 1) phi = 3
    else if(theta.X > 1 & theta.Y <= 1) phi = 4
    
    return(phi)
}


####
####
Gen_CaseMat <- function(r, theta1, theta2)
{
    J <- length(theta1)
    
    case4 <- which((theta1-1)^2+(theta2-1)^2 <= r^2)
    case3 <- setdiff(which(theta1<=1+r & theta1>1-r & theta2<=1+r & theta2>1-r), case4)
    case2 <- setdiff(union(which(theta1<=1+r & theta1>1-r), which(theta2<=1+r & theta2>1-r)), union(case3, case4))
    case1 <- setdiff(1:J, union(union(case2, case3), case4))
    
    
    case3_123 <- intersect(which(theta1 <= 1 & theta2 > 1), case3)
    case3_234 <- intersect(which(theta1 <= 1 & theta2 <= 1), case3)
    case3_341 <- intersect(which(theta1 > 1 & theta2 <= 1), case3)
    case3_412 <- intersect(which(theta1 > 1 & theta2 > 1), case3)
    
    case2_12 <- intersect(which(theta2 > 1+r), case2)
    case2_23 <- intersect(which(theta1 <= 1-r), case2)
    case2_34 <- intersect(which(theta2 <= 1-r), case2)
    case2_41 <- intersect(which(theta1 > 1+r), case2)
    
    case1_1 <- intersect(which(theta1 > 1 & theta2 > 1), case1)
    case1_2 <- intersect(which(theta1 <= 1 & theta2 > 1), case1)
    case1_3 <- intersect(which(theta1 <= 1 & theta2 <= 1), case1)
    case1_4 <- intersect(which(theta1 > 1 & theta2 <= 1), case1)
    
    nCase.vec <- rep(NA, J)
    
    nCase.vec[case1] <- 1
    nCase.vec[case2] <- 2
    nCase.vec[case3] <- 3
    nCase.vec[case4] <- 4
    
    case.mat <- matrix(0, 4, J)
    for(i in case4) case.mat[,i] <-  1:4
    for(i in case3_123) case.mat[1:3,i] <-  1:3
    for(i in case3_234) case.mat[1:3,i] <-  c(2,3,4)
    for(i in case3_341) case.mat[1:3,i] <-  c(1,3,4)
    for(i in case3_412) case.mat[1:3,i] <-  c(1,2,4)
    for(i in case2_12) case.mat[1:2,i] <-  1:2
    for(i in case2_23) case.mat[1:2,i] <-  c(2,3)
    for(i in case2_34) case.mat[1:2,i] <-  c(3,4)
    for(i in case2_41) case.mat[1:2,i] <-  c(1,4)
    for(i in case1_1) case.mat[1,i] <-  1
    for(i in case1_2) case.mat[1,i] <-  2
    for(i in case1_3) case.mat[1,i] <-  3
    for(i in case1_4) case.mat[1,i] <-  4
    
    list(case.mat=case.mat, nCase.vec=nCase.vec, case1=case1, case2=case2, case3=case3, case4=case4)
}


####
####
PEL_biv <- function(phi.post, case.mat, nCase.vec, phi.ini)
{
    nS = dim(phi.post)[1]
    J = dim(phi.post)[2]
    
    val <- .C("PEL_biv_C",
    phi_post = as.double(as.matrix(phi.post)),
    case_mat = as.double(as.matrix(case.mat)),
    nCase_vec = as.double(nCase.vec),
    phi_ini = as.double(phi.ini),
    nS = as.integer(nS),
    J = as.integer(J),
    nMin = as.integer(1),
    minPEL = as.double(0),
    iniPEL = as.double(0),
    inx_mat = as.double(rep(0, 100*(J))),
    nUpdate = as.double(rep(0, J)))
    
    Phi.hat <- matrix(val$inx_mat, nrow = 100, byrow = TRUE)
    
    list(nMin = val$nMin, Phi.hat=Phi.hat, minPEL = val$minPEL, iniPEL = val$iniPEL, nUpdate=val$nUpdate)
}

