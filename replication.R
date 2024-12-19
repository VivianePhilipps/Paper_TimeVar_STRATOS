##################################################################
#### Replication script for the simulation part of the paper #####
##################################################################


### load R packages
library(splines)
library(mvtnorm)
library(survival)
library(lcmm)
library(JM) #installed with install_github("CecileProust-Lima/JM")


### define the seed (run this script multiple times with different seeds to do multiple replicates of a scenario)
seed <- 692308

### definition of the scenario ###

## change the 6 following values to switch between the different scenarios presented in the paper :
trend <- "linear" # or "quad"
asso <- 0.2       # weak associaiton
error <- 1        # low measurement error
surv <- "surv1"   # for high survival. Alternatives : "surv2" for medium survivial or "surv3" for low survival
pMCAR <- 0        # no intermittent missing data
corRE <- FALSE    # independent random effects

if(trend == "linear") {
    if(corRE == FALSE) {
        prmsim <- c(weib1 = 0.01, weib2 = 2, asso = 0.2, beta0 = 0, beta1 = 1, var0 = 1, cov = 0, var1 = 0.5, lin0 = 0, lin1 = 1, sigma = 1)
    } else {
        prmsim <- c(weib1 = 0.01, weib2 = 2, asso = 0.2, beta0 = 0, beta1 = 1, var0 = 2, cov = -0.12, var1 = 0.05, lin0 = 0, lin1 = 1, sigma = 2)
    }
} else {
    ## quadratic
    prmsim <- c(weib1 = 0.01, weib2 = 2, asso = 0.2, beta0 = 0, beta1 = -4, beta2 = 0.6, var0 = 1, cov01 = 0, var1 = 0.5, cov02 = 0, cov12 = 0, var2 = 0.1, lin0 = 0, lin1 = 1, sigma = 1)
}

prmsim["asso"] <- asso
prmsim["sigma"] <- error

if(surv == "surv1"){ prmsim["weib1"] <- 0.01; prmsim["weib2"] <- 1}
if(surv == "surv2"){ prmsim["weib1"] <- 0.01; prmsim["weib2"] <- 2}
if(surv == "surv3"){ prmsim["weib1"] <- 0.1; prmsim["weib2"] <- 1.5}


### get the simJM function to generate the data
source("simJM.R")


### drawRE function : draw random effects from their posterior distribution (used for MI method)
drawRE <- function(m)
{
    predRE <- m$predRE
    varRE <- m$varRE
    nea <- ncol(predRE)-1
    ns <- nrow(predRE)
    
    REdraw <- sapply(1:ns, function(i, nea, predRE, varRE){
        v <- matrix(0, nea, nea)
        v[upper.tri(v, diag = TRUE)] <- varRE[(i - 1) * nea * (nea + 1) / 2 + 1:(nea * (nea + 1) / 2)]
        v <- t(v)
        v[upper.tri(v,diag = TRUE)] <- varRE[(i - 1) * nea * (nea + 1) / 2 + 1:(nea * (nea + 1) / 2)]
        rmvnorm(1, mean = as.numeric(predRE[i, -1]), sigma = v)
    },
    nea = nea, predRE = predRE, varRE = varRE)
    
    return(data.frame(ID = predRE[, 1], t(REdraw)))
}

### bootRC function : bootstrap for two stage methods (RC, PE-RC and MI)
boot2step <- function(m, d, method = "RC")
{
    Bb <- rmvnorm(1, estimates(m), vcov(m))
    if(m$N[3]>0)
    {
        ch <- matrix(0, sum(m$idea), sum(m$idea))
        ch[upper.tri(ch, diag = TRUE)] <- Bb[m$N[1] + m$N[2] + 1:m$N[3]]
        v <- t(ch) %*% ch
        Bb[m$N[1] + m$N[2] + 1:m$N[3]] <- v[upper.tri(v, diag = TRUE)]
    }
    
    z <- m$call
    z$maxiter <- 0
    z$B <- as.numeric(Bb)
    z$verbose <- FALSE
    mb <- eval(z)

    if(method == "RC")
    {
        ## use predicted RE (RC, PE-RC)
        db <- merge(d, mb$predRE)
        old <- colnames(d)
        old[which(old == "T")] <- "t"
        colnames(db) <- c(old, paste("u", 0:(sum(m$idea)-1), sep = ""))
        dbtimes <- sort(unique(with(db, t)))
        db3 <- survSplit(Surv(T0,t,D)~., data = db, cut = dbtimes)
        
        EF <- model.matrix(as.formula(m$call$fixed[-2]), data = db3)
        EA <- model.matrix(as.formula(m$call$random), data = db3)
        
        db3$Yt <- EF %*% mb$best[1:m$N[2]] + apply(db3[, paste("u",0:(sum(m$idea)-1),sep = ""),drop = FALSE]*EA, 1, sum)
        
        mbrc <- coxph(Surv(T0,t,D)~Yt, data = db3, control = coxph.control(timefix = FALSE))

        res <- c(coef(mbrc), as.numeric(vcov(mbrc)))
    }
    else
    {
        ## draw RE from their conditional distribution (MI)
        REdraws <- drawRE(mb)
        db2 <- merge(d, REdraws)
        old <- colnames(d)
        old[which(old == "T")] <- "t"
        colnames(db2) <- c(old, paste("u",0:(sum(m$idea)-1),sep = ""))
        dbtimes2 <- sort(unique(with(db2, t)))
        db32 <- survSplit(Surv(T0,t,D)~., data = db2, cut = dbtimes2)
        db32$event2 <- db32$D 
        
        EF2 <- model.matrix(as.formula(m$call$fixed[-2]), data = db32)
        EA2 <- model.matrix(as.formula(m$call$random), data = db32)
        
        db32$Yt <- EF2 %*% mb$best[1:m$N[2]] + apply(db32[,paste("u",0:(sum(m$idea)-1),sep = ""),drop = FALSE]*EA2, 1, sum)
        
        mbrc2 <- coxph(Surv(T0,t,D)~Yt, data = db32, control = coxph.control(timefix = FALSE))

        res <- c(coef(mbrc2), as.numeric(vcov(mbrc2)))
    }
    
        return(res)
}



vcov.default <- function(object,...){return(object$vcov)}



### Data generation ###
if(trend == "linear"){
    dsimtot <- simJM(n = 500, t.ecart = 2, t.marge = 1, t.max = 10, 
                     beta = prmsim[3+1:2],
                     nRE = 2, B.random = prmsim[5 + 1:3],
                     sigma = prmsim[11],
                     seuils = list(linearY1 = prmsim[8 + 1:2]), 
                     weibull = prmsim[1:2], 
                     sharedtype = "value", association = prmsim[3],
                     seed = seed, NAafterEvent = FALSE, pMCAR = pMCAR)
} else {
        dsimtot <- simJM(n = 500, t.ecart = 2, t.marge = 1, t.max = 10, 
                         beta = prmsim[3 + 1:3],
                         nRE = 3, B.random = prmsim[6 + 1:6],
                         sigma = prmsim[15],
                         seuils = list(linearY1 = prmsim[12 + 1:2]), 
                         weibull = prmsim[1:2], 
                         sharedtype = "value", association = prmsim[3],
                         seed = seed, NAafterEvent = FALSE,
                         model = list(fixed = ~ t + I(t^2), random = ~ t + I(t^2)))
}


### rounding to avoid errors in jointModel
dsimtot$t <- ceiling(dsimtot$t*100000)/100000
dsimtot$T <- ceiling(dsimtot$T*10000)/10000
dsimtot$Y1 <- as.numeric(dsimtot$Y1)

### remove data after the event occurence
afterT <- which(dsimtot$t > dsimtot$T)
if(length(afterT)) dsim <- dsimtot[-afterT,]

### survival data
dsim0 <- dsim[!duplicated(dsim$ID), c("ID", "T0", "T", "D")]

### sample description
tmax <- max(dsim$t)
tdiff <- diff(dsim$t)
tdiff <- mean(tdiff[which(tdiff > 0)])
ns <- length(unique(dsim$ID))
nevt <- length(which(dsim0$D == 1))
nmes <- mean(as.numeric(table(dsim$ID[which(!is.na(dsim$Y1))])))
nmes0 <- mean(as.numeric(table(dsim$ID[which(!is.na(dsim$Y1) & dsim$D == 0)])), na.rm=TRUE)
nmes1 <- mean(as.numeric(table(dsim$ID[which(!is.na(dsim$Y1) & dsim$D == 1)])), na.rm=TRUE)


### LOCF ###

dsim1 <- tmerge(dsim0, dsim0, id = ID, tstart = T0, tstop = T, Dinterv = event(T,D)) 
dsim2 <- tmerge(dsim1, dsim[, c("ID", "t", "Y1")], id = ID, Yt = tdc(t,Y1))
        
mlocf <- coxph(Surv(tstart, tstop, Dinterv) ~ Yt, data = dsim2)        
##summary(mlocf)


### Regression calibration (RC) ###
if(trend=="linear"){
    mY1 <- hlme(Y1 ~ t, random = ~ t, subject = "ID", data = dsim, verbose = FALSE)
} else {
        mY1 <- hlme(Y1 ~ t + I(t^2), random = ~ t + I(t^2), subject = "ID", data = dsim, verbose = FALSE)            
}

if(mY1$conv==1){
    ## bootstrap for the first step
    boot_mrc <- replicate(500, boot2step(mY1, dsim0, method = "RC"))
    coefm <- mean(boot_mrc[1,]) # mean estimate
    varm <- mean(boot_mrc[2,])  # mean asymptotic variance
    vare <- var(boot_mrc[1,]) # empirical variance 
    v_mrc <- varm + (1+1/500)*vare # corrected variance
    
    mrc <- list(coefficients = coefm, vcov = v_mrc, conv = 1)   
} else {
    mrc <- list(coefficients = NA, vcov = NA, conv = mY1$conv)
}

        
### Regression calibration with post-event information (PE-RC)
if(trend == "linear"){
    mY1_ext <- hlme(Y1 ~ t, random = ~ t, subject = "ID", data = dsimtot, verbose = FALSE)
} else {
        mY1_ext <- hlme(Y1 ~ t + I(t^2), random = ~ t + I(t^2), subject = "ID", data = dsimtot, verbose = FALSE)
}

if(mY1_ext$conv==1){
    ## bootstrap 
    boot_mrc_ext <- replicate(500, boot2step(mY1_ext, dsim0, method = "RC"))
    coefm_ext <- mean(boot_mrc_ext[1,]) 
    varm_ext <- mean(boot_mrc_ext[2,])
    vare_ext <- var(boot_mrc_ext[1,])
    v_mrc_ext <- varm_ext + (1+1/500)*vare_ext
    
    mrc_ext <- list(coefficients = coefm_ext, vcov = v_mrc_ext, conv = 1) 
} else {
    mrc_ext <- list(coefficients = NA, vcov = NA, conv = mY1_ext$conv)
}

       
### Multiple imputation (MI)
cumhz <- basehaz(coxph(Surv(T,D) ~ 1, data = dsim0)) # cumulative hazard (Nelson-Aalen)
dsim00 <- merge(dsim0, cumhz, by.x = "T", by.y = "time", sort = FALSE)
dsim00 <- dsim00[,c("ID", "T0", "T", "D", "hazard")]

dsim4 <- merge(dsim, dsim00[, c("ID", "hazard")], sort = FALSE)
dsim4 <- dsim4[order(dsim4$ID, dsim4$t), ]
dsim4$event2 <- unlist(tapply(dsim4$D, dsim4$ID, function(x) {c(rep(0, length(x) - 1), x[length(x)])} ))

if(trend=="linear"){
    mY1_mijm <- hlme(Y1 ~ t + hazard + event2, random = ~ t, subject = "ID", data = dsim4, verbose = FALSE)
} else {
    mY1_mijm <- hlme(Y1 ~ t + I(t^2) + hazard + event2, random = ~ t + I(t^2), subject = "ID", data = dsim4, verbose = FALSE)
}

if(mY1_mijm$conv==1){
    boot_mijm <- replicate(500, boot2step(mY1_mijm, dsim00, method = "MI"))    
    coefm_mijm2 <- mean(boot_mijm[1,]) 
    varm_mijm2 <- mean(boot_mijm[2,])
    vare_mijm2 <- var(boot_mijm[1,]) 
    v_mrc_mijm2 <- varm_mijm2 + (1+1/500)*vare_mijm2
    
    mrc_mijm <- list(coefficients = coefm_mijm2, vcov = v_mrc_mijm2, conv = 1)
} else {
    mrc_mijm <- list(coefficients = NA, vcov = NA, conv = mY1_mijm$conv)
}



### Joint Model (JM)
if(trend=="linear"){
    mlme <- lme(Y1 ~ t, random = ~ 1 +t | ID, data = dsim)
    GHk <- 15
    jj <- 2
} else {
    GHk <- 9
    jj <- c(2, 3, 6)
        mlme <- lme(Y1 ~ t + I(t^2), random = ~ t + I(t^2) | ID, data = dsim, control = lmeControl(opt = "optim"))
}
if(!exists("mlme")){
    cat("lme failed \n")
    next;
}
        
mcox <- coxph(Surv(T,D) ~ 1, data = dsim0, x = TRUE)
mJM <- try(jointModel(mlme, mcox, timeVar = "t", method = "spline-PH-aGH", control = list(optimizer = "mla",GHk = GHk,lng.in.kn = 1,iter.EM = 0,iter.qN = 100,nproc = 1), parameterization = "value"), silent = FALSE)


if(inherits(mJM, "try-error")){ # try other initial values
    
    if(trend=="linear"){
        Binit <- list(betas=c(Intercept=0, slope=0), sigma=1, D=matrix(c(1,0,0,1), 2,2, dimnames=list(c("Intercept", "slope"),c("Intercept", "slope"))), alpha=0, gammas.bs=rep(0.1,5))
        BNA <- list(betas=c(Intercept=NA, slope=NA), sigma=NA, alpha=NA, gammas.bs=rep(NA,5), D=matrix(c(1,0,0,1), 2,2, dimnames=list(c("Intercept", "slope"),c("Intercept", "slope"))))
        jj <- 2
    } else {
        Binit <- list(betas=c(Intercept=0, t1=0, t2=0), sigma=1, D=matrix(c(1,0,0,0,1,0,0,0,1), 3,3, dimnames=list(c("Intercept", "t1", "t2"),c("Intercept", "t1", "t2"))), alpha=0, gammas.bs=rep(0.1,5)) 
        BNA <- list(betas=c(Intercept=NA, t1=NA, t2=NA), sigma=NA, alpha=NA, gammas.bs=rep(NA,5), D=matrix(c(1,0,0,0,1,0,0,0,1), 3,3, dimnames=list(c("Intercept", "t1", "t2"),c("Intercept", "t1", "t2")))) 
        jj <- c(2, 3, 6)
    }
    mJM <- try(jointModel(mlme, mcox, timeVar = "t", method = "spline-PH-aGH", control = list(optimizer = "mla", GHk = GHk, lng.in.kn = 1, iter.EM = 0, iter.qN = 100, nproc = 1), parameterization = "value", init = Binit))

    if(inherits(mJM, "try-error")) mJM <- list(convergence=99, coefficients=BNA, Hessian=matrix(NA, length(unlist(BNA))-length(jj), length(unlist(BNA))-length(jj)))
}

mjm <- list(conv = as.numeric(mJM$convergence == 0),
            coefficients = c(unlist(mJM$coefficients[1:4]), as.numeric(chol(mJM$coefficients$D))[-jj]),
            vcov = solve(mJM$Hessian))


cat("OK for seed=",seed,"\n")


### results
res <- matrix(c(as.numeric(mlocf$info[4]==0), coef(mlocf), sqrt(diag(vcov(mlocf))),
                mrc$conv, coef(mrc), sqrt(vcov(mrc)),
                mrc_ext$conv, coef(mrc_ext), sqrt(vcov(mrc_ext)),
                mrc_mijm$conv, coef(mrc_mijm), sqrt(vcov(mrc_mijm)),
                mjm$conv, coef(mjm), sqrt(diag(vcov(mjm))),
                seed, tmax, tdiff, ns, nevt, nmes, nmes0, nmes1), ncol = 1)




if(trend=="linear") {
noms <- c("convlocf", "prmlocf", "sdlocf", "convrc", "prmrc", "sdrc", "convrcext", "prmrcext", "sdrcext", "convrcmijm", "prmrcmijm", "sdrcmijm",  "convjm", "beta0", "beta1", "sigma", "alpha", "b1", "b2", "b3", "b4", "b5", "chol0", "chol01", "chol1", paste("sd",c("beta0", "beta1", "sigma", "alpha", "b1", "b2", "b3", "b4", "b5", "chol0", "chol01", "chol1"),sep=""), "seed", "tmax", "tdiff", "ns", "nevt", "nmes", "nmes0", "nmes1") ## trend=linear
} else {
noms <- c("convlocf", "prmlocf", "sdlocf", "convrc", "prmrc", "sdrc", "convrcext", "prmrcext", "sdrcext", "convrcmijm", "prmrcmijm","sdrcmijm", "convjm", "beta0", "beta1", "beta2", "sigma", "alpha", "b1", "b2", "b3", "b4", "b5", "chol0", "chol01", "chol1", "chol02", "chol12", "chol2", paste("sd",c("beta0", "beta1", "beta2", "sigma", "alpha", "b1", "b2", "b3", "b4", "b5", "chol0", "chol01", "chol1", "chol02", "chol12", "chol2"),sep=""), "seed", "tmax", "tdiff", "ns", "nevt", "nmes", "nmes0", "nmes1") ## trend=splines ou quad
}

if(length(noms) == nrow(res)) rownames(res) <- noms


### save
write.table(res,file=paste("simu_", trend, "_pMCAR", pMCAR, "_", surv, "_locf_rc_rcext_rcmijm_jm_prmsim_", paste(prmsim,collapse="_"), "_seed", seed, ".txt", sep = ""), sep = "\t")


quit("no")
