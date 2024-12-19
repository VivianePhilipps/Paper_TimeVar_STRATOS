################################################################################
#                                                                              #
#                                   simJM                                      #
#                                                                              #
################################################################################

## Generates a dataset according to a joint shared random effect model

## Generates longitudinal data with fixed effects for Gaussian, curvilinear or ordinal outcomes and an event linked to the longitudinal part through the random effets or the current value of the latent process.

## Parameters :
## beta : the coefficients for the fixed effects (intercept firtst, then time, the binary covariate, the interaction between time and the binary covariate, and the continuous covariate)
## t.ecart : the gap between two visits
## t.marge : the variability in the visits spacing (simulated uniformly between -abs(t.marge) and abs(t.marge))
## t.max : the maximum measurement time
## Xbin : the proportion of success for the binary (Bernoulli) covariate, named X is the returned dataset
## Xcont : the mean and sd for the continuous (Gaussian) covariate , named Xc is the returned dataset
## n : the number of subjects
## B.random : the upper part of the random effects variance matrix
## nRE : the number of random effects
## dropout : a vector specifying the proportion of dropout and the miminum time of follow-up (time of dropout is simulated uniformly between dropout[2]and t.max)
## pMCAR : proportion of missing values with a completely at random mechanism
## pMAR : proportion of missing values with a missing at random mechanism
## sigma : the sd of the (Gaussian) error measurement
## seuils : a  named list specifying the thresholds for ordinal outcomes, the splines parameters for curvilinear outcomes or the linear parameters for Gaussian outcomes. The names shoukd begin with "thresh" for ordinal outcomes, "spline" for curvilinear outcomes and "linear" for Gaussian outcomes.
## modalites : a list specifying the levels for ordinal outcomes or the splines nodes for a curvilinear outcomes. Should be NULL for Gaussian outcomes.
## nomsY : the name of the outcomes in the return dataset. By default, the outcomes are named Y1, Y2, etc.
## nEvent : number of competing events. Default to 1.
## weibull : the Weibull parameters for the baseline hazard of the simulated time-to-event.
## piecewise : a list containing the parameters and nodes for the piecewise constant baseline hazard of the simulated time-to-event.
## Xsurv : the proportion of success for the binary (Bernoulli) covariate included in the proportional hazard model. To use the same covariate as in the longitudinal part, use Xsurv="Xbin" or Xsurv="Xcont".
## betaSurv : the effect of the covariate in the survival model.
## sharedtype : "re" for an association through the random effects (default), "value" for an association through the current value of the longitudinal process.
## association : the parameter(s) associated to the random effects (or current value) in the survival model.
## seed : the seed
## entry : for left truncated data, an expression to simulated entry time
## model : a list of formula for the fixed and random parts. By default, fixed=~t*X+Xc and random=~t.
## NAafterEvent : should the longitudinal outcome(s) be missing at measurement times following the time-to-event? Default to TRUE.

## Returns :
## A sample of n subjects including the following columns :
##  - ID : the subjects identifiers
##  - t : the measurement time
##  - platent : the fixed and random parts of the longitudinal model (Xbeta + Zu)
##  - X :  the binary covariate
##  - Xc : the continuous covariate
##  - Y1, Y2, etc : the longitudinal outcomes
##  - T0 : the entry time
##  - Xs : the covariate for the survival model
##  - T : the time-to-event
##  - D :  the indicator of event

## Depends : uniroot, rmvnorm

## Author : Viviane Philipps, Tiphaine Saulnier, Anais Rouanet

## Example :
##     simJM(beta=c(0,-1,0.5,0), # the fixed effects
##            t.ecart=2, t.marge=0.1, t.max=10, # 10 years follow-up with visits approx. every 2 year (between 1.9 years and 2.1 years for the first visit after baseline, between 3.9 and 4.1 for the second, ...)
##            Xbin=0.5, # binary covariate with proportion=50%
##            n=100, # 100 subjects
##            B.random=c(1,0,0.8), # variance of 1 and 0.8 for the random intercept and random effect of time respectively
##            nRE=2, # two random effects (intercept and time)
##            sigma=1, # the measurement error has variance 1
##            seuils=list(linearY1=c(0,1)), # we simulated only one outcome that is Gaussian
##            weibull=c(0.015,2), # the Weibull baseline hazard
##            Xsurv="Xbin", betaSurv=0.3, # the binary covariate is also included in the survival part with an effect of 0.3
##            association=c(-0.2,-0.4)) # the random intercept has an effect of -0.2 in the survial model, the random effect of time has an effect of -0.4.

simJM <- function(beta,t.ecart,t.marge,t.max,Xbin=NULL,Xcont=NULL,n,B.random,nRE,dropout=NULL,pMCAR,pMAR,sigma,seuils,modalites,nomsY,
                   nEvent=1,weibull=NULL,piecewise=NULL,Xsurv=0,betaSurv=rep(0,nEvent),sharedtype="re",association=rep(0,nEvent*nRE),
                   seed, entry=NULL, model=NULL,NAafterEvent=TRUE)
{
    
####   simulation des donnees    ####

    if(missing(seed)) seed <- round(abs(rnorm(1,mean=5)),5)*100000
    print(paste("seed = ",seed,sep=""))
    set.seed(seed)
    options(warn=-1)
    on.exit(options(warn=0))
    
    if(missing(nomsY))
    {
        nomsY <- paste("Y",1:length(sigma),sep="")
    }
    ny <- length(nomsY) 

    ## si model est un modele estime, prendre les formules
    if(inherits(model, c("hlme","lcmm","multlcmm","jointLPM","Jointlcmm")))
    {
        modele <- list(fixed=formula(paste("~",model$call$fixed[3])),
                       random=formula(paste("~",model$call$random[2])),
                       tsurvdizaine=FALSE)
    }
    else
    {
        modele <- model
        if(is.null(modele$tsurvdizaine)) modele$tsurvdizaine <- FALSE
    }
    
    ## temps inverse si t.max negatif
    t.inv <- FALSE
    if(t.max<0) t.inv <- TRUE
    t.max <- abs(t.max)
    t.ecart <- abs(t.ecart)

    
    ## simuler les donnees pour n sujets
    res <- NULL
    i <- 1
    while(i<=n)
    {

        ## simuler les temps, selon une loi uniforme et avec ecart moyen entre chq visite
        t <- 0
        nbt <- 1
        while(t[nbt]<t.max)
        {
            t <- c(t,t.ecart*nbt+runif(1,-abs(t.marge),abs(t.marge)))
            nbt <- nbt+1
        }
        if(t[nbt]>t.max) t <- t[-nbt]
        
        ## censure
        tcensure <- t.max
        if(!is.null(dropout))
        {
            idcensure <- rbinom(n=1,size=1,prob=dropout[1])
            if(idcensure==1)
            {
                tcensure <- runif(1,abs(dropout[2]),t.max)
                if(any(t>tcensure)) t <- t[-which(t>tcensure)]
            }
        }

        ## inverser le temps
        if(t.inv==TRUE) t <- -t


        ## simuler le temps d'entree et l'ajouter aux temps simules
        age0 <- 0
        if(!is.null(entry))
        {
            age0 <- eval(match.call()$entry)
        }
        t <- t + age0


        ## simuler les effets aleatoires
        varRE <- matrix(0,nrow=nRE,ncol=nRE)
        if(length(B.random)==nRE)
        {
            diag(varRE) <- B.random
        }
        else
        {
            varRE[upper.tri(varRE,diag=TRUE)] <- B.random
            varRE <- t(varRE)
            varRE[upper.tri(varRE,diag=TRUE)] <- B.random
        }
        brandom <- as.numeric(mvtnorm::rmvnorm(1, mean=rep(0,nRE), sigma=varRE))

        
        ##initialisation du data frame
        donnees <- data.frame(1,t=t)

        
        ## simuler la variable expl
        if(length(Xbin))
        {
            X <- rbinom(1,size=1,prob=Xbin[1])
            
            ## ajouter X dans donnees
            donnees <- cbind(donnees,X)
            
            ## ajouter interaction
            donnees <- cbind(donnees,donnees[,2]*donnees[,3])
        }


        Xc <- NULL
        if(length(Xcont))
        {
            Xc <- rnorm(1,mean=Xcont[1],sd=Xcont[2])
            
            ## ajouter X dans donnees
            donnees <- cbind(donnees,Xc)
        }
        
        ## survie
        ##TS: deplacement partie survie ici car impacte brandom dans le cas du niveau courant partage
        if(!is.null(weibull) | !is.null(piecewise))
        {
            if(Xsurv=="Xbin")
            {
                xs <- X
            }
            else
            {
                if(Xsurv=="Xcont")
                {
                    xs <- Xc
                }
                else
                {
                    xs <- rbinom(1,size=1,prob=Xsurv)
                }
            }

            nweib <- length(weibull)/nEvent
            
            Tevt <- rep(NA,nEvent)
            Devt <- rep(NA,nEvent)
            for (j in 1:nEvent){
                
                expSurv <- exp(xs*betaSurv[j])
                
                if(sharedtype == "re"){    # shared random effects
                    
                    if(!is.null(weibull))
                    {
                        weibull_j <- weibull[(j-1)*nweib+1:nweib]
                        if(nweib==2)
                        {
                            ## attention prm doivent etre tels que H(t) = w1 * t^w2 (ie parametrisation de logscale=TRUE dans lcmm)
                            tempsSurv <- rweibull(1, shape=weibull_j[2], scale=(weibull_j[1]*expSurv*exp(t(association[(j-1)*nRE+1:nRE])%*%brandom))^(-1/weibull_j[2]))
                            #cat("tempsSurv=",tempsSurv,"\n")
                        }
                        else
                        {
                            ## ici prm tels que H(t) = (w1*t)^w2
                            unif <- runif(1)
                            tempsSurv <- ((-log(unif)/(expSurv*exp(t(association[(j-1)*nRE+1:nRE])%*%brandom)))^(1/weibull_j[2]))/weibull_j[1] + weibull_j[3]
                        }
                    }
                    
                    if(!is.null(piecewise))
                    {
                        zl <- piecewise[[paste("nodes",j,sep="")]]
                        bl <- piecewise[[paste("brisq",j,sep="")]]
                                                
                        surv <- function(t,zl,bl,expS,p)
                        {
                            j <- which.max(zl[which(zl<=t)])
                            if(j==1) som <- 0 
                            else som <- sum(bl[1:(j-1)]*(zl[2:j]-zl[1:(j-1)]))
                            
                            if(j<length(zl)) surv <- exp(-(som+bl[j]*(t-zl[j]))*expS)
                            else surv <- exp(-som*expS)
                            
                            return(surv-p)
                        }
                        
                        unif <- runif(1)
                        zero <- try(uniroot(surv,interval=c(zl[1],zl[length(zl)]),
                                            zl=zl,bl=bl,
                                            expS=expSurv*exp(t(association[(j-1)*nRE+1:nRE])%*%brandom),
                                            p=unif),
                                    silent=TRUE)
                        if(class(zero)=="try-error") tempsSurv <- Inf
                        else tempsSurv <- zero$root
                        
                    }
                }
                else if(sharedtype == "value"){    # shared latent process current level    #TS
                    # code inspired from simulateJM() fct from package JM simulating data with survival model adjusted on the unique outcome current level
                    # principle : compute cumulative risk fct values according to several times of event and stop when equals a randomly generated uniform value 
                    
                    if(!is.null(piecewise))
                        stop("piecewise with sharedtype=value not programmed yet")
                    
                    if(!is.null(weibull)){
                        if(nweib != 2) stop("3 parameters Weibull risk not programmed yet")
                        weibull_j <- weibull[(j-1)*nweib+1:nweib]
                        
                        u <- runif(1)  #TS: sim 1 valeur U[0,1] par sujet
                        
                    
                        ###################################################################################
                        # function to compute the inverse survival / cumulative risk function
                        invS <- function (t, u) {
                            
                            TD <- function (v) {
                                # function to compute the (time-dependent) latent process current level for patient i at time v, multiplicated by the association prm
                                dd <- data.frame(int = rep(1, length(v)))
                                dd$t <- pmax(v, 0)
                                if(length(Xbin))
                                {
                                    dd$X <- rep(X, length(t))
                                    dd$interactionXt <- dd$X*dd$t
                                }
                                if(!is.null(Xc)) dd$Xc <- Xc
                                
                                if(is.null(modele$fixed))
                                {
                                    ## par defaut, fixed=~t*X+Xc
                                    XXbeta <- apply(dd,1,function(x) sum(x*beta))
                                }else{
                                    XX <- model.matrix(object=modele$fixed, data=dd) 
                                    XXbeta <- XX %*% beta
                                }
                                    
                                if(is.null(modele$random))
                                {
                                    XXrandom <- matrix(1,nrow=nrow(dd),ncol=1)
                                    if(nRE==2) XXrandom <- cbind(1,dd[,2])
                                    
                                    ZZu <- XXrandom %*% brandom
                                }else{
                                    ZZ <- model.matrix(modele$random, data=dd) 
                                    ZZu <- ZZ %*% brandom
                                }
                                                                    
                                out <- (XXbeta + ZZu) * association[j]  # multiplied by the association prm    
                                if(any(is.na(out))){
                                    message("Pb Ypredict = NA \n")
                                    browser()
                                }
                                
                                return(out)
                            }
                            
                            h <- function (s) {
                                # instant risk fct for patient i at time s
                                TD.i <- TD(s)
                                
                                    #print("Make sure weibull[1] = lambda^alpha and weibull[2] = alpha")
                                instant <- weibull_j[1] * weibull_j[2] * s^(weibull_j[2]-1) * exp(TD.i) * expSurv
                                
                                return(instant)
                            }
                            return(integrate(h, lower = 0, upper = t)$value + log(u))  # cumulative risk fct value + log(uniform value)
                        }
                        ###################################################################################

                        risqcum_min <- invS(age0, u=u) 
                        risqcum_max <- invS(age0+t.max, u=u)
                        if(risqcum_min > 0)
                        {
                            tempsSurv <- age0-1 # event before entry
                        }
                        else
                        {
                            if(risqcum_max < 0)
                            {
                                tempsSurv <- age0+t.max+1 # event after the end of the study
                            }
                            else
                            {
                                ## find event time between 0 and t.max
                                tempsSurv <- uniroot(invS, interval = age0+c(1e-05, t.max), u = u)$root
                            }
                        }
                       
                        
                        ## temp=1000
                        ## while(inherits(Root, "try-error") & (temp<1000)) { 
                        ##     cat("resampling random effects indiv: ",i, " ",temp,"\n")
                        ##     cat("brandom =",brandom,"\n")
                        ##     cat("log(u) =",log(u),"\n")
                        ##     inth <- rep(NA,11)
                        ##     for(k in 0:10) inth[k] <- invS(k,u)
                        ##     x11();plot(0:10, inth, type="l",ylim=c(-3,1))
                        ##     ##TS: si pas de tps trouve precedemment, re-essai en retirant aleatoirement les EAs du sujet, et continue jusqu a ce qu il trouve un tps correspondant
                        ##     brandom <- as.numeric(mvtnorm::rmvnorm(1, mean=rep(0,nRE), sigma=varRE))
                            
                        ##     Root <- try(expr = uniroot(invS, 
                        ##                                interval = c(1e-05, 10*t.max), 
                        ##                                u = u)$root, 
                        ##                 silent = TRUE)
                        ##     temp = temp +1
                        ## }
                        ## if(inherits(Root, "try-error"))
                        ## {
                        ##     tempsSurv <- as.numeric(2*t.max) #stop("unable to generate a survival time")
                        ## }
                        ## else
                        ## {
                        ##     tempsSurv <- Root  
                        ## }
                    }
                }

                if(tempsSurv < age0) break #on sort de la boucle j

                Tevt[j] <- tempsSurv                
            
            }#fin boucle j (1:nEvent)
            
            if(tempsSurv < age0) next # on passe au sujet suivant (boucle i)
           
            if(modele$tsurvdizaine==FALSE)
            {
                T <- min(Tevt,age0+tcensure, max(t))
                whichD <- which.min(c(Tevt,age0+tcensure, max(t)))
            }
            else
            {
                T <- min(Tevt,(age0+tcensure)/10, max(t)/10)
                whichD <- which.min(c(Tevt,(age0+tcensure)/10))
            }

            if(whichD > nEvent){
                D <- 0
            }else{
                D <- whichD
            }
        } ## fin survie
        
        ## longitudinale
        ## calculer la valeur du processus latent
        ##fonction pour passer de H(Y) a Y
        transfinv <- function(ytilde,k)
        {
            
            nam <- names(seuils)
            nam <- tolower(nam)
            namk <- substr(nam[k],1,6)
            linktype <- pmatch(namk,c("linear","spline","thresh"))
            
            yy <- NA
            
            if(linktype==3) ## thresholds
            {
                v <- c(ytilde,lapply(seuils[[k]], function(x){c(x[1], x[1]+cumsum(x[-1]^2))}))
                indic <- c(1,rep(0,length(seuils[[k]])))
                indic <- indic[order(v)]
                pos <- which(indic==1)
                yy <- modalites[[k]][pos]
            }
            
            if(linktype==1) ##linear
            {
                yy <- seuils[[k]][2]*ytilde+seuils[[k]][1]
            }
            
            if(linktype==2) ##splines
            {
                z <- modalites[[k]] #spline nodes
                ff <- function(x,hy,z,b){transfo_spl(x,z,b)-hy}
                yy <- uniroot(f=ff,lower=z[1],upper=z[length(z)],
                              hy=ytilde,z=z,b=seuils[[k]])$root
            }
            return(yy)
        }


        if(is.null(modele$fixed))
        {
            ## par defaut, fixed=~t*X+Xc
            Xbeta <- apply(donnees,1,function(x) sum(x*beta))
        }
        else
        {
            Xb <- model.matrix(object=modele$fixed, data=donnees)
            Xbeta <- Xb %*% beta
        }
        
        if(is.null(modele$random))
        {
            Xrandom <- matrix(1,nrow=nrow(donnees),ncol=1)
            if(nRE==2) Xrandom <- cbind(1,donnees[,2])
            
            Zu <- Xrandom %*% brandom
        }
        else
        {
            Zu <- model.matrix(modele$random, data=donnees) %*% brandom
        }
        platent <- Xbeta + Zu
        
        if(any(is.na(platent)))
            message("Pb platent = NA  \n")
        
        ## simuler les erreurs de mesures
        erreurs <- matrix(sapply(sigma,rnorm,n=length(t),mean=0),nrow=length(t),ncol=ny)
            
        ## calculer Lambda+erreurs
        ytransf <- matrix(sweep(erreurs,STATS=platent,FUN="+",MARGIN=1),nrow=length(t),ncol=ny)
            
        ## prendre l'inverse pour avoir des simulations de Y
        y <- NULL
        for(k in 1:ny)
        {
            yk <- sapply(ytransf[,k],transfinv, k = k) 
            y <- c(y,yk)
        } 
         
        ## garder t,X,y,age0
        d <- data.frame(rep(i,length(t)),t,platent)

        
        if((!length(Xbin)) & (!length(Xcont)))
        {
            d <- cbind(d,matrix(y,nrow=length(t),ncol=ny),age0)
            colnames(d) <- c("ID","t","platent",nomsY,"T0")
        }
        if((length(Xbin)) & (!length(Xcont)))
        {
            d <- cbind(d,X,matrix(y,nrow=length(t),ncol=ny),age0)
            colnames(d) <- c("ID","t","platent","X",nomsY,"T0")
        }
        if((!length(Xbin)) & (length(Xcont)))
        {
            d <- cbind(d,Xc,matrix(y,nrow=length(t),ncol=ny),age0)
            colnames(d) <- c("ID","t","platent","Xcont",nomsY,"T0")
        }
        if((length(Xbin)) & (length(Xcont)))
        {
            d <- cbind(d,X,Xc,matrix(y,nrow=length(t),ncol=ny),age0)
            colnames(d) <- c("ID","t","platent","X","Xcont",nomsY,"T0")
        }
        
        ## combiner donnees longitudinales et de survie ds un meme dataframe
        if(!is.null(weibull) || !is.null(piecewise)){
            oldnames <- colnames(d)
            d <- cbind(d,xs,T,D)
            colnames(d) <- c(oldnames,"Xs","T","D")
            
            if(NAafterEvent)
            {
                ## mettre Y=NA si t > T
                jj <- which(d$t > (d$T * (1 + 9*as.numeric(modele$tsurvdizaine))))
                if(length(jj)) d[jj,nomsY] <- NA
            }
        }
        
        i <- i+1
        res <- rbind(res,d)
        
    }
    

### introduction de donnees manquantes
    ## introduire les NA en MCAR
    if(!missing(pMCAR))
    {
        ##nbmestot <- length(unlist(res[,nomsY])) ## NA can be everywhere
        nbmestot <- length(unlist(res[which(res[,"t"]>0),nomsY])) ## NA only on t>0
        isna <- rbinom(n=nbmestot,size=1,prob=pMCAR)
        idMCAR <- which(isna==1)
        y <- unlist(res[which(res[,"t"]>0),nomsY])
        y[idMCAR] <- NA
        y <- matrix(y,ncol=length(nomsY))
        data <- res
        data[which(res[,"t"]>0),nomsY] <- y
        
    }
    else
    {
        data <- res
    }
    
    ## on ne peut pas avoir du MCAR et du MAR
    ## si les 2 sont specifies, on n'aura que du MAR au final
    
    ## introduire les NA en MAR
    if(!missing(pMAR))
    {
        q <- quantile(res[,"platent"],probs=0.25,names=FALSE)
        nbmestot <- length(unlist(res[which(res[,"platent"]<=q),nomsY]))
        isna <- rbinom(n=nbmestot,size=1,prob=4*pMAR)
        idMAR <- which(isna==1)
        y <- unlist(res[which(res[,"platent"]<=q),nomsY])
        y[idMAR] <- NA
        y <- matrix(y,ncol=length(nomsY))
        data <- res
        data[which(res[,"platent"]<=q),nomsY] <- y
    }
    
    return(data)
}

