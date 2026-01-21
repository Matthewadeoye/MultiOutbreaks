#Inference Method 1 -- Smoothing
SMOOTHING_INFERENCE<- function(y, e_it, Modeltype, adjmat, step_sizes, num_iteration = 15000, sdBs=0.03, sdGs=0.05, sdLambdas=0.03, sdCop=0.1, y_total=NULL){
  start_time <- Sys.time()
  ndept <- nrow(e_it)
  time <- ncol(e_it)
  nstrain<- dim(y)[3]
  nstate<- 2^nstrain
  Bits<- encodeBits(nstrain)
  gh <- statmod::gauss.quad(30, kind = "hermite")

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  SumYk_vec<- numeric(nstrain)
  SumYk_vec[1]<- sum(y[,,1], na.rm = T)
  sumY<- ifelse(is.na(y[,,1]),0,y[,,1])
  for(k in 2:nstrain){
    sumY<- sumY + ifelse(is.na(y[,,k]),0,y[,,k])
    SumYk_vec[k]<- sum(y[,,k], na.rm = T)
  }

  if(is.null(y_total)){
    sumY1<- y[,,1]
    for(k in 2:nstrain){
      sumY1<- sumY1 + y[,,1]
    }
    y_total<- sumY1
  }else{
    y_total<- y_total
  }

  crudeResults<- crudeEst(sumY, e_it)
  crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)

  initGs<- gtools::rdirichlet(nstate, rep(1, nstate))
  initstateD<- stationarydistArma_cpp(initGs)[ncol(initGs)]

  Model<- ifelse(Modeltype>0,1,0)
  if(Modeltype %in% c(0,1,2,7)){
    n_copParams<- 0
    n_factloadings<- 0
  }else if(Modeltype %in% c(3,4)){
    n_copParams<- (nstrain*(nstrain-1))/2
    n_factloadings<- nstrain
  }else if(Modeltype %in% c(5,6)){
    n_copParams<- 1
    n_factloadings<- 0
  }

  if(Modeltype %in% c(0,1,3,5)){
    num_Gammas<- 2
    MC_chain<- matrix(NA, nrow=num_iteration, ncol=num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+n_copParams)
    MC_chain[1,]<- c(runif(num_Gammas), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_factloadings), rep(0.1, n_copParams))
  }else if(Modeltype %in% c(2,4,6)){
    num_Gammas<- 2 * nstrain
    MC_chain<- matrix(NA, nrow=num_iteration, ncol=num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+n_copParams)
    MC_chain[1,]<- c(runif(num_Gammas), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_factloadings), rep(0.1, n_copParams))
  }else if(Modeltype==7){
    num_Gammas<- nstate * nstate
    MC_chain<- matrix(NA, nrow=num_iteration, ncol=num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+n_copParams)
    MC_chain[1,]<- c(as.numeric(t(initGs)), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_factloadings), rep(0.1, n_copParams))
  }

  Q_r<- MC_chain[1,num_Gammas + 1] * RW2PrecMat
  Q_s<- MC_chain[1,num_Gammas + 2] * RW1PrecMat
  Q_u<- MC_chain[1,num_Gammas + 3] * R

  #Compute gradients
  if(Modeltype %in% c(1,2,5,6)){
    JointTPM<- Multipurpose_JointTransitionMatrix(MC_chain[1,1:num_Gammas], nstrain, MC_chain[1,num_Gammas+3+time+12+ndept+nstrain+nstrain], Modeltype)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
  }else if(Modeltype %in% c(3,4)){
    JointTPM<- Multipurpose_JointTransitionMatrix2(MC_chain[1,1:num_Gammas], nstrain, MC_chain[1,num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)], Modeltype, gh)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
    if(any(!is.finite(JointTPM))) JointTPM<- initGs
  }else if(Modeltype == 7){
    JointTPM<- initGs
  }

  if(Model == 0) JointTPM<- matrix(0, nstate, nstate)

  Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[1, num_Gammas+3+(1:time)], s=MC_chain[1, num_Gammas+3+time+(1:12)], u=MC_chain[1, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total = y_total)
  likelihoodcurrent<- Allquantities$loglike
  priorcurrentRcomps<- randomwalk2_cpp(MC_chain[1, num_Gammas+3+(1:time)], MC_chain[1, num_Gammas+1])
  priorcurrentScomps<- seasonalComp2_cpp(MC_chain[1, num_Gammas+3+time+(1:12)], MC_chain[1, num_Gammas+2], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1_cpp(MC_chain[1, num_Gammas+3+time+12+(1:ndept)], MC_chain[1, num_Gammas+3], R, rankdef)

  grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

  deltaP<- 1
  RMdelta<- 1/(0.234*(1-0.234))

  propShape<- 0.01+SumYk_vec

  for (i in 2:num_iteration) {

    MC_chain[i,num_Gammas+1]<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, num_Gammas+3+(1:time)]) %*% RW2PrecMat %*% MC_chain[i-1, num_Gammas+3+(1:time)])/2)
    MC_chain[i,num_Gammas+2]<- rgamma(1, shape = 1 + 11/2, rate = 0.001 + (t(MC_chain[i-1, num_Gammas+3+time+(1:12)]) %*% RW1PrecMat %*% MC_chain[i-1, num_Gammas+3+time+(1:12)])/2)
    MC_chain[i,num_Gammas+3]<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)])/2)

    Q_r<- MC_chain[i,num_Gammas+1] * RW2PrecMat
    Q_s<- MC_chain[i,num_Gammas+2] * RW1PrecMat
    Q_u<- MC_chain[i,num_Gammas+3] * R

    current_r <- MC_chain[i-1, num_Gammas+3+(1:time)]
    current_s <- MC_chain[i-1, num_Gammas+3+time+(1:12)]
    current_u <- MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)]

    #Update s
    Mmatcs<- as.numeric(grad_current$cov_s %*% grad_current$grad_s)

    eps_s <- rnorm(12)
    proposedScomps <- as.numeric(MC_chain[i-1, num_Gammas+3+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs + step_sizes$s * chol(grad_current$cov_s) %*% eps_s)
    proposedScomps<- proposedScomps - mean(proposedScomps)

    Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=current_r, s=proposedScomps, u=current_u, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total = y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- Allquantities$loglike
    q_prop <- mvnfast::dmvn(proposedScomps, mu = MC_chain[i-1, num_Gammas+3+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs, sigma = grad_current$cov_s * step_sizes$s^2, log = TRUE)
    q_curr <- mvnfast::dmvn(MC_chain[i-1, num_Gammas+3+time+(1:12)], mu = proposedScomps + 0.5 * step_sizes$s^2 * Mmatps, sigma = grad_proposed$cov_s * step_sizes$s^2, log = TRUE)
    priorproposedScomps<- seasonalComp2_cpp(proposedScomps, MC_chain[i, num_Gammas+2], RW1PrecMat)

    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop
    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
      MC_chain[i, num_Gammas+3+time+(1:12)]<- proposedScomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentScomps<- priorproposedScomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, num_Gammas+3+time+(1:12)]<- MC_chain[i-1, num_Gammas+3+time+(1:12)]
    }

    #Update r
    eps_r <- rnorm(time)
    Mmatrc<- as.numeric(grad_current$cov_r %*% grad_current$grad_r)

    proposedRcomps <- as.numeric(current_r + 0.5 * step_sizes$r^2 * Mmatrc + step_sizes$r * chol(grad_current$cov_r) %*% eps_r)
    proposedRcomps<- proposedRcomps - mean(proposedRcomps)

    Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=proposedRcomps, s=MC_chain[i, num_Gammas+3+time+(1:12)], u=current_u, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total = y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcomps, mu = current_r + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(current_r, mu = proposedRcomps + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike
    priorproposedRcomps <- randomwalk2_cpp(proposedRcomps, MC_chain[i, num_Gammas+1])

    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop

    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
      MC_chain[i, num_Gammas+3+(1:time)] <- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentRcomps<- priorproposedRcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, num_Gammas+3+(1:time)]<- MC_chain[i-1, num_Gammas+3+(1:time)]
    }

    #Update u
    eps_u <- rnorm(ndept)

    proposedUcomps <- as.numeric(MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u + step_sizes$u * eps_u)
    proposedUcomps<- proposedUcomps - mean(proposedUcomps)

    #Pseudo-Gibbs A_k's proposal
    propRate<- as.numeric(Allquantities$poisMean4GibbsUpdate) + 0.01/exp(-15)
    proposedAks<- log(rgamma(nstrain, shape = propShape, rate = propRate))

    Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=proposedUcomps, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=proposedAks, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total = y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    likelihoodproposed<- Allquantities$loglike

    q_prop <- sum(dnorm(proposedUcomps, mean = MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u, sd = step_sizes$u, log = TRUE))
    q_curr <- sum(dnorm(MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)], mean = proposedUcomps + 0.5 * step_sizes$u^2 * grad_proposed$grad_u, sd = step_sizes$u, log = TRUE))

    priorproposedUcomps<- logIGMRF1_cpp(proposedUcomps, MC_chain[i, num_Gammas+3], R, rankdef)
    priorcurrentUcomps<- logIGMRF1_cpp(MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)], MC_chain[i, num_Gammas+3], R, rankdef)

    proposalcurrentAks <- sum(dgamma(exp(MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]),
                                     shape=propShape, rate=propRate, log=TRUE) + MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]) #Jacobian
    proposalproposedAks <- sum(dgamma(exp(proposedAks),
                                      shape=propShape, rate=propRate, log=TRUE) + proposedAks)  #Jacobian

    priorcurrentAks <- sum(dgamma(exp(MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]),
                                  shape=0.01, rate=0.01/exp(-15), log=TRUE) + MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])
    priorproposedAks <- sum(dgamma(exp(proposedAks),
                                   shape=0.01, rate=0.01/exp(-15), log=TRUE) +  proposedAks)

    log_alpha_u_a <- (likelihoodproposed + priorproposedUcomps + q_curr + priorproposedAks + proposalcurrentAks
                      - likelihoodcurrent - priorcurrentUcomps - q_prop - priorcurrentAks - proposalproposedAks)

    if (is.finite(log_alpha_u_a) && log(runif(1)) < log_alpha_u_a){
      MC_chain[i, num_Gammas+3+time+12+(1:ndept)]<- proposedUcomps
      MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- proposedAks
      likelihoodcurrent<- likelihoodproposed
      grad_current<- grad_proposed
    }else{
      MC_chain[i, num_Gammas+3+time+12+(1:ndept)]<- MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)]
      MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
    }

    if(Model == 0){
      MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)]
      MC_chain[i, 1:2]<- MC_chain[i-1, 1:2]
    }else{
      #proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain)))

      proposedB <- rnorm(nstrain, mean = MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain))
      if(any(proposedB < 0)){
        MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)]
      }else{
      priorcurrentB<- sum(dgamma(MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

      Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=proposedB, Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total = y_total)
      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedB
                     - likelihoodcurrent - priorcurrentB)

      #print(paste("mh.ratioB = ", mh.ratio))

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- proposedB
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
      }
      else{
        MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)]
      }
      sdBs<- sdBs * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))
    }
      if(Modeltype %in% c(1,2)){
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))

        proposedLambdas<- rep(0, nstrain)
        JointTPM1<- Multipurpose_JointTransitionMatrix(proposedGs, nstrain, proposedLambdas, Modeltype)
        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total = y_total)
        grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

        likelihoodproposed<- Allquantities$loglike

        mh.ratio<- exp(likelihoodproposed + priorproposedGs
                       - likelihoodcurrent - priorcurrentGs)

        #print(mh.ratio)

        if(!is.na(mh.ratio) && runif(1) < mh.ratio){
          MC_chain[i, 1:num_Gammas]<- proposedGs
          likelihoodcurrent<- likelihoodproposed
          grad_current<- grad_proposed
          JointTPM<- JointTPM1
        }
        else{
          MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
        }
        sdGs<- sdGs * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))
      }else if(Modeltype %in% c(3,4)){
        #Transition probabilities update
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))

        JointTPM1<- Multipurpose_JointTransitionMatrix2(proposedGs, nstrain, MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)], Modeltype, gh)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        if(any(!is.finite(JointTPM1))){
          MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
        }else{

          Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total = y_total)
          grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

          likelihoodproposed<- Allquantities$loglike

          mh.ratioGC<- exp(likelihoodproposed + priorproposedGs
                           - likelihoodcurrent - priorcurrentGs)

          #print(mh.ratioGC)

          if(!is.na(mh.ratioGC) && runif(1) < mh.ratioGC){
            MC_chain[i, 1:num_Gammas]<- proposedGs
            likelihoodcurrent<- likelihoodproposed
            grad_current<- grad_proposed
            JointTPM<- JointTPM1
          }
          else{
            MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
          }
          sdGs<- sdGs * exp((RMdelta/i) * (min(mh.ratioGC, 1) - 0.234))
        }

        #Factor loadings update
        proposedLambdas<- rnorm(nstrain,mean=MC_chain[i-1,num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)], sd=rep(sdLambdas, nstrain))
        #proposedLambdas<- ifelse(proposedLambdas<1, proposedLambdas, 2-proposedLambdas)

        priorcurrentLambdas<- sum(dunif(MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)], min = rep(-1, n_factloadings), max = rep(1, n_factloadings), log=TRUE))
        priorproposedLambdas<- sum(dunif(proposedLambdas, min = rep(-1, n_factloadings), max = rep(1, n_factloadings), log=TRUE))

        JointTPM1<- Multipurpose_JointTransitionMatrix2(MC_chain[i, 1:num_Gammas], nstrain, proposedLambdas, Modeltype, gh)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        if(any(!is.finite(JointTPM1)) || any(abs(proposedLambdas) > 1)){
          MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]
        }else{

          Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total = y_total)
          grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

          likelihoodproposed<- Allquantities$loglike

          mh.ratioGC<- exp(likelihoodproposed + priorproposedLambdas
                           - likelihoodcurrent - priorcurrentLambdas)

          #print(mh.ratioGC)

          if(!is.na(mh.ratioGC) && runif(1) < mh.ratioGC){
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- proposedLambdas
            likelihoodcurrent<- likelihoodproposed
            grad_current<- grad_proposed
            JointTPM<- JointTPM1
          }
          else{
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]
          }
          sdLambdas<- sdLambdas * exp((RMdelta/i) * (min(mh.ratioGC, 1) - 0.234))
        }
      }else if(Modeltype %in% c(5,6)){
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))

        if(nstrain==2){
          proposedcopPs<- rnorm(1,mean=MC_chain[i-1, ncol(MC_chain)], sd=sdCop)
          JointTPM1<- Multipurpose_JointTransitionMatrix(proposedGs, nstrain, proposedcopPs, Modeltype)
        }else{
          proposedcopPs<- rnorm(1,mean=log(MC_chain[i-1, ncol(MC_chain)]), sd=sdCop)
          proposalcurrentcop<- dnorm(log(MC_chain[i-1, ncol(MC_chain)]), mean=proposedcopPs, sd=sdCop, log = T) + MC_chain[i-1, ncol(MC_chain)]
          proposalproposedcop<- dnorm(proposedcopPs, mean=log(MC_chain[i-1, ncol(MC_chain)]), sd=sdCop, log = T) + exp(proposedcopPs)
          JointTPM1<- Multipurpose_JointTransitionMatrix(proposedGs, nstrain, exp(proposedcopPs), Modeltype)
        }

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)
        if(any(!is.finite(JointTPM1))){
          MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_copParams)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_copParams)]
          MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
        }else{
          JointTPM1 <- JointTPM1 / rowSums(JointTPM1)

          Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,  y_total = y_total)
          grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

          likelihoodproposed<- Allquantities$loglike

          if(nstrain==2){
            mh.ratio<- exp(likelihoodproposed + priorproposedGs
                           - likelihoodcurrent - priorcurrentGs)
          }else{
            mh.ratio<- exp(likelihoodproposed + priorproposedGs + proposalcurrentcop
                           - likelihoodcurrent - priorcurrentGs - proposalproposedcop)
          }

          #print(mh.ratio)

          if(!is.na(mh.ratio) && runif(1) < mh.ratio){
            MC_chain[i, 1:num_Gammas]<- proposedGs
            likelihoodcurrent<- likelihoodproposed
            grad_current<- grad_proposed
            JointTPM<- JointTPM1
            MC_chain[i, ncol(MC_chain)]<- ifelse(nstrain==2,proposedcopPs,exp(proposedcopPs))
          }
          else{
            MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
            MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
          }
          sdGs<- sdGs * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))
          sdCop<- sdCop * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))
        }
      }else if(Modeltype==7){

        for(n in 1:nstate){
          index<- nstate * (n-1) + 1

          JointTPM[n, ] <- gtools::rdirichlet(1, rep(1, nstate) + deltaP * MC_chain[i-1, (index:(n*nstate))])

          proposalproposedGs<-  log(gtools::ddirichlet(JointTPM[n, ], MC_chain[i-1, (index:(n*nstate))]))
          proposalcurrentproposedGs<- log(gtools::ddirichlet(MC_chain[i-1, (index:(n*nstate))], JointTPM[n, ]))

          priorcurrentGs<- log(gtools::ddirichlet(MC_chain[i-1, (index:(n*nstate))], rep(1, nstate)))
          priorproposedGs<- log(gtools::ddirichlet(JointTPM[n, ], rep(1, nstate)))

          if(n == nstate){
            Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1, y_total = y_total)
            grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)
          }else{
            Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=0, y_total = y_total)
          }

          likelihoodproposed<- Allquantities$loglike

          mh.ratio<- exp(likelihoodproposed + priorproposedGs + proposalcurrentproposedGs
                         - likelihoodcurrent - priorcurrentGs - proposalproposedGs)

          #print(mh.ratio)

          if(!is.na(mh.ratio) && runif(1) < mh.ratio){
            MC_chain[i, (index:(n*nstate))]<- as.numeric(JointTPM[n, ])
            likelihoodcurrent<- likelihoodproposed
            grad_current<- grad_proposed
            deltaP<- max(0, deltaP-3)
          }
          else{
            MC_chain[i, (index:(n*nstate))]<- MC_chain[i-1, (index:(n*nstate))]
            deltaP<- deltaP + 1
          }
        }
      }
    }
    if(i %% 1000 == 0) cat("Iteration:", i, "\n")
  }

  if(Modeltype %in% c(0, 1)){
    colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep="")))
  }else if(Modeltype == 2){
    colnames(MC_chain) <- paste(c(paste0(rep(c("G12", "G21"), nstrain), "Strain", rep(1:nstrain,each=2)), "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep="")))
  }else if(Modeltype == 3){
    #Derive pair-correlations from factor-loadings
    copInd<- 1
    for(w in 1:(n_factloadings - 1)){
      for (z in (w + 1):n_factloadings){
        MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+copInd] <- MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+w] * MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+z]
        copInd <- copInd + 1
      }
    }
    colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), paste("FactorLoading", 1:nstrain, sep =""), paste("copulaParam", 1:n_copParams, sep="")))
  }else if(Modeltype == 4){
    #Derive pair-correlations from factor-loadings
    copInd<- 1
    for(w in 1:(n_factloadings - 1)){
      for (z in (w + 1):n_factloadings){
        MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+copInd] <- MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+w] * MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+z]
        copInd <- copInd + 1
      }
    }
    colnames(MC_chain) <- paste(c(paste0(rep(c("G12", "G21"), nstrain), "Strain", rep(1:nstrain,each=2)), "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), paste("FactorLoading", 1:nstrain, sep =""),paste("copulaParam", 1:n_copParams, sep="")))
  }else if(Modeltype == 5){
    colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), paste("copulaParam", 1:n_copParams, sep ="")))
  }else if(Modeltype == 6){
    colnames(MC_chain) <- paste(c(paste0(rep(c("G12", "G21"), nstrain), "Strain", rep(1:nstrain,each=2)), "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), paste("copulaParam", 1:n_copParams, sep ="")))
  }else if(Modeltype == 7){
    colnames(MC_chain) <- paste(c(paste0(rep("G_", nstate*nstate), rep(1:nstate, each=nstate), ",", 1:nstate), "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep="")))
  }
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}


#Inference Method 2 -- Sampling
FFBS_INFERENCE<- function(y, e_it, Modeltype, adjmat, step_sizes, num_iteration = 15000, sdBs=0.03, sdGs=0.05, sdLambdas=0.03, sdCop=0.1, y_total=NULL){
  start_time <- Sys.time()
  ndept <- nrow(e_it)
  time <- ncol(e_it)
  nstrain<- dim(y)[3]
  nstate<- 2^nstrain
  Bits<- encodeBits(nstrain)
  gh <- statmod::gauss.quad(30, kind = "hermite")

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  SumYk_vec<- numeric(nstrain)
  SumYk_vec[1]<- sum(y[,,1], na.rm = T)
  sumY<- ifelse(is.na(y[,,1]),0,y[,,1])
  for(k in 2:nstrain){
    sumY<- sumY + ifelse(is.na(y[,,k]),0,y[,,k])
    SumYk_vec[k]<- sum(y[,,k], na.rm = T)
  }

  if(is.null(y_total)){
    sumY1<- y[,,1]
    for(k in 2:nstrain){
      sumY1<- sumY1 + y[,,1]
    }
    y_total<- sumY1
  }else{
    y_total<- y_total
  }

  crudeResults<- crudeEst(sumY, e_it)
  crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)

  initGs<- gtools::rdirichlet(nstate, rep(1, nstate))
  initstateD<- stationarydistArma_cpp(initGs)[ncol(initGs)]

  Model<- ifelse(Modeltype>0,1,0)
  if(Modeltype %in% c(0,1,2,7)){
    n_copParams<- 0
    n_factloadings<- 0
  }else if(Modeltype %in% c(3,4)){
    n_copParams<- (nstrain*(nstrain-1))/2
    n_factloadings<- nstrain
  }else if(Modeltype %in% c(5,6)){
    n_copParams<- 1
    n_factloadings<- 0
  }

  if(Modeltype %in% c(0,1,3,5)){
    num_Gammas<- 2
    MC_chain<- matrix(NA, nrow=num_iteration, ncol=num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+n_copParams)
    MC_chain[1,]<- c(runif(num_Gammas), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0.1, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_factloadings), rep(0.1, n_copParams))
  }else if(Modeltype %in% c(2,4,6)){
    num_Gammas<- 2 * nstrain
    MC_chain<- matrix(NA, nrow=num_iteration, ncol=num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+n_copParams)
    MC_chain[1,]<- c(runif(num_Gammas), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0.1, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_factloadings), rep(0.1, n_copParams))
  }else if(Modeltype==7){
    num_Gammas<- nstate * nstate
    MC_chain<- matrix(NA, nrow=num_iteration, ncol=num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+n_copParams)
    MC_chain[1,]<- c(as.numeric(t(initGs)), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0.1, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_factloadings), rep(0.1, n_copParams))
  }

  Q_r<- MC_chain[1,num_Gammas + 1] * RW2PrecMat
  Q_s<- MC_chain[1,num_Gammas + 2] * RW1PrecMat
  Q_u<- MC_chain[1,num_Gammas + 3] * R

  #Compute gradients
  if(Modeltype %in% c(1,2,5,6)){
    JointTPM<- Multipurpose_JointTransitionMatrix(MC_chain[1,1:num_Gammas], nstrain, MC_chain[1,num_Gammas+3+time+12+ndept+nstrain+nstrain], Modeltype)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
  }else if(Modeltype %in% c(3,4)){
    JointTPM<- Multipurpose_JointTransitionMatrix2(MC_chain[1,1:num_Gammas], nstrain, MC_chain[1,num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)], Modeltype, gh)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
    if(any(!is.finite(JointTPM))) JointTPM<- initGs
  }else if(Modeltype == 7){
    JointTPM<- initGs
  }

  if(Model == 0) JointTPM<- matrix(0, nstate, nstate)

  Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[1, num_Gammas+3+(1:time)], s=MC_chain[1, num_Gammas+3+time+(1:12)], u=MC_chain[1, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total=y_total)
  likelihoodcurrent<- Allquantities$loglike
  priorcurrentRcomps<- randomwalk2_cpp(MC_chain[1, num_Gammas+3+(1:time)], MC_chain[1, num_Gammas+1])
  priorcurrentScomps<- seasonalComp2_cpp(MC_chain[1, num_Gammas+3+time+(1:12)], MC_chain[1, num_Gammas+2], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1_cpp(MC_chain[1, num_Gammas+3+time+12+(1:ndept)], MC_chain[1, num_Gammas+3], R, rankdef)

  grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

  deltaP<- 1
  RMdelta<- 1/(0.234*(1-0.234))

  for (i in 2:num_iteration) {

    MC_chain[i,num_Gammas+1]<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, num_Gammas+3+(1:time)]) %*% RW2PrecMat %*% MC_chain[i-1, num_Gammas+3+(1:time)])/2)
    MC_chain[i,num_Gammas+2]<- rgamma(1, shape = 1 + 11/2, rate = 0.001 + (t(MC_chain[i-1, num_Gammas+3+time+(1:12)]) %*% RW1PrecMat %*% MC_chain[i-1, num_Gammas+3+time+(1:12)])/2)
    MC_chain[i,num_Gammas+3]<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)])/2)

    Q_r<- MC_chain[i,num_Gammas+1] * RW2PrecMat
    Q_s<- MC_chain[i,num_Gammas+2] * RW1PrecMat
    Q_u<- MC_chain[i,num_Gammas+3] * R

    current_r <- MC_chain[i-1, num_Gammas+3+(1:time)]
    current_s <- MC_chain[i-1, num_Gammas+3+time+(1:12)]
    current_u <- MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)]

    #Update s
    Mmatcs<- as.numeric(grad_current$cov_s %*% grad_current$grad_s)

    eps_s <- rnorm(12)
    proposedScomps <- as.numeric(MC_chain[i-1, num_Gammas+3+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs + step_sizes$s * chol(grad_current$cov_s) %*% eps_s)
    proposedScomps<- proposedScomps - mean(proposedScomps)

    Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=current_r, s=proposedScomps, u=current_u, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total=y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- Allquantities$loglike
    q_prop <- mvnfast::dmvn(proposedScomps, mu = MC_chain[i-1, num_Gammas+3+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs, sigma = grad_current$cov_s * step_sizes$s^2, log = TRUE)
    q_curr <- mvnfast::dmvn(MC_chain[i-1, num_Gammas+3+time+(1:12)], mu = proposedScomps + 0.5 * step_sizes$s^2 * Mmatps, sigma = grad_proposed$cov_s * step_sizes$s^2, log = TRUE)
    priorproposedScomps<- seasonalComp2_cpp(proposedScomps, MC_chain[i, num_Gammas+2], RW1PrecMat)

    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop
    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
      MC_chain[i, num_Gammas+3+time+(1:12)]<- proposedScomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentScomps<- priorproposedScomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, num_Gammas+3+time+(1:12)]<- MC_chain[i-1, num_Gammas+3+time+(1:12)]
    }

    #Update r
    eps_r <- rnorm(time)
    Mmatrc<- as.numeric(grad_current$cov_r %*% grad_current$grad_r)

    proposedRcomps <- as.numeric(current_r + 0.5 * step_sizes$r^2 * Mmatrc + step_sizes$r * chol(grad_current$cov_r) %*% eps_r)
    proposedRcomps<- proposedRcomps - mean(proposedRcomps)

    Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=proposedRcomps, s=MC_chain[i, num_Gammas+3+time+(1:12)], u=current_u, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total=y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcomps, mu = current_r + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(current_r, mu = proposedRcomps + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike
    priorproposedRcomps <- randomwalk2_cpp(proposedRcomps, MC_chain[i, num_Gammas+1])

    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop

    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
      MC_chain[i, num_Gammas+3+(1:time)] <- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentRcomps<- priorproposedRcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, num_Gammas+3+(1:time)]<- MC_chain[i-1, num_Gammas+3+(1:time)]
    }

    #Update u
    eps_u <- rnorm(ndept)

    proposedUcomps <- as.numeric(MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u + step_sizes$u * eps_u)
    proposedUcomps<- proposedUcomps - mean(proposedUcomps)

    Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=proposedUcomps, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total=y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    likelihoodproposed<- Allquantities$loglike

    q_prop <- sum(dnorm(proposedUcomps, mean = MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u, sd = step_sizes$u, log = TRUE))
    q_curr <- sum(dnorm(MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)], mean = proposedUcomps + 0.5 * step_sizes$u^2 * grad_proposed$grad_u, sd = step_sizes$u, log = TRUE))

    priorproposedUcomps<- logIGMRF1_cpp(proposedUcomps, MC_chain[i, num_Gammas+3], R, rankdef)
    priorcurrentUcomps<- logIGMRF1_cpp(MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)], MC_chain[i, num_Gammas+3], R, rankdef)

    log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop
    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
      MC_chain[i, num_Gammas+3+time+12+(1:ndept)]<- proposedUcomps
      likelihoodcurrent<- likelihoodproposed
      #priorcurrentUcomps<- priorproposedUcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, num_Gammas+3+time+12+(1:ndept)]<- MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)]
    }

    if(Model == 0){
      MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)]
      MC_chain[i, 1:2]<- MC_chain[i-1, 1:2]
    }else{
      #proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain)))

      proposedB <- rnorm(nstrain, mean = MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain))

      if(any(proposedB < 0)){
        MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)]
      }else{

      priorcurrentB<- sum(dgamma(MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

      Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=proposedB, Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total=y_total)
      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedB #+ proposalcurrentB
                     - likelihoodcurrent - priorcurrentB) #- proposalproposedB
      #print(paste("mh.ratioB = ", mh.ratio))

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- proposedB
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
      }
      else{
        MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)]
      }
      sdBs<- sdBs * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))
    }
      if(Modeltype %in% c(1,2)){
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))

        proposedLambdas<- rep(0, nstrain)
        JointTPM1<- Multipurpose_JointTransitionMatrix(proposedGs, nstrain, proposedLambdas, Modeltype)
        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total=y_total)
        grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

        likelihoodproposed<- Allquantities$loglike

        mh.ratio<- exp(likelihoodproposed + priorproposedGs
                       - likelihoodcurrent - priorcurrentGs)

        #print(mh.ratio)

        if(!is.na(mh.ratio) && runif(1) < mh.ratio){
          MC_chain[i, 1:num_Gammas]<- proposedGs
          likelihoodcurrent<- likelihoodproposed
          grad_current<- grad_proposed
          JointTPM<- JointTPM1
        }
        else{
          MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
        }
        sdGs<- sdGs * exp((2/i) * (min(mh.ratio, 1) - 0.234))
      }else if(Modeltype %in% c(3,4)){
        #Transition probabilities update
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))

        JointTPM1<- Multipurpose_JointTransitionMatrix2(proposedGs, nstrain, MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)], Modeltype, gh)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        if(any(!is.finite(JointTPM1))){
          MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
        }else{

          Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total=y_total)
          grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

          likelihoodproposed<- Allquantities$loglike

          mh.ratioGC<- exp(likelihoodproposed + priorproposedGs
                           - likelihoodcurrent - priorcurrentGs)

          #print(mh.ratioGC)

          if(!is.na(mh.ratioGC) && runif(1) < mh.ratioGC){
            MC_chain[i, 1:num_Gammas]<- proposedGs
            likelihoodcurrent<- likelihoodproposed
            grad_current<- grad_proposed
            JointTPM<- JointTPM1
          }
          else{
            MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
          }
          sdGs<- sdGs * exp((RMdelta/i) * (min(mh.ratioGC, 1) - 0.234))
        }

        #Factor loadings update
        proposedLambdas<- rnorm(nstrain,mean=MC_chain[i-1,num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)], sd=rep(sdLambdas, nstrain))
        #proposedLambdas<- ifelse(proposedLambdas<1, proposedLambdas, 2-proposedLambdas)

        priorcurrentLambdas<- sum(dunif(MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)], min = rep(-1, n_factloadings), max = rep(1, n_factloadings), log=TRUE))
        priorproposedLambdas<- sum(dunif(proposedLambdas, min = rep(-1, n_factloadings), max = rep(1, n_factloadings), log=TRUE))

        JointTPM1<- Multipurpose_JointTransitionMatrix2(MC_chain[i, 1:num_Gammas], nstrain, proposedLambdas, Modeltype, gh)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        if(any(!is.finite(JointTPM1)) || any(abs(proposedLambdas) > 1)){
          MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]
        }else{

          Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total=y_total)
          grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

          likelihoodproposed<- Allquantities$loglike

          mh.ratioGC<- exp(likelihoodproposed + priorproposedLambdas
                           - likelihoodcurrent - priorcurrentLambdas)

          #print(mh.ratioGC)

          if(!is.na(mh.ratioGC) && runif(1) < mh.ratioGC){
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- proposedLambdas
            likelihoodcurrent<- likelihoodproposed
            grad_current<- grad_proposed
            JointTPM<- JointTPM1
          }
          else{
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]
          }
          sdLambdas<- sdLambdas * exp((RMdelta/i) * (min(mh.ratioGC, 1) - 0.234))
        }
      }else if(Modeltype %in% c(5,6)){
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = rep(2,num_Gammas), shape2 = rep(2,num_Gammas), log=TRUE))

        if(nstrain==2){
          proposedcopPs<- rnorm(1,mean=MC_chain[i-1, ncol(MC_chain)], sd=sdCop)
          JointTPM1<- Multipurpose_JointTransitionMatrix(proposedGs, nstrain, proposedcopPs, Modeltype)
        }else{
          proposedcopPs<- rnorm(1,mean=log(MC_chain[i-1, ncol(MC_chain)]), sd=sdCop)
          proposalcurrentcop<- dnorm(log(MC_chain[i-1, ncol(MC_chain)]), mean=proposedcopPs, sd=sdCop, log = T) + MC_chain[i-1, ncol(MC_chain)]
          proposalproposedcop<- dnorm(proposedcopPs, mean=log(MC_chain[i-1, ncol(MC_chain)]), sd=sdCop, log = T) + exp(proposedcopPs)
          JointTPM1<- Multipurpose_JointTransitionMatrix(proposedGs, nstrain, exp(proposedcopPs), Modeltype)
        }

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)
        if(any(!is.finite(JointTPM1))){
          MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
          MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
        }else{
          #JointTPM1 <- JointTPM1 / rowSums(JointTPM1)

          Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, y_total=y_total)
          grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

          likelihoodproposed<- Allquantities$loglike

          if(nstrain==2){
            mh.ratio<- exp(likelihoodproposed + priorproposedGs
                           - likelihoodcurrent - priorcurrentGs)
          }else{
            mh.ratio<- exp(likelihoodproposed + priorproposedGs + proposalcurrentcop
                           - likelihoodcurrent - priorcurrentGs - proposalproposedcop)
          }

          #print(mh.ratio)

          if(!is.na(mh.ratio) && runif(1) < mh.ratio){
            MC_chain[i, 1:num_Gammas]<- proposedGs
            likelihoodcurrent<- likelihoodproposed
            grad_current<- grad_proposed
            JointTPM<- JointTPM1
            MC_chain[i, ncol(MC_chain)]<- ifelse(nstrain==2,proposedcopPs,exp(proposedcopPs))
          }
          else{
            MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
            MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
          }
          sdGs<- sdGs * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))
          sdCop<- sdCop * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))
        }
      }else if(Modeltype==7){

        for(n in 1:nstate){
          index<- nstate * (n-1) + 1

          JointTPM[n, ] <- gtools::rdirichlet(1, rep(1, nstate) + deltaP * MC_chain[i-1, (index:(n*nstate))])

          proposalproposedGs<-  log(gtools::ddirichlet(JointTPM[n, ], MC_chain[i-1, (index:(n*nstate))]))
          proposalcurrentproposedGs<- log(gtools::ddirichlet(MC_chain[i-1, (index:(n*nstate))], JointTPM[n, ]))

          priorcurrentGs<- log(gtools::ddirichlet(MC_chain[i-1, (index:(n*nstate))], rep(1, nstate)))
          priorproposedGs<- log(gtools::ddirichlet(JointTPM[n, ], rep(1, nstate)))

          if(n == nstate){
            Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1, y_total=y_total)
            grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)
          }else{
            Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=0, y_total=y_total)
          }

          likelihoodproposed<- Allquantities$loglike

          mh.ratio<- exp(likelihoodproposed + priorproposedGs + proposalcurrentproposedGs
                         - likelihoodcurrent - priorcurrentGs - proposalproposedGs)

          #print(mh.ratio)

          if(!is.na(mh.ratio) && runif(1) < mh.ratio){
            MC_chain[i, (index:(n*nstate))]<- as.numeric(JointTPM[n, ])
            likelihoodcurrent<- likelihoodproposed
            grad_current<- grad_proposed
            deltaP<- max(0, deltaP-3)
          }
          else{
            MC_chain[i, (index:(n*nstate))]<- MC_chain[i-1, (index:(n*nstate))]
            deltaP<- deltaP + 1
          }
        }
      }
    }
    #Gibbs A_k's update
    if(all(is.finite(as.numeric(Allquantities$poisMean4GibbsUpdate)))){
      MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- log(rgamma(nstrain, shape = 0.01+SumYk_vec, rate = as.numeric(Allquantities$poisMean4GibbsUpdate) + 0.01/exp(-15)))
    }else{
      MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
    }

    if(i %% 1000 == 0) cat("Iteration:", i, "\n")
  }

  if(Modeltype %in% c(0, 1)){
    colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep="")))
  }else if(Modeltype == 2){
    colnames(MC_chain) <- paste(c(paste0(rep(c("G12", "G21"), nstrain), "Strain", rep(1:nstrain,each=2)), "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep="")))
  }else if(Modeltype == 3){
    #Derive pair-correlations from factor-loadings
    copInd<- 1
    for(w in 1:(n_factloadings - 1)){
      for (z in (w + 1):n_factloadings){
        MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+copInd] <- MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+w] * MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+z]
        copInd <- copInd + 1
      }
    }
    colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), paste("FactorLoading", 1:nstrain, sep =""), paste("copulaParam", 1:n_copParams, sep="")))
  }else if(Modeltype == 4){
    #Derive pair-correlations from factor-loadings
    copInd<- 1
    for(w in 1:(n_factloadings - 1)){
      for (z in (w + 1):n_factloadings){
        MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+copInd] <- MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+w] * MC_chain[, num_Gammas+3+time+12+ndept+nstrain+nstrain+z]
        copInd <- copInd + 1
      }
    }
    colnames(MC_chain) <- paste(c(paste0(rep(c("G12", "G21"), nstrain), "Strain", rep(1:nstrain,each=2)), "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), paste("FactorLoading", 1:nstrain, sep =""),paste("copulaParam", 1:n_copParams, sep="")))
  }else if(Modeltype == 5){
    colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), paste("copulaParam", 1:n_copParams, sep ="")))
  }else if(Modeltype == 6){
    colnames(MC_chain) <- paste(c(paste0(rep(c("G12", "G21"), nstrain), "Strain", rep(1:nstrain,each=2)), "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), paste("copulaParam", 1:n_copParams, sep ="")))
  }else if(Modeltype == 7){
    colnames(MC_chain) <- paste(c(paste0(rep("G_", nstate*nstate), rep(1:nstate, each=nstate), ",", 1:nstate), "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep="")))
  }
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}
