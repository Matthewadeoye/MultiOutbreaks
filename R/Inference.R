#Inference Method 1 -- Smoothing
SMOOTHING_INFERENCE<- function(y, e_it, Modeltype, adjmat, step_sizes = list("r"=0.3,"s"=0.3,"u"=0.3),
                               num_iteration = 15000, sdBs=0.003, sdGs=0.005, sdLambdas=0.003, sdCop=0.0002,
                               y_total=NULL, RM_Gs=TRUE, RM_Bs=TRUE, RM_Cop=TRUE, RM_Lambdas=TRUE,
                               burn_in=30000, beta_param1=1, beta_param2=1){
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
      sumY1<- sumY1 + y[,,k]
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
    MC_chain[1,]<- c(rep(0.1,num_Gammas), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_factloadings), rep(0.5, n_copParams))
    shape1params<- rep(c(1,6), num_Gammas)
    shape2params<- rep(c(11,6),num_Gammas)
  }else if(Modeltype %in% c(2,4,6)){
    num_Gammas<- 2 * nstrain
    MC_chain<- matrix(NA, nrow=num_iteration, ncol=num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+n_copParams)
    MC_chain[1,]<- c(rep(0.1,num_Gammas), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_factloadings), rep(0.5, n_copParams))
    shape1params<- rep(c(1,6), num_Gammas)
    shape2params<- rep(c(11,6),num_Gammas)
  }else if(Modeltype==7){
    num_Gammas<- nstate * nstate
    MC_chain<- matrix(NA, nrow=num_iteration, ncol=num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+n_copParams)
    MC_chain[1,]<- c(as.numeric(t(initGs)), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_factloadings), rep(0.1, n_copParams))
  }

  Q_r<- MC_chain[1,num_Gammas + 1] * RW2PrecMat
  Q_s<- MC_chain[1,num_Gammas + 2] * RW1PrecMat
  Q_u<- MC_chain[1,num_Gammas + 3] * R

  Qstz_r<- QRstz_basis(time)
  Qstz_s<- QRstz_basis(12)
  Qstz_u<- QRstz_basis(ndept)

  #Compute gradients
  if(Modeltype %in% c(1,2,5,6)){
    JointTPM<- Multipurpose_JointTransitionMatrix(MC_chain[1,1:num_Gammas], nstrain, MC_chain[1,num_Gammas+3+time+12+ndept+nstrain+nstrain], Modeltype)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
  }else if(Modeltype %in% c(3,4)){
    eta <- atanh(MC_chain[1,num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)])
    JointTPM<- Multipurpose_JointTransitionMatrix2(MC_chain[1,1:num_Gammas], nstrain, MC_chain[1,num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)], Modeltype, gh)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
    if(any(!is.finite(JointTPM))) JointTPM<- initGs
  }else if(Modeltype == 7){
    OutbreakPrior_ExpectationMatrix<- matrix(c(11/12,1/12,6/12,6/12), nrow = 2, byrow = T)
    Dirichlet_Prior<- 12 * JointTransitionMatrix_arma_cpp(OutbreakPrior_ExpectationMatrix, nstrain)
    JointTPM<- initGs
  }

  if(Model == 0) JointTPM<- matrix(0, nstate, nstate)

  Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[1, num_Gammas+3+(1:time)], s=MC_chain[1, num_Gammas+3+time+(1:12)], u=MC_chain[1, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
  likelihoodcurrent<- Allquantities$loglike
  priorcurrentRcomps<- randomwalk2_cpp(MC_chain[1, num_Gammas+3+(1:time)], MC_chain[1, num_Gammas+1])
  priorcurrentScomps<- seasonalComp2_cpp(MC_chain[1, num_Gammas+3+time+(1:12)], MC_chain[1, num_Gammas+2], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1_cpp(MC_chain[1, num_Gammas+3+time+12+(1:ndept)], MC_chain[1, num_Gammas+3], R, rankdef)

  grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
  current_rtrans<- t(Qstz_r) %*% MC_chain[1, num_Gammas+3+(1:time)]
  current_strans<- t(Qstz_s) %*% MC_chain[1, num_Gammas+3+time+(1:12)]
  current_utrans<- t(Qstz_u) %*% MC_chain[1, num_Gammas+3+time+12+(1:ndept)]

  deltaP<- 1
  RMdelta<- 1/(0.234*(1-0.234))  #for Betas, Gammas
  RMLdelta<- 1/(0.44*(1-0.44))   #for factor_loadings
  sdLambdas<- rep(sdLambdas, nstrain)
  sdLambdasJoint<- sdLambdas
  sdGsJoint<- sdGs

  if(Modeltype %in% c(1,2)){
    zigmaJ <- diag(0.0001, num_Gammas+nstrain+nstrain)
  }else if(Modeltype %in% c(3,4)){
    zigmaJ <- diag(0.0001, num_Gammas+nstrain+nstrain+n_factloadings)
  }else if(Modeltype %in% c(5,6)){
    zigmaJ <- diag(0.0001, num_Gammas+nstrain+nstrain+n_copParams)
  }else if(Modeltype==7){
    zigmaJ <- diag(0.0001, nstrain+nstrain)
  }
  acc<- 0

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

    eps_s <- rnorm(11)
    proposedScompst <- as.numeric(current_strans + 0.5 * step_sizes$s^2 * Mmatcs + step_sizes$s * chol(grad_current$cov_s) %*% eps_s)

    proposedScomps<- Qstz_s %*% proposedScompst
    Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=current_r, s=proposedScomps, u=current_u, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- Allquantities$loglike
    q_prop <- mvnfast::dmvn(proposedScompst, mu = current_strans + 0.5 * step_sizes$s^2 * Mmatcs, sigma = grad_current$cov_s * step_sizes$s^2, log = TRUE)
    q_curr <- mvnfast::dmvn(t(current_strans), mu = proposedScompst + 0.5 * step_sizes$s^2 * Mmatps, sigma = grad_proposed$cov_s * step_sizes$s^2, log = TRUE)
    priorproposedScomps<- seasonalComp2_cpp(proposedScomps, MC_chain[i, num_Gammas+2], RW1PrecMat)

    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop
    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
      MC_chain[i, num_Gammas+3+time+(1:12)]<- proposedScomps
      likelihoodcurrent<- likelihoodproposed
      current_strans<- proposedScompst
      priorcurrentScomps<- priorproposedScomps
      grad_current<- grad_proposed
      Allquantities<- Allquantities
    }else{
      MC_chain[i, num_Gammas+3+time+(1:12)]<- MC_chain[i-1, num_Gammas+3+time+(1:12)]
    }

    #Update r
    eps_r <- rnorm(time-1)
    Mmatrc<- as.numeric(grad_current$cov_r %*% grad_current$grad_r)

    proposedRcompst <- as.numeric(current_rtrans + 0.5 * step_sizes$r^2 * Mmatrc + step_sizes$r * chol(grad_current$cov_r) %*% eps_r)

    proposedRcomps<- Qstz_r %*% proposedRcompst
    Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=proposedRcomps, s=MC_chain[i, num_Gammas+3+time+(1:12)], u=current_u, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcompst, mu = current_rtrans + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(t(current_rtrans), mu = proposedRcompst + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike
    priorproposedRcomps <- randomwalk2_cpp(proposedRcomps, MC_chain[i, num_Gammas+1])

    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop

    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
      MC_chain[i, num_Gammas+3+(1:time)] <- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
      current_rtrans<- proposedRcompst
      priorcurrentRcomps<- priorproposedRcomps
      grad_current<- grad_proposed
      Allquantities<- Allquantities
    }else{
      MC_chain[i, num_Gammas+3+(1:time)]<- MC_chain[i-1, num_Gammas+3+(1:time)]
    }

    #Update u
    eps_u <- rnorm(ndept-1)
    Mmatuc<- as.numeric(grad_current$cov_u %*% grad_current$grad_u)

    proposedUcompst <- as.numeric(current_utrans + 0.5 * step_sizes$u^2 * Mmatuc + step_sizes$u * chol(grad_current$cov_u) %*% eps_u)

    proposedUcomps<- Qstz_u %*% proposedUcompst
    Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=proposedUcomps, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)

    Mmatup<- as.numeric(grad_proposed$cov_u %*% grad_proposed$grad_u)

    q_prop <- mvnfast::dmvn(proposedUcompst, mu = current_utrans + 0.5 * step_sizes$u^2 * Mmatuc, sigma = grad_current$cov_u * step_sizes$u^2, log = TRUE)
    q_curr <- mvnfast::dmvn(t(current_utrans), mu = proposedUcompst + 0.5 * step_sizes$u^2 * Mmatup, sigma = grad_proposed$cov_u * step_sizes$u^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike

    priorproposedUcomps<- logIGMRF1_cpp(proposedUcomps, MC_chain[i, num_Gammas+3], R, rankdef)

    log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop
    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
      MC_chain[i, num_Gammas+3+time+12+(1:ndept)]<- proposedUcomps
      current_utrans<- proposedUcompst
      likelihoodcurrent<- likelihoodproposed
      priorcurrentUcomps<- priorproposedUcomps
      grad_current<- grad_proposed
      Allquantities<- Allquantities
    }else{
      MC_chain[i, num_Gammas+3+time+12+(1:ndept)]<- MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)]
    }

    #Pseudo-Gibbs A_k's proposal
    propRate<- as.numeric(Allquantities$poisMean4GibbsUpdate) + 0.01/exp(-15)
    proposedAks<- log(rgamma(nstrain, shape = propShape, rate = propRate))

    Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=proposedAks, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)

    likelihoodproposed<- Allquantities$loglike

    proposalcurrentAks <- sum(dgamma(exp(MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]),
                                     shape=propShape, rate=propRate, log=TRUE) + MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]) #Jacobian
    proposalproposedAks <- sum(dgamma(exp(proposedAks),
                                      shape=propShape, rate=propRate, log=TRUE) + proposedAks)  #Jacobian

    priorcurrentAks <- sum(dgamma(exp(MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]),
                                  shape=0.01, rate=0.01/exp(-15), log=TRUE) + MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])
    priorproposedAks <- sum(dgamma(exp(proposedAks),
                                   shape=0.01, rate=0.01/exp(-15), log=TRUE) +  proposedAks)

    mh.ratioAk <- (likelihoodproposed + priorproposedAks + proposalcurrentAks
                   - likelihoodcurrent - priorcurrentAks - proposalproposedAks)

    if (is.finite(mh.ratioAk) && log(runif(1)) < mh.ratioAk){
      MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- proposedAks
      likelihoodcurrent<- likelihoodproposed
      grad_current<- grad_proposed
      Allquantities<- Allquantities
    }else{
      MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
    }

    if(Model == 0){
      MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)]
      MC_chain[i, 1:2]<- MC_chain[i-1, 1:2]
    }else{
      proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain)))

      priorcurrentB<- sum(dgamma(MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

      Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=proposedB, Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedB
                     - likelihoodcurrent - priorcurrentB)

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- proposedB
        likelihoodcurrent<- likelihoodproposed
      }
      else{
        MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)]
      }
      if(RM_Bs && i<burn_in && !is.na(mh.ratio)) {sdBs= sdBs * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))}

      if(Modeltype %in% c(1,2)){
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = shape1params, shape2 = shape2params, log=TRUE))

        proposedLambdas<- rep(0, nstrain)
        JointTPM1<- Multipurpose_JointTransitionMatrix(proposedGs, nstrain, proposedLambdas, Modeltype)
        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
        grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)

        likelihoodproposed<- Allquantities$loglike

        mh.ratio<- exp(likelihoodproposed + priorproposedGs
                       - likelihoodcurrent - priorcurrentGs)

        if(!is.na(mh.ratio) && runif(1) < mh.ratio){
          MC_chain[i, 1:num_Gammas]<- proposedGs
          likelihoodcurrent<- likelihoodproposed
          grad_current<- grad_proposed
          JointTPM<- JointTPM1
        }
        else{
          MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
        }
        if(RM_Gs && i<burn_in && !is.na(mh.ratio)) {sdGs= sdGs * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))}
      }else if(Modeltype %in% c(3,4)){
        #Transition probabilities update
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = shape1params, shape2 = shape2params, log=TRUE))

        currentL<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)]
        if(currentL[1]<0) currentL= -currentL

        JointTPM1<- Multipurpose_JointTransitionMatrix2(proposedGs, nstrain, currentL, Modeltype, gh)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        if(any(!is.finite(JointTPM1))){
          MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
          if(RM_Gs && i<burn_in) {sdGs= sdGs * exp((RMdelta/i) * (0 - 0.234))}
          sdGs<- max(sdGs, 1e-6)
        }else{

          Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)

          likelihoodproposed<- Allquantities$loglike

          mh.ratioGC<- exp(likelihoodproposed + priorproposedGs
                           - likelihoodcurrent - priorcurrentGs)

          if(!is.na(mh.ratioGC) && runif(1) < mh.ratioGC){
            MC_chain[i, 1:num_Gammas]<- proposedGs
            likelihoodcurrent<- likelihoodproposed
            JointTPM<- JointTPM1
          }
          else{
            MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
          }
          if(RM_Gs && i<burn_in && !is.na(mh.ratioGC)) {sdGs= sdGs * exp((RMdelta/i) * (min(mh.ratioGC, 1) - 0.234))}
        }

        #FactorLoadings joint update
        eta_prop <- rnorm(n_factloadings, eta, sdLambdasJoint)
        LAMBDAS_prop <- tanh(eta_prop)
        if(!is.na(LAMBDAS_prop[1]) && LAMBDAS_prop[1] < 0){
          LAMBDAS_prop<- -LAMBDAS_prop
        }

        JointTPM1<- Multipurpose_JointTransitionMatrix2(MC_chain[i,1:num_Gammas], nstrain, LAMBDAS_prop, Modeltype, gh)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        if(any(!is.finite(JointTPM1))){
          MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]
          if(RM_Lambdas && i<burn_in) {sdLambdasJoint= sdLambdasJoint * exp((RMdelta/i) * (0 - 0.234))}
          sdLambdasJoint<- max(sdLambdasJoint, 1e-6)
        }else{

          Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)

          likelihoodproposed<- Allquantities$loglike

          priorcurrentLAMBDA<- sum(log_GDP(eta, 3, 1))
          priorproposedLAMBDA<- sum(log_GDP(eta_prop, 3, 1))

          mh.ratioGC<- exp(likelihoodproposed + priorproposedLAMBDA
                           - likelihoodcurrent - priorcurrentLAMBDA)

          if(!is.na(mh.ratioGC) && runif(1) < mh.ratioGC){
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- LAMBDAS_prop
            eta<- eta_prop
            likelihoodcurrent<- likelihoodproposed
            JointTPM<- JointTPM1
          }
          else{
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]
          }
          if(RM_Lambdas && i<burn_in && !is.na(mh.ratioGC)) {sdLambdasJoint= sdLambdasJoint * exp((RMdelta/i) * (min(mh.ratioGC, 1) - 0.234))}
        }
        Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)
        grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
        likelihoodcurrent<- Allquantities$loglike
      }else if(Modeltype %in% c(5,6)){
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = shape1params, shape2 = shape2params, log=TRUE))

        JointTPM1<- Multipurpose_JointTransitionMatrix(proposedGs, nstrain, MC_chain[i-1, ncol(MC_chain)], Modeltype)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)
        if(any(!is.finite(JointTPM1))){
          MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
          if(RM_Gs && i<burn_in) {sdGs= sdGs * exp((RMdelta/i) * (0 - 0.234))}
          sdGs<- max(sdGs, 1e-6)
        }else{

          Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)

          likelihoodproposed<- Allquantities$loglike

          mh.ratio<- exp(likelihoodproposed + priorproposedGs
                           - likelihoodcurrent - priorcurrentGs)

          if(!is.na(mh.ratio) && runif(1) < mh.ratio){
            MC_chain[i, 1:num_Gammas]<- proposedGs
            likelihoodcurrent<- likelihoodproposed
            JointTPM<- JointTPM1
          }
          else{
            MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
          }
          if(RM_Gs && i<burn_in && !is.na(mh.ratio)) {sdGs= sdGs * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))}
        }

          #simple random-walk update of Phi
        proposedcopPs<- rnorm(1,mean=MC_chain[i-1, ncol(MC_chain)], sd=sdCop)

        if(nstrain != 2 && proposedcopPs < 0){
          MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
          if(RM_Cop && i<burn_in) {sdCop= sdCop * exp((RMLdelta/i) * (0 - 0.44))}
          sdCop<- max(sdCop, 1e-6)
        }else{

        JointTPM1<- Multipurpose_JointTransitionMatrix(MC_chain[i, 1:num_Gammas], nstrain, proposedcopPs, Modeltype)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)
        if(any(!is.finite(JointTPM1))){
          MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
          if(RM_Cop && i<burn_in) {sdCop= sdCop * exp((RMLdelta/i) * (0 - 0.44))}
          sdCop<- max(sdCop, 1e-6)
        }else{

          Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0, Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)

          likelihoodproposed<- Allquantities$loglike

          if(nstrain==2){
            priorproposedcopPs<- dnorm(proposedcopPs, mean = 0, sd=100, log = TRUE)
            priorcurrentcopPs<- dnorm(MC_chain[i-1, ncol(MC_chain)], mean=0, sd=100, log = TRUE)
          }else{
            priorproposedcopPs<- dexp(proposedcopPs, rate=0.5, log = TRUE)
            priorcurrentcopPs<- dexp(MC_chain[i-1, ncol(MC_chain)], rate = 0.5, log = TRUE)
          }

          mh.ratio<- exp(likelihoodproposed + priorproposedcopPs - likelihoodcurrent - priorcurrentcopPs)

          if(!is.na(mh.ratio) && runif(1) < mh.ratio){
            likelihoodcurrent<- likelihoodproposed
            JointTPM<- JointTPM1
            MC_chain[i, ncol(MC_chain)]<- proposedcopPs
          }
          else{
            MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
          }
          if(RM_Cop && i<burn_in && !is.na(mh.ratio)) {sdCop= sdCop * exp((RMLdelta/i) * (min(mh.ratio, 1) - 0.44))}
          }
        }
        Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
        grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
        likelihoodcurrent<- Allquantities$loglike
      }else if(Modeltype==7){

        for(n in 1:nstate){
          index<- nstate * (n-1) + 1

          JointTPM1<- JointTPM
          state_n_Prior<- Dirichlet_Prior[n, ]
          JointTPM1[n, ]<- gtools::rdirichlet(1, state_n_Prior + deltaP * MC_chain[i-1, (index:(n*nstate))])

          proposalproposedGs<-  log(gtools::ddirichlet(JointTPM1[n, ], MC_chain[i-1, (index:(n*nstate))]))
          proposalcurrentproposedGs<- log(gtools::ddirichlet(MC_chain[i-1, (index:(n*nstate))], JointTPM1[n, ]))

          priorcurrentGs<- log(gtools::ddirichlet(MC_chain[i-1, (index:(n*nstate))], state_n_Prior))
          priorproposedGs<- log(gtools::ddirichlet(JointTPM1[n, ], state_n_Prior))

          Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)

          likelihoodproposed<- Allquantities$loglike

          mh.ratio<- exp(likelihoodproposed + priorproposedGs + proposalcurrentproposedGs
                         - likelihoodcurrent - priorcurrentGs - proposalproposedGs)

          if(!is.na(mh.ratio) && runif(1) < mh.ratio){
            MC_chain[i, (index:(n*nstate))]<- as.numeric(JointTPM1[n, ])
            JointTPM<- JointTPM1
            likelihoodcurrent<- likelihoodproposed
            deltaP<- max(0, deltaP-3)
          }
          else{
            MC_chain[i, (index:(n*nstate))]<- MC_chain[i-1, (index:(n*nstate))]
            deltaP<- deltaP + 1
          }
        }
        JointTPM<- matrix(MC_chain[i, 1:num_Gammas], nrow = nstate, ncol = nstate, byrow = TRUE)
        Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
        grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
        likelihoodcurrent<- Allquantities$loglike
      }
      #Joint update of gammas, intercepts, and regression parameters
      if(Modeltype %in% c(1,2)){
        #Adapting zigmaJ
        if(i==2000){
          acc<- 0
          optconstantJ<- 2.38^2/(num_Gammas+nstrain+nstrain)
          lambdaJ<- 1
          epsilonJ<- 1e-6
          half_hist<- 0.5*i
          XnJ<- cbind(MC_chain[half_hist:i, 1:num_Gammas], MC_chain[half_hist:i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[half_hist:i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])
          XnbarJ <- colMeans(XnJ)
          zigmaJ <- cov(XnJ) + epsilonJ * diag(1, num_Gammas+nstrain+nstrain)
          zigmaJ<- optconstantJ * zigmaJ
        } else if (i > 5){
          currentJcomps<- c(MC_chain[i, 1:num_Gammas], MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])
          proposedJcomps<- mvnfast::rmvn(1, mu = currentJcomps, sigma = zigmaJ)

          if(any(proposedJcomps[1:num_Gammas]<0) || any(proposedJcomps[1:num_Gammas]>1)){
            MC_chain[i,1:num_Gammas]<- MC_chain[i,1:num_Gammas]
            MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
          }else{
          JointTPM1<- Multipurpose_JointTransitionMatrix(proposedJcomps[1:num_Gammas], nstrain, MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_copParams)], Modeltype)

          JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
          JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)
          if(any(!is.finite(JointTPM1))){
            MC_chain[i,1:num_Gammas]<- MC_chain[i,1:num_Gammas]
            MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
          }else{
          priorcurrentGs<- sum(dbeta(MC_chain[i,1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))
          priorproposedGs<- sum(dbeta(proposedJcomps[1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))

          priorcurrentB<- sum(dgamma(MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
          priorproposedB<- sum(dgamma(proposedJcomps[num_Gammas+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

          priorcurrentAks <- sum(dgamma(exp(MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]),
                                        shape=0.01, rate=0.01/exp(-15), log=TRUE) + MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])
          priorproposedAks <- sum(dgamma(exp(proposedJcomps[num_Gammas+nstrain+(1:nstrain)]),
                                         shape=0.01, rate=0.01/exp(-15), log=TRUE) +  proposedJcomps[num_Gammas+nstrain+(1:nstrain)])

          proposalproposedJcomps<- mvnfast::dmvn(proposedJcomps, mu = currentJcomps, sigma = zigmaJ, log = TRUE)
          proposalcurrentJcomps<- mvnfast::dmvn(currentJcomps, mu = proposedJcomps, sigma = zigmaJ, log = TRUE)

          Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=proposedJcomps[num_Gammas+(1:nstrain)], Bits=Bits, a_k=proposedJcomps[num_Gammas+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
          grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
          likelihoodproposed<- Allquantities$loglike

          mh.ratioJ<- exp(likelihoodproposed + priorproposedGs + priorproposedB + priorproposedAks + proposalcurrentJcomps
                          - likelihoodcurrent - priorcurrentGs - priorcurrentB - priorcurrentAks - proposalproposedJcomps)

          if(!is.na(mh.ratioJ) && runif(1) < mh.ratioJ){
            MC_chain[i,1:num_Gammas]<- proposedJcomps[1:num_Gammas]
            MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- proposedJcomps[num_Gammas+(1:nstrain)]
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- proposedJcomps[num_Gammas+nstrain+(1:nstrain)]
            grad_current<- grad_proposed
            JointTPM<- JointTPM1
            likelihoodcurrent<- likelihoodproposed
            acc<- acc + 1
          }else{
            MC_chain[i,1:num_Gammas]<- MC_chain[i,1:num_Gammas]
            MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
          }
          if(i>2000 && i<burn_in){
          currentJcomps<- c(MC_chain[i, 1:num_Gammas], MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])

          XnbarPrevJ <- XnbarJ
          XnbarJ <- (i*XnbarJ + currentJcomps)/(i+1)
          zigmaJ <- ((i-1)*zigmaJ + tcrossprod(currentJcomps) + i*tcrossprod(XnbarPrevJ) - (i+1)*tcrossprod(XnbarJ) + epsilonJ*diag(1,num_Gammas+nstrain+nstrain))/i
          #Robbins Munro tuning
          lambdaJ<- lambdaJ * exp((2/max(1, i-2000)) * (min(mh.ratioJ, 1) - 0.234))
          zigmaJ<- lambdaJ* optconstantJ * zigmaJ
              }
            }
          }
        }
      }else if(Modeltype %in% c(3,4)){
        #Adapting zigmaJ
        if(i==2000){
          acc<- 0
          optconstantJ<- 2.38^2/(num_Gammas+nstrain+nstrain+n_factloadings)
          lambdaJ<- 1
          epsilonJ<- 1e-6
          half_hist<- 0.5*i
          eta_matrix <- atanh(MC_chain[half_hist:i,num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)])
          XnJ<- cbind(MC_chain[half_hist:i, 1:num_Gammas], MC_chain[half_hist:i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[half_hist:i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], eta_matrix)
          XnbarJ <- colMeans(XnJ)
          zigmaJ <- cov(XnJ) + epsilonJ * diag(1, num_Gammas+nstrain+nstrain+n_factloadings)
          zigmaJ<- optconstantJ * zigmaJ
        } else if (i > 5){
          currentJcomps<- c(MC_chain[i, 1:num_Gammas], MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], eta)
          proposedJcomps<- mvnfast::rmvn(1, mu = currentJcomps, sigma = zigmaJ)
          LAMBDAS_prop <- tanh(proposedJcomps[num_Gammas+nstrain+nstrain+(1:n_factloadings)])
          if(!is.na(LAMBDAS_prop[1]) && LAMBDAS_prop[1] < 0){
            LAMBDAS_prop<- -LAMBDAS_prop
            }

          if(any(proposedJcomps[1:num_Gammas]<0) || any(proposedJcomps[1:num_Gammas]>1)){
            MC_chain[i,1:num_Gammas]<- MC_chain[i,1:num_Gammas]
            MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]
          }else{
            JointTPM1<- Multipurpose_JointTransitionMatrix2(proposedJcomps[1:num_Gammas], nstrain, LAMBDAS_prop, Modeltype, gh)

            JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
            JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)
            if(any(!is.finite(JointTPM1))){
              MC_chain[i,1:num_Gammas]<- MC_chain[i,1:num_Gammas]
              MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]
              MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
              MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]
            }else{
              priorcurrentGs<- sum(dbeta(MC_chain[i,1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))
              priorproposedGs<- sum(dbeta(proposedJcomps[1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))

              priorcurrentB<- sum(dgamma(MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
              priorproposedB<- sum(dgamma(proposedJcomps[num_Gammas+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

              priorcurrentAks <- sum(dgamma(exp(MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]),
                                            shape=0.01, rate=0.01/exp(-15), log=TRUE) + MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])
              priorproposedAks <- sum(dgamma(exp(proposedJcomps[num_Gammas+nstrain+(1:nstrain)]),
                                             shape=0.01, rate=0.01/exp(-15), log=TRUE) +  proposedJcomps[num_Gammas+nstrain+(1:nstrain)])

              priorcurrentLAMBDA<- sum(log_GDP(eta, 3, 1))
              priorproposedLAMBDA<- sum(log_GDP(proposedJcomps[num_Gammas+nstrain+nstrain+(1:n_factloadings)], 3, 1))

              proposalproposedJcomps<- mvnfast::dmvn(proposedJcomps, mu = currentJcomps, sigma = zigmaJ, log = TRUE)
              proposalcurrentJcomps<- mvnfast::dmvn(currentJcomps, mu = proposedJcomps, sigma = zigmaJ, log = TRUE)

              Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=proposedJcomps[num_Gammas+(1:nstrain)], Bits=Bits, a_k=proposedJcomps[num_Gammas+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
              grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
              likelihoodproposed<- Allquantities$loglike

              mh.ratioJ<- exp(likelihoodproposed + priorproposedGs + priorproposedB + priorproposedAks +  priorproposedLAMBDA + proposalcurrentJcomps
                              - likelihoodcurrent - priorcurrentGs - priorcurrentB - priorcurrentAks - priorcurrentLAMBDA - proposalproposedJcomps)

              if(!is.na(mh.ratioJ) && runif(1) < mh.ratioJ){
                MC_chain[i,1:num_Gammas]<- proposedJcomps[1:num_Gammas]
                MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- proposedJcomps[num_Gammas+(1:nstrain)]
                MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- proposedJcomps[num_Gammas+nstrain+(1:nstrain)]
                MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- LAMBDAS_prop
                eta<- proposedJcomps[num_Gammas+nstrain+nstrain+(1:n_factloadings)]
                grad_current<- grad_proposed
                JointTPM<- JointTPM1
                likelihoodcurrent<- likelihoodproposed
                acc<- acc + 1
              }else{
                MC_chain[i,1:num_Gammas]<- MC_chain[i,1:num_Gammas]
                MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]
                MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
                MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]
              }
              if(i>2000 && i<burn_in){
              currentJcomps<- c(MC_chain[i, 1:num_Gammas], MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], eta)
              XnbarPrevJ <- XnbarJ
              XnbarJ <- (i*XnbarJ + currentJcomps)/(i+1)
              zigmaJ <- ((i-1)*zigmaJ + tcrossprod(currentJcomps) + i*tcrossprod(XnbarPrevJ) - (i+1)*tcrossprod(XnbarJ) + epsilonJ*diag(1,num_Gammas+nstrain+nstrain+n_factloadings))/i
              #Robbins Munro tuning
              lambdaJ<- lambdaJ * exp((2/max(1, i-2000)) * (min(mh.ratioJ, 1) - 0.234))
              zigmaJ<- lambdaJ* optconstantJ * zigmaJ
              }
            }
          }
        }
      }else if(Modeltype %in% c(5,6)){
        #Adapting zigmaJ
        if(i==2000){
          acc<- 0
          optconstantJ<- 2.38^2/(num_Gammas+nstrain+nstrain+1)
          lambdaJ<- 1
          epsilonJ<- 1e-6
          half_hist<- 0.5*i
          XnJ<- cbind(MC_chain[half_hist:i, 1:num_Gammas], MC_chain[half_hist:i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[half_hist:i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], MC_chain[half_hist:i, ncol(MC_chain)])
          XnbarJ <- colMeans(XnJ)
          zigmaJ <- cov(XnJ) + epsilonJ * diag(1, num_Gammas+nstrain+nstrain+1)
          zigmaJ<- optconstantJ * zigmaJ
        } else if (i > 5){
          currentJcomps<- c(MC_chain[i, 1:num_Gammas], MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], MC_chain[i, ncol(MC_chain)])
          proposedJcomps<- mvnfast::rmvn(1, mu = currentJcomps, sigma = zigmaJ)

          if(any(proposedJcomps[1:num_Gammas]<0) || any(proposedJcomps[1:num_Gammas]>1) || (nstrain > 2 && proposedJcomps[length(proposedJcomps)]<0)){
            MC_chain[i,1:num_Gammas]<- MC_chain[i,1:num_Gammas]
            MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
            MC_chain[i, ncol(MC_chain)]<- MC_chain[i, ncol(MC_chain)]
          }else{
            JointTPM1<- Multipurpose_JointTransitionMatrix(proposedJcomps[1:num_Gammas], nstrain, proposedJcomps[length(proposedJcomps)], Modeltype)

            JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
            JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)
            if(any(!is.finite(JointTPM1))){
              MC_chain[i,1:num_Gammas]<- MC_chain[i,1:num_Gammas]
              MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]
              MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
              MC_chain[i, ncol(MC_chain)]<- MC_chain[i, ncol(MC_chain)]
            }else{
              priorcurrentGs<- sum(dbeta(MC_chain[i,1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))
              priorproposedGs<- sum(dbeta(proposedJcomps[1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))

              priorcurrentB<- sum(dgamma(MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
              priorproposedB<- sum(dgamma(proposedJcomps[num_Gammas+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

              priorcurrentAks <- sum(dgamma(exp(MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]),
                                            shape=0.01, rate=0.01/exp(-15), log=TRUE) + MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])
              priorproposedAks <- sum(dgamma(exp(proposedJcomps[num_Gammas+nstrain+(1:nstrain)]),
                                             shape=0.01, rate=0.01/exp(-15), log=TRUE) +  proposedJcomps[num_Gammas+nstrain+(1:nstrain)])

              if(nstrain==2){
                priorproposedcopPs<- dnorm(proposedJcomps[length(proposedJcomps)], mean = 0, sd=100, log = TRUE)
                priorcurrentcopPs<- dnorm(MC_chain[i, ncol(MC_chain)], mean=0, sd=100, log = TRUE)
              }else{
                priorproposedcopPs<- dexp(proposedJcomps[length(proposedJcomps)], rate=0.5, log = TRUE)
                priorcurrentcopPs<- dexp(MC_chain[i, ncol(MC_chain)], rate = 0.5, log = TRUE)
              }

              proposalproposedJcomps<- mvnfast::dmvn(proposedJcomps, mu = currentJcomps, sigma = zigmaJ, log = TRUE)
              proposalcurrentJcomps<- mvnfast::dmvn(currentJcomps, mu = proposedJcomps, sigma = zigmaJ, log = TRUE)

              Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=proposedJcomps[num_Gammas+(1:nstrain)], Bits=Bits, a_k=proposedJcomps[num_Gammas+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
              grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
              likelihoodproposed<- Allquantities$loglike

              mh.ratioJ<- exp(likelihoodproposed + priorproposedGs + priorproposedB + priorproposedAks + priorproposedcopPs + proposalcurrentJcomps
                              - likelihoodcurrent - priorcurrentGs - priorcurrentB - priorcurrentAks - priorcurrentcopPs - proposalproposedJcomps)

              if(!is.na(mh.ratioJ) && runif(1) < mh.ratioJ){
                MC_chain[i,1:num_Gammas]<- proposedJcomps[1:num_Gammas]
                MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- proposedJcomps[num_Gammas+(1:nstrain)]
                MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- proposedJcomps[num_Gammas+nstrain+(1:nstrain)]
                MC_chain[i, ncol(MC_chain)]<- proposedJcomps[length(proposedJcomps)]
                grad_current<- grad_proposed
                JointTPM<- JointTPM1
                likelihoodcurrent<- likelihoodproposed
                acc<- acc + 1
              }else{
                MC_chain[i,1:num_Gammas]<- MC_chain[i,1:num_Gammas]
                MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]
                MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
                MC_chain[i, ncol(MC_chain)]<- MC_chain[i, ncol(MC_chain)]
              }
              if(i>2000 && i<burn_in){
              currentJcomps<- c(MC_chain[i, 1:num_Gammas], MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], MC_chain[i, ncol(MC_chain)])

              XnbarPrevJ <- XnbarJ
              XnbarJ <- (i*XnbarJ + currentJcomps)/(i+1)
              zigmaJ <- ((i-1)*zigmaJ + tcrossprod(currentJcomps) + i*tcrossprod(XnbarPrevJ) - (i+1)*tcrossprod(XnbarJ) + epsilonJ*diag(1,num_Gammas+nstrain+nstrain+1))/i
              #Robbins Munro tuning
              lambdaJ<- lambdaJ * exp((2/max(1, i-2000)) * (min(mh.ratioJ, 1) - 0.234))
              zigmaJ<- lambdaJ* optconstantJ * zigmaJ
              }
            }
          }
        }
      }else if(Modeltype == 7){
        #Adapting zigmaJ
        if(i==2000){
          acc<- 0
          optconstantJ<- 2.38^2/(nstrain+nstrain)
          lambdaJ<- 1
          epsilonJ<- 1e-6
          half_hist<- 0.5*i
          XnJ<- cbind(MC_chain[half_hist:i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[half_hist:i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])
          XnbarJ <- colMeans(XnJ)
          zigmaJ <- cov(XnJ) + epsilonJ * diag(1, nstrain+nstrain)
          zigmaJ<- optconstantJ * zigmaJ
        } else if (i > 5){
          currentJcomps<- c(MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])
          proposedJcomps<- mvnfast::rmvn(1, mu = currentJcomps, sigma = zigmaJ)

          priorcurrentB<- sum(dgamma(MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
          priorproposedB<- sum(dgamma(proposedJcomps[1:nstrain], shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

          priorcurrentAks <- sum(dgamma(exp(MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]),
                                            shape=0.01, rate=0.01/exp(-15), log=TRUE) + MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])
          priorproposedAks <- sum(dgamma(exp(proposedJcomps[nstrain+(1:nstrain)]),
                                             shape=0.01, rate=0.01/exp(-15), log=TRUE) +  proposedJcomps[nstrain+(1:nstrain)])

          proposalproposedJcomps<- mvnfast::dmvn(proposedJcomps, mu = currentJcomps, sigma = zigmaJ, log = TRUE)
          proposalcurrentJcomps<- mvnfast::dmvn(currentJcomps, mu = proposedJcomps, sigma = zigmaJ, log = TRUE)

          Allquantities<- SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=proposedJcomps[1:nstrain], Bits=Bits, a_k=proposedJcomps[nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
          grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
          likelihoodproposed<- Allquantities$loglike

          mh.ratioJ<- exp(likelihoodproposed + priorproposedGs + priorproposedB + priorproposedAks + proposalcurrentJcomps
                              - likelihoodcurrent - priorcurrentGs - priorcurrentB - priorcurrentAks - proposalproposedJcomps)

              if(!is.na(mh.ratioJ) && runif(1) < mh.ratioJ){
                MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- proposedJcomps[1:nstrain]
                MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- proposedJcomps[nstrain+(1:nstrain)]
                grad_current<- grad_proposed
                likelihoodcurrent<- likelihoodproposed
                acc<- acc + 1
              }else{
                MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]
                MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]
              }
              if(i>2000 && i<burn_in){
              currentJcomps<- c(MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)])

              XnbarPrevJ <- XnbarJ
              XnbarJ <- (i*XnbarJ + currentJcomps)/(i+1)
              zigmaJ <- ((i-1)*zigmaJ + tcrossprod(currentJcomps) + i*tcrossprod(XnbarPrevJ) - (i+1)*tcrossprod(XnbarJ) + epsilonJ*diag(1,nstrain+nstrain))/i
              #Robbins Munro tuning
              lambdaJ<- lambdaJ * exp((2/max(1, i-2000)) * (min(mh.ratioJ, 1) - 0.234))
              zigmaJ<- lambdaJ* optconstantJ * zigmaJ
            }
        }
      }
    }
    if(i %% 1000 == 0) cat("Iteration:", i, "\n")
  }
  if(Modeltype %in% (1:7) && burn_in>2000 && num_iteration>2000) print(acc/(num_iteration-2000))

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
FFBS_INFERENCE<- function(y, e_it, Modeltype, adjmat, step_sizes = list("r"=0.3,"s"=0.3,"u"=0.3),
                          num_iteration = 15000, sdBs=0.03, sdGs=0.05, sdLambdas=0.03, sdCop=0.002,
                          y_total=NULL, RM_Gs=TRUE, RM_Bs=TRUE, RM_Cop=TRUE, RM_Lambdas=TRUE,
                          burn_in=5000){
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
      sumY1<- sumY1 + y[,,k]
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
    shape1params<- rep(c(1,2), num_Gammas)
    shape2params<- rep(c(11,2),num_Gammas)
  }else if(Modeltype %in% c(2,4,6)){
    num_Gammas<- 2 * nstrain
    MC_chain<- matrix(NA, nrow=num_iteration, ncol=num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+n_copParams)
    MC_chain[1,]<- c(runif(num_Gammas), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0.1, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_factloadings), rep(0.1, n_copParams))
    shape1params<- rep(c(1,2), num_Gammas)
    shape2params<- rep(c(11,2),num_Gammas)
  }else if(Modeltype==7){
    num_Gammas<- nstate * nstate
    MC_chain<- matrix(NA, nrow=num_iteration, ncol=num_Gammas+3+time+12+ndept+nstrain+nstrain+n_factloadings+n_copParams)
    MC_chain[1,]<- c(as.numeric(t(initGs)), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0.1, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_factloadings), rep(0.1, n_copParams))
  }

  Q_r<- MC_chain[1,num_Gammas + 1] * RW2PrecMat
  Q_s<- MC_chain[1,num_Gammas + 2] * RW1PrecMat
  Q_u<- MC_chain[1,num_Gammas + 3] * R

  Qstz_r<- QRstz_basis(time)
  Qstz_s<- QRstz_basis(12)
  Qstz_u<- QRstz_basis(ndept)

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

  Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[1, num_Gammas+3+(1:time)], s=MC_chain[1, num_Gammas+3+time+(1:12)], u=MC_chain[1, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)
  likelihoodcurrent<- Allquantities$loglike
  priorcurrentRcomps<- randomwalk2_cpp(MC_chain[1, num_Gammas+3+(1:time)], MC_chain[1, num_Gammas+1])
  priorcurrentScomps<- seasonalComp2_cpp(MC_chain[1, num_Gammas+3+time+(1:12)], MC_chain[1, num_Gammas+2], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1_cpp(MC_chain[1, num_Gammas+3+time+12+(1:ndept)], MC_chain[1, num_Gammas+3], R, rankdef)

  grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
  current_rtrans<- t(Qstz_r) %*% MC_chain[1, num_Gammas+3+(1:time)]
  current_strans<- t(Qstz_s) %*% MC_chain[1, num_Gammas+3+time+(1:12)]
  current_utrans<- t(Qstz_u) %*% MC_chain[1, num_Gammas+3+time+12+(1:ndept)]

  deltaP<- 1
  RMdelta<- 1/(0.234*(1-0.234))  #for Betas, Gammas
  RMLdelta<- 1/(0.44*(1-0.44))   #for factor_loadings
  sdLambdas<- rep(sdLambdas, nstrain)
  sdLambdasJoint<- sdLambdas

  for (i in 2:num_iteration) {
    current_r <- MC_chain[i-1, num_Gammas+3+(1:time)]
    current_s <- MC_chain[i-1, num_Gammas+3+time+(1:12)]
    current_u <- MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)]

    MC_chain[i,num_Gammas+1]<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(current_r) %*% RW2PrecMat %*% current_r)/2)
    MC_chain[i,num_Gammas+2]<- rgamma(1, shape = 1 + 11/2, rate = 0.001 + (t(current_s) %*% RW1PrecMat %*% current_s)/2)
    MC_chain[i,num_Gammas+3]<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(current_u) %*% R %*% current_u)/2)

    Q_r<- MC_chain[i,num_Gammas+1] * RW2PrecMat
    Q_s<- MC_chain[i,num_Gammas+2] * RW1PrecMat
    Q_u<- MC_chain[i,num_Gammas+3] * R

    #Update s
    Mmatcs<- as.numeric(grad_current$cov_s %*% grad_current$grad_s)

    eps_s <- rnorm(11)
    proposedScompst <- as.numeric(current_strans + 0.5 * step_sizes$s^2 * Mmatcs + step_sizes$s * chol(grad_current$cov_s) %*% eps_s)

    proposedScomps<- Qstz_s %*% proposedScompst
    Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=current_r, s=proposedScomps, u=current_u, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- Allquantities$loglike
    q_prop <- mvnfast::dmvn(proposedScompst, mu = current_strans + 0.5 * step_sizes$s^2 * Mmatcs, sigma = grad_current$cov_s * step_sizes$s^2, log = TRUE)
    q_curr <- mvnfast::dmvn(t(current_strans), mu = proposedScompst + 0.5 * step_sizes$s^2 * Mmatps, sigma = grad_proposed$cov_s * step_sizes$s^2, log = TRUE)
    priorproposedScomps<- seasonalComp2_cpp(proposedScomps, MC_chain[i, num_Gammas+2], RW1PrecMat)

    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop
    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
      MC_chain[i, num_Gammas+3+time+(1:12)]<- proposedScomps
      likelihoodcurrent<- likelihoodproposed
      current_strans<- proposedScompst
      priorcurrentScomps<- priorproposedScomps
      grad_current<- grad_proposed
      Allquantities<- Allquantities
    }else{
      MC_chain[i, num_Gammas+3+time+(1:12)]<- MC_chain[i-1, num_Gammas+3+time+(1:12)]
    }

    #Update r
    eps_r <- rnorm(time-1)
    Mmatrc<- as.numeric(grad_current$cov_r %*% grad_current$grad_r)

    proposedRcompst <- as.numeric(current_rtrans + 0.5 * step_sizes$r^2 * Mmatrc + step_sizes$r * chol(grad_current$cov_r) %*% eps_r)

    proposedRcomps<- Qstz_r %*% proposedRcompst
    Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=proposedRcomps, s=MC_chain[i, num_Gammas+3+time+(1:12)], u=current_u, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1, Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcompst, mu = current_rtrans + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(t(current_rtrans), mu = proposedRcompst + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike
    priorproposedRcomps <- randomwalk2_cpp(proposedRcomps, MC_chain[i, num_Gammas+1])

    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop

    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
      MC_chain[i, num_Gammas+3+(1:time)] <- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
      current_rtrans<- proposedRcompst
      priorcurrentRcomps<- priorproposedRcomps
      grad_current<- grad_proposed
      Allquantities<- Allquantities
    }else{
      MC_chain[i, num_Gammas+3+(1:time)]<- MC_chain[i-1, num_Gammas+3+(1:time)]
    }

    #Update u
    eps_u <- rnorm(ndept-1)
    Mmatuc<- as.numeric(grad_current$cov_u %*% grad_current$grad_u)

    proposedUcompst <- as.numeric(current_utrans + 0.5 * step_sizes$u^2 * Mmatuc + step_sizes$u * chol(grad_current$cov_u) %*% eps_u)

    proposedUcomps<- Qstz_u %*% proposedUcompst
    Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=proposedUcomps, jointTPM=JointTPM, B=MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)

    Mmatup<- as.numeric(grad_proposed$cov_u %*% grad_proposed$grad_u)

    q_prop <- mvnfast::dmvn(proposedUcompst, mu = current_utrans + 0.5 * step_sizes$u^2 * Mmatuc, sigma = grad_current$cov_u * step_sizes$u^2, log = TRUE)
    q_curr <- mvnfast::dmvn(t(current_utrans), mu = proposedUcompst + 0.5 * step_sizes$u^2 * Mmatup, sigma = grad_proposed$cov_u * step_sizes$u^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike

    priorproposedUcomps<- logIGMRF1_cpp(proposedUcomps, MC_chain[i, num_Gammas+3], R, rankdef)

    log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop
    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
      MC_chain[i, num_Gammas+3+time+12+(1:ndept)]<- proposedUcomps
      current_utrans<- proposedUcompst
      likelihoodcurrent<- likelihoodproposed
      priorcurrentUcomps<- priorproposedUcomps
      grad_current<- grad_proposed
      Allquantities<- Allquantities
    }else{
      MC_chain[i, num_Gammas+3+time+12+(1:ndept)]<- MC_chain[i-1, num_Gammas+3+time+12+(1:ndept)]
    }

    if(Model == 0){
      MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)]
      MC_chain[i, 1:2]<- MC_chain[i-1, 1:2]
    }else{
      proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain)))

      priorcurrentB<- sum(dgamma(MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

      Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=proposedB, Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedB
                     - likelihoodcurrent - priorcurrentB)
      #print(paste("mh.ratioB = ", mh.ratio))

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- proposedB
        likelihoodcurrent<- likelihoodproposed
      }
      else{
        MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+(1:nstrain)]
      }
      if(RM_Bs && i<burn_in && !is.na(mh.ratio)) {sdBs= sdBs * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))}

      if(Modeltype %in% c(1,2)){
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = shape1params, shape2 = shape2params, log=TRUE))

        proposedLambdas<- rep(0, nstrain)
        JointTPM1<- Multipurpose_JointTransitionMatrix(proposedGs, nstrain, proposedLambdas, Modeltype)
        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)
        grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)

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
        if(RM_Gs && i<burn_in && !is.na(mh.ratio)) {sdGs= sdGs * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))}
      }else if(Modeltype %in% c(3,4)){
        #Transition probabilities update
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = shape1params, shape2 = shape2params, log=TRUE))

        JointTPM1<- Multipurpose_JointTransitionMatrix2(proposedGs, nstrain, MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)], Modeltype, gh)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        if(any(!is.finite(JointTPM1))){
          MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
          if(RM_Gs && i<burn_in) {sdGs= sdGs * exp((RMdelta/i) * (0 - 0.234))}
          sdGs <- max(sdGs, 1e-6)
        }else{

          Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)

          likelihoodproposed<- Allquantities$loglike

          mh.ratioGC<- exp(likelihoodproposed + priorproposedGs
                           - likelihoodcurrent - priorcurrentGs)

          #print(mh.ratioGC)

          if(!is.na(mh.ratioGC) && runif(1) < mh.ratioGC){
            MC_chain[i, 1:num_Gammas]<- proposedGs
            likelihoodcurrent<- likelihoodproposed
            JointTPM<- JointTPM1
          }
          else{
            MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
          }
          if(RM_Gs && i<burn_in && !is.na(mh.ratioGC)) {sdGs= sdGs * exp((RMdelta/i) * (min(mh.ratioGC, 1) - 0.234))}
        }

        #Factor loadings update
        LAMBDAS<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]

        for(l in 1:n_factloadings){

        LAMBDAS_prop <- LAMBDAS
        eta <- atanh(LAMBDAS[l])
        eta_prop <- rnorm(1, eta, sdLambdas[l])
        LAMBDAS_prop[l] <- tanh(eta_prop)
        LAMBDAS_current <- LAMBDAS[l]

        JointTPM1<- Multipurpose_JointTransitionMatrix2(MC_chain[i, 1:num_Gammas], nstrain, LAMBDAS_prop, Modeltype, gh)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        if(any(!is.finite(JointTPM1))){
          MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+l]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+l]
          if(RM_Lambdas && i<burn_in) {sdLambdas[l]= sdLambdas[l] * exp((RMLdelta/i) * (0 - 0.44))}
          sdLambdas[l] <- max(sdLambdas[l], 1e-6)
        }else{

          Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)

          likelihoodproposed<- Allquantities$loglike

          mh.ratioL<- exp(likelihoodproposed + log(1 - LAMBDAS_prop[l]^2)
                          - likelihoodcurrent - log(1 - LAMBDAS_current^2))

          if(!is.na(mh.ratioL) && runif(1) < mh.ratioL){
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+l]<- LAMBDAS_prop[l]
            LAMBDAS[l]<- LAMBDAS_prop[l]
            likelihoodcurrent<- likelihoodproposed
            JointTPM<- JointTPM1
          }
          else{
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+l]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+l]
            LAMBDAS[l]<- MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+nstrain+l]
          }
          if(RM_Lambdas && i<burn_in && !is.na(mh.ratioL)) {sdLambdas[l]= sdLambdas[l] * exp((RMLdelta/i) * (min(mh.ratioL, 1) - 0.44))}
            }
        }

        #Joint FactorLoadings update when posterior is less sensitive to individual update
        eta <- atanh(MC_chain[i,num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)])
        eta_prop <- rnorm(n_factloadings, eta, sdLambdasJoint)
        LAMBDAS_prop <- tanh(eta_prop)
        LAMBDAS_current <- MC_chain[i,num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]

        JointTPM1<- Multipurpose_JointTransitionMatrix2(MC_chain[i,1:num_Gammas], nstrain, LAMBDAS_prop, Modeltype, gh)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

        if(any(!is.finite(JointTPM1))){
          MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]
          if(RM_Lambdas && i<burn_in) {sdLambdasJoint= sdLambdasJoint * exp((RMdelta/i) * (0 - 0.234))}
          sdLambdasJoint<- max(sdLambdasJoint, 1e-6)
        }else{

          Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total = y_total)

          likelihoodproposed<- Allquantities$loglike

          mh.ratioGC<- exp(likelihoodproposed + sum(log(1 - LAMBDAS_prop^2))
                           - likelihoodcurrent - sum(log(1 - LAMBDAS_current^2)))

          if(!is.na(mh.ratioGC) && runif(1) < mh.ratioGC){
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- LAMBDAS_prop
            likelihoodcurrent<- likelihoodproposed
            JointTPM<- JointTPM1
          }
          else{
            MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]<- MC_chain[i, num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:n_factloadings)]
          }
          if(RM_Lambdas && i<burn_in && !is.na(mh.ratioGC)) {sdLambdasJoint= sdLambdasJoint * exp((RMdelta/i) * (min(mh.ratioGC, 1) - 0.234))}
        }

        Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)
        grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
        likelihoodcurrent<- Allquantities$loglike
      }else if(Modeltype %in% c(5,6)){
        proposedGs<- abs(rnorm(num_Gammas,mean=MC_chain[i-1,1:num_Gammas], sd=rep(sdGs, num_Gammas)))
        proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

        priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:num_Gammas], shape1 = shape1params, shape2 = shape2params, log=TRUE))
        priorproposedGs<- sum(dbeta(proposedGs, shape1 = shape1params, shape2 = shape2params, log=TRUE))

       JointTPM1<- Multipurpose_JointTransitionMatrix(proposedGs, nstrain, MC_chain[i-1, ncol(MC_chain)], Modeltype)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)
        if(any(!is.finite(JointTPM1))){
          MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
          if(RM_Gs && i<burn_in) {sdGs= sdGs * exp((RMdelta/i) * (0 - 0.234))}
          sdGs <- max(sdGs, 1e-6)
        }else{

          Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)

          likelihoodproposed<- Allquantities$loglike

            mh.ratio<- exp(likelihoodproposed + priorproposedGs
                           - likelihoodcurrent - priorcurrentGs)

          if(!is.na(mh.ratio) && runif(1) < mh.ratio){
            MC_chain[i, 1:num_Gammas]<- proposedGs
            likelihoodcurrent<- likelihoodproposed
            JointTPM<- JointTPM1
            }
          else{
            MC_chain[i, 1:num_Gammas]<- MC_chain[i-1,1:num_Gammas]
          }
          if(RM_Gs && i<burn_in && !is.na(mh.ratio)) {sdGs= sdGs * exp((RMdelta/i) * (min(mh.ratio, 1) - 0.234))}
        }

        proposedcopPs<- rnorm(1,mean=MC_chain[i-1, ncol(MC_chain)], sd=sdCop)

        if(nstrain != 2 && proposedcopPs < 0){
          MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
          if(RM_Cop && i<burn_in) {sdCop= sdCop * exp((RMLdelta/i) * (0 - 0.44))}
          sdCop <- max(sdCop, 1e-6)
        }else{

        JointTPM1<- Multipurpose_JointTransitionMatrix(MC_chain[i, 1:num_Gammas], nstrain, proposedcopPs, Modeltype)

        JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
        JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)
        if(any(!is.finite(JointTPM1))){
          MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
          if(RM_Cop && i<burn_in) {sdCop= sdCop * exp((RMLdelta/i) * (0 - 0.44))}
          sdCop <- max(sdCop, 1e-6)
        }else{

          Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)

          likelihoodproposed<- Allquantities$loglike

          mh.ratio<- exp(likelihoodproposed - likelihoodcurrent)

          if(!is.na(mh.ratio) && runif(1) < mh.ratio){
            likelihoodcurrent<- likelihoodproposed
            JointTPM<- JointTPM1
            MC_chain[i, ncol(MC_chain)]<- proposedcopPs
          }
          else{
            MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
          }
          if(RM_Cop && i<burn_in && !is.na(mh.ratio)) {sdCop= sdCop * exp((RMLdelta/i) * (min(mh.ratio, 1) - 0.44))}
          }
        }
        Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)
        grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
        likelihoodcurrent<- Allquantities$loglike
      }else if(Modeltype==7){

        for(n in 1:nstate){
          index<- nstate * (n-1) + 1

          JointTPM1<- JointTPM
          JointTPM1[n, ] <- gtools::rdirichlet(1, rep(1/nstate, nstate) + deltaP * MC_chain[i-1, (index:(n*nstate))])

          proposalproposedGs<-  log(gtools::ddirichlet(JointTPM1[n, ], MC_chain[i-1, (index:(n*nstate))]))
          proposalcurrentproposedGs<- log(gtools::ddirichlet(MC_chain[i-1, (index:(n*nstate))], JointTPM1[n, ]))

          priorcurrentGs<- log(gtools::ddirichlet(MC_chain[i-1, (index:(n*nstate))], rep(1/nstate, nstate)))
          priorproposedGs<- log(gtools::ddirichlet(JointTPM1[n, ], rep(1/nstate, nstate)))

          Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM1, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)

          likelihoodproposed<- Allquantities$loglike

          mh.ratio<- exp(likelihoodproposed + priorproposedGs + proposalcurrentproposedGs
                         - likelihoodcurrent - priorcurrentGs - proposalproposedGs)

          #print(mh.ratio)

          if(!is.na(mh.ratio) && runif(1) < mh.ratio){
            MC_chain[i, (index:(n*nstate))]<- as.numeric(JointTPM1[n, ])
            JointTPM<- JointTPM1
            likelihoodcurrent<- likelihoodproposed
            deltaP<- max(0, deltaP-3)
          }
          else{
            MC_chain[i, (index:(n*nstate))]<- MC_chain[i-1, (index:(n*nstate))]
            deltaP<- deltaP + 1
          }
        }
        JointTPM<- matrix(MC_chain[i, 1:num_Gammas], nrow = nstate, ncol = nstate, byrow = TRUE)
        Allquantities<- FFBSgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, num_Gammas+3+(1:time)], s=MC_chain[i, num_Gammas+3+time+(1:12)], u=MC_chain[i, num_Gammas+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, num_Gammas+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)
        grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s, cov_u=Allquantities$cov_u)
        likelihoodcurrent<- Allquantities$loglike
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
