modelevidenceLP <- function(theta, dataset) {

  kappaR <- theta["kappa_r"]
  kappaS <- theta["kappa_s"]
  kappaU <- theta["kappa_u"]
  Q_r<- kappaR * dataset[["RW2PrecMat"]]
  Q_s<- kappaS * dataset[["RW1PrecMat"]]
  Q_u<- kappaU * dataset[["R"]]

  r   <- theta[startsWith(names(theta),"r")]
  s   <- theta[startsWith(names(theta), "s")]
  u   <- theta[startsWith(names(theta), "u")]
  B   <- theta[startsWith(names(theta), "B")]
  a_k <- theta[startsWith(names(theta), "a")]

  OutbreakPrior_ExpectationMatrix<- matrix(c(11/12,1/12,6/12,6/12), nrow = 2, byrow = T)
  Dirichlet_Prior<- 12 * JointTransitionMatrix_arma_cpp(OutbreakPrior_ExpectationMatrix, dataset$nstrain)

  #gammas
  if(dataset$Modeltype %in% c(1,2,3,4,5,6,7)){
    Gammas <- theta[startsWith(names(theta), "G")]
  }

  log_prior <- 0

  if(dataset$Modeltype==0){
    B<- rep(0, dataset$nstrain)
    jointTPM<- matrix(0, dataset$nstate, dataset$nstate)
  }else if(dataset$Modeltype %in% c(1,2,5,6)){
    copulaParam <- theta[startsWith(names(theta), "c")]
    JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, dataset$nstrain, copulaParam, dataset$Modeltype)
    if(dataset$Modeltype %in% c(5,6) && dataset$nstrain==2){
      log_prior <- log_prior + dnorm(copulaParam, mean = 0, sd=100, log = TRUE)
    }else if(dataset$Modeltype %in% c(5,6) && dataset$nstrain>2){
      log_prior <- log_prior + dexp(copulaParam, rate=0.5, log = TRUE)
    }
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
    jointTPM<- JointTPM
  }else if(dataset$Modeltype %in% c(3,4)){
    copulaParam <- theta[startsWith(names(theta), "F")]
    JointTPM<- Multipurpose_JointTransitionMatrix2(Gammas, dataset$nstrain, copulaParam, dataset$Modeltype, dataset$gh)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
    jointTPM<- JointTPM
  }else if(dataset$Modeltype ==7){
    copulaParam <- theta[startsWith(names(theta), "c")]
    JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, dataset$nstrain, copulaParam, dataset$Modeltype)
    jointTPM <- pmax(JointTPM, 1e-8)
    jointTPM <- jointTPM / rowSums(jointTPM)
    jointTPM <- jointTPM
  }

  log_lik <- SMOOTHINGgradmultstrainLoglikelihood_cpp(y = dataset$y,e_it = dataset$e_it,nstrain = dataset$nstrain,r = r,s = s,u = u,jointTPM = jointTPM,B = B,Bits = dataset$Bits,a_k = a_k,Model = dataset$Modeltype,Q_r = Q_r,Q_s = Q_s,Q_u = Q_u,gradients = 0,Qstz_r = dataset$Qstz_r,Qstz_s = dataset$Qstz_s,Qstz_u = dataset$Qstz_u,y_total = dataset$y_total)$loglike

  log_prior <- log_prior + dgamma(kappaR,1,0.0001,log=TRUE) + dgamma(kappaS,1,0.001, log=TRUE) + dgamma(kappaU,1,0.01,  log=TRUE)

  log_prior <- log_prior + randomwalk2_cpp(r, kappaR) + seasonalComp2_cpp(s, kappaS, dataset$SMat) +logIGMRF1_cpp(u, kappaU, dataset$R, dataset$rankdef)

  if(dataset$Modeltype>0){
    log_prior <- log_prior + sum(dgamma(B, shape=2, rate=2, log=TRUE))
  }

  # a_k prior (with Jacobian)
    log_prior <- log_prior + sum(dgamma(exp(a_k), shape=0.01, rate=0.01/exp(-15),log=TRUE) + a_k)

  if(dataset$Modeltype %in% c(1,2,3,4,5,6)){
    log_prior <- log_prior + sum(dbeta(Gammas,dataset$shape1params,dataset$shape2params,log=TRUE))
  }else if(dataset$Modeltype==7){
    for(n in 1:dataset$nstate){
    log_prior <- log_prior + log(gtools::ddirichlet(jointTPM[n,], Dirichlet_Prior[n, ]))
    }
  }

  return(log_lik + log_prior)
}


ModelEvidenceBridgeSamplingPackage <- function(y, e_it, adjmat, Modeltype, inf.object, y_total=NULL, cores=1){
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
  SMat<- RW1PrecMat

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  Qstz_r<- QRstz_basis(time)
  Qstz_s<- QRstz_basis(12)
  Qstz_u<- QRstz_basis(ndept)

  if(is.null(y_total)){
    sumY1<- y[,,1]
    for(k in 2:nstrain){
      sumY1<- sumY1 + y[,,k]
    }
    y_total<- sumY1
  }else{
    y_total<- y_total
  }

  if(Modeltype %in% c(0,1,3,5,7)){
    num_Gammas<- 2
    shape1params<- rep(c(1,6), num_Gammas)
    shape2params<- rep(c(11, 6),num_Gammas)
  }else if(Modeltype %in% c(2,4,6)){
    num_Gammas<- 2 * nstrain
    shape1params<- rep(c(1,6), num_Gammas)
    shape2params<- rep(c(11, 6),num_Gammas)
  }

  dataset<- list(
    y=y,
    e_it=e_it,
    Modeltype=Modeltype,
    nstrain=nstrain,
    nstate=nstate,
    Bits=Bits,
    RW1PrecMat=RW1PrecMat,
    RW2PrecMat=RW2PrecMat,
    R=R,
    Qstz_r=Qstz_r,
    Qstz_s=Qstz_s,
    Qstz_u=Qstz_u,
    y_total=y_total,
    SMat=SMat,
    rankdef=rankdef,
    shape1params=shape1params,
    shape2params=shape2params,
    gh=gh)

  if(Modeltype==0){
    inf.object<- inf.object[, !(startsWith(colnames(inf.object), "G") |
                                                 startsWith(colnames(inf.object), "B"))]
  }else if(Modeltype %in% c(3, 4)){
    inf.object<- inf.object[, !(startsWith(colnames(inf.object), "c"))]
  }

  Model<- ifelse(Modeltype>0,1,0)
  paramNames<- colnames(inf.object)
  samples <- as.matrix(inf.object)
  colnames(samples)<- paramNames

  if(is.null(colnames(samples))){
    stop("Posterior samples must have column names.")
  }

  lb <- rep(-Inf, ncol(samples))
  ub <- rep(Inf, ncol(samples))
  names(lb) <- paramNames
  names(ub) <- paramNames

  # positivity constraints
  lb[startsWith(names(lb), "k")]<- 0
if(Modeltype %in% c(1,2,7)){
  lb[startsWith(names(lb), "B")]<- 0
  lb[startsWith(names(lb), "G")]<- 0
  ub[startsWith(names(ub), "G")]<- 1
}else if(Modeltype %in% c(3,4)){
  lb[startsWith(names(lb), "B")]<- 0
  lb[startsWith(names(lb), "G")]<- 0
  ub[startsWith(names(ub), "G")]<- 1
  lb[startsWith(names(lb), "F")]<- -1
  ub[startsWith(names(ub), "F")]<- 1
}else if(Modeltype %in% c(5,6) && nstrain>2){
  lb[startsWith(names(lb), "B")]<- 0
  lb[startsWith(names(lb), "G")]<- 0
  ub[startsWith(names(ub), "G")]<- 1
  lb[startsWith(names(lb), "c")]<- 0
}else if(Modeltype %in% c(5,6)){
  lb[startsWith(names(lb), "B")]<- 0
  lb[startsWith(names(lb), "G")]<- 0
  ub[startsWith(names(ub), "G")]<- 1
}

  bridge_result <- bridgesampling::bridge_sampler(samples = samples,log_posterior = modelevidenceLP, data = dataset, lb = lb, ub = ub, cores = cores)

  return(bridge_result)
}


#Model evidence via Importance sampling method
ModelEvidence<- function(y, e_it, adjmat, Modeltype, inf.object, num_samples = 5000, y_total=NULL){

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  ndept<- dim(y)[1]
  time<- dim(y)[2]
  nstrain<- dim(y)[3]
  nstate<- 2^nstrain
  Bits<- encodeBits(nstrain)
  gh <- statmod::gauss.quad(30, kind = "hermite")

  if(is.null(y_total)){
    sumY1<- y[,,1]
    for(k in 2:nstrain){
      sumY1<- sumY1 + y[,,k]
    }
    y_total<- sumY1
  }else{
    y_total<- y_total
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

  SMat<- matrix(0, nrow=12, ncol=12)
  SMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  SMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  SMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  SMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  SMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    SMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  Qstz_r<- QRstz_basis(time)
  Qstz_s<- QRstz_basis(12)
  Qstz_u<- QRstz_basis(ndept)

  if(Modeltype==0){inf.object<- inf.object[, !(startsWith(colnames(inf.object), "G") |
                                                 startsWith(colnames(inf.object), "B"))]}
  Model<- ifelse(Modeltype>0,1,0)

  mu<- colMeans(inf.object)
  varcov<- cov(inf.object)
  eps <- 1e-6
  varcov <- varcov + diag(eps, nrow(varcov))
#  varcov <- as.matrix(Matrix::nearPD(varcov)$mat)
  cond<- check_cholesky(varcov)
  a<- numeric(num_samples)

  if(Modeltype==0){
    num_Gammas<- 0
    JointTPM<- matrix(0, nstate, nstate)
  }else if(Modeltype %in% c(1,3,5)){
    num_Gammas<- 2
  }else if(Modeltype %in% c(2,4,6)){
    num_Gammas<- 2 * nstrain
  }else if(Modeltype==7){
    num_Gammas<- nstate * nstate
  }
  shape1params<- rep(c(1, 2), num_Gammas)
  shape2params<- rep(c(11, 2), num_Gammas)

  if(cond){
    theta <- matrix(nrow = 1, ncol = length(mu))
    class(theta) <- "numeric"
    for(i in 1:num_samples){
      mvnfast::rmvt(n=1, mu = mu, sigma = varcov, df = 3, A = theta)
      if(Modeltype==0){
        if(any(theta[1:3]<=0)){
          a[i]<- NA
        }else{
          kappaR<- theta[1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[2]
          Q_s <- kappaS * SMat
          kappaU<- theta[3]
          Q_u<- kappaU * R
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          B<- 0
          intercepts<- theta[3+time+12+ndept+(1:nstrain)]
          a[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)$loglike +
            sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
            randomwalk2_cpp(r, kappaR) +
            seasonalComp2_cpp(s, kappaS, SMat) +
            logIGMRF1_cpp(u, kappaU, R, rankdef) -
            mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(5,6) && nstrain>2){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0) || theta[length(theta)]<0){
          a[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, nstrain, theta[length(theta)], Modeltype)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          a[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)$loglike +
                    sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                    sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                    randomwalk2_cpp(r, kappaR) +
                    seasonalComp2_cpp(s, kappaS, SMat) +
                    logIGMRF1_cpp(u, kappaU, R, rankdef) +
                    sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE)) +
                    ifelse(Modeltype %in% c(5,6), dnorm(theta[length(theta)], mean = 0, sd=100, log = TRUE), 0)
                  - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(3,4)){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0) || any(abs(theta[num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)])>1)){
          a[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix2(Gammas, nstrain, theta[num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)], Modeltype, gh)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          a[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)$loglike +
            sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
            sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
            randomwalk2_cpp(r, kappaR) +
            seasonalComp2_cpp(s, kappaS, SMat) +
            logIGMRF1_cpp(u, kappaU, R, rankdef) +
            sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
          - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(1,2,5,6,7)){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0)){
          a[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, nstrain, theta[length(theta)], Modeltype)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          a[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)$loglike +
                    sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                    sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                    randomwalk2_cpp(r, kappaR) +
                    seasonalComp2_cpp(s, kappaS, SMat) +
                    logIGMRF1_cpp(u, kappaU, R, rankdef) +
                    sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE)) +
                    ifelse(Modeltype %in% c(5,6), dexp(copulaParam, rate=0.5, log = TRUE), 0)
                  - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }
    }
  }else{
    for(i in 1:num_samples){
      theta<- mvtnorm::rmvt(n=1, delta = mu, sigma = varcov, df = 3)
      if(Modeltype==0){
        if(any(theta[1:3]<=0)){
          a[i]<- NA
        }else{
          kappaR<- theta[1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[2]
          Q_s <- kappaS * SMat
          kappaU<- theta[3]
          Q_u<- kappaU * R
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          B<- 0
          intercepts<- theta[3+time+12+ndept+(1:nstrain)]
          a[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)$loglike +
                    sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                    randomwalk2_cpp(r, kappaR) +
                    seasonalComp2_cpp(s, kappaS, SMat) +
                    logIGMRF1_cpp(u, kappaU, R, rankdef) -
                    mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(5,6) && nstrain>2){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0) || theta[length(theta)]<0){
          a[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, nstrain, theta[length(theta)], Modeltype)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          a[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)$loglike +
                    sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                    sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                    randomwalk2_cpp(r, kappaR) +
                    seasonalComp2_cpp(s, kappaS, SMat) +
                    logIGMRF1_cpp(u, kappaU, R, rankdef) +
                    sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
                  - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(3,4)){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0) || any(abs(theta[num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)])>1)){
          a[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix2(Gammas, nstrain, theta[num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)], Modeltype, gh)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          a[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)$loglike +
                    sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                    sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                    randomwalk2_cpp(r, kappaR) +
                    seasonalComp2_cpp(s, kappaS, SMat) +
                    logIGMRF1_cpp(u, kappaU, R, rankdef) +
                    sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
                  - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(1,2,5,6,7)){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0)){
          a[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, nstrain, theta[length(theta)], Modeltype)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          a[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u, y_total=y_total)$loglike +
                    sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                    sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                    randomwalk2_cpp(r, kappaR) +
                    seasonalComp2_cpp(s, kappaS, SMat) +
                    logIGMRF1_cpp(u, kappaU, R, rankdef) +
                    sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
                  - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }
    }
  }
  validDensities <- a[!is.na(a)]
  MarginalLikelihood<- log(1) - log(length(validDensities)) + logSumExp_cpp2(validDensities)
  print(paste(length(validDensities), "valid samples used"))
  return(MarginalLikelihood)
}


#Model evidence via Bridge sampling method
ModelEvidence_Bridge<- function(y, e_it, adjmat, Modeltype,inf.object, num_samples = 50000, y_total=NULL){
  logPy<- ModelEvidence(y=y, e_it = e_it, adjmat = adjmat, Modeltype = Modeltype,inf.object = inf.object, num_samples = num_samples, y_total = y_total)
  print(logPy)

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  ndept<- dim(y)[1]
  time<- dim(y)[2]
  nstrain<- dim(y)[3]
  nstate<- 2^nstrain
  Bits<- encodeBits(nstrain)
  gh <- statmod::gauss.quad(30, kind = "hermite")

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  SMat<- matrix(0, nrow=12, ncol=12)
  SMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  SMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  SMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  SMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  SMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    SMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
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

  Qstz_r<- QRstz_basis(time)
  Qstz_s<- QRstz_basis(12)
  Qstz_u<- QRstz_basis(ndept)

  if(Modeltype==0){inf.object<- inf.object[, !(startsWith(colnames(inf.object), "G") |
                                                 startsWith(colnames(inf.object), "B"))]}
  Model<- ifelse(Modeltype>0,1,0)

  post.s<- nrow(inf.object)/2
  s1<- post.s/(num_samples+post.s)
  s2<- num_samples/(num_samples+post.s)

  mu<- colMeans(inf.object[1:post.s, ])
  varcov<- cov(inf.object[1:post.s, ])
  eps <- 1e-6
  varcov <- varcov + diag(eps, nrow(varcov))
  cond<- check_cholesky(varcov)

  logDensProposal<- numeric(num_samples)
  logDensPosterior<- numeric(post.s)

  if(Modeltype==0){
    num_Gammas<- 0
    JointTPM<- matrix(0, nstate, nstate)
  }else if(Modeltype %in% c(1,3,5)){
    num_Gammas<- 2
  }else if(Modeltype %in% c(2,4,6)){
    num_Gammas<- 2 * nstrain
  }else if(Modeltype==7){
    num_Gammas<- nstate * nstate
  }
  shape1params<- rep(c(1, 2), num_Gammas)
  shape2params<- rep(c(11, 2), num_Gammas)

  if(cond){
    theta <- matrix(nrow = 1, ncol = length(mu))
    class(theta) <- "numeric"
    for(i in 1:num_samples){
      mvnfast::rmvt(n=1, mu = mu, sigma = varcov, df = 3, A = theta)
      if(Model==0){
        if(any(theta[1:3]<=0)){
          logDensProposal[i]<- NA
        }else{
          kappaR<- theta[1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[2]
          Q_s <- kappaS * SMat
          kappaU<- theta[3]
          Q_u<- kappaU * R
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          B<- 0
          intercepts<- theta[3+time+12+ndept+(1:nstrain)]
          logDensProposal[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                    sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                    randomwalk2_cpp(r, kappaR) +
                    seasonalComp2_cpp(s, kappaS, SMat) +
                    logIGMRF1_cpp(u, kappaU, R, rankdef) -
                    mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(5,6) && nstrain>2){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0) || theta[length(theta)]<0){
          logDensProposal[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, nstrain, theta[length(theta)], Modeltype)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          logDensProposal[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                    sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                    sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                    randomwalk2_cpp(r, kappaR) +
                    seasonalComp2_cpp(s, kappaS, SMat) +
                    logIGMRF1_cpp(u, kappaU, R, rankdef) +
                    sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
                  - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(3,4)){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0) || any(abs(theta[num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)])>1)){
          logDensProposal[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix2(Gammas, nstrain, theta[num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)], Modeltype, gh)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          logDensProposal[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                                  sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                                  sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                                  randomwalk2_cpp(r, kappaR) +
                                  seasonalComp2_cpp(s, kappaS, SMat) +
                                  logIGMRF1_cpp(u, kappaU, R, rankdef) +
                                  sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
                                - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(1,2,5,6,7)){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0)){
          logDensProposal[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, nstrain, theta[length(theta)], Modeltype)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          logDensProposal[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                                  sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                                  sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                                  randomwalk2_cpp(r, kappaR) +
                                  seasonalComp2_cpp(s, kappaS, SMat) +
                                  logIGMRF1_cpp(u, kappaU, R, rankdef) +
                                  sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
                                - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }
    }
    for(j in (post.s+1):nrow(inf.object)){
      theta<- as.numeric(inf.object[j, ])
      if(Model==0){
        kappaR<- theta[1]
        Q_r<- kappaR * RW2PrecMat
        kappaS<- theta[2]
        Q_s <- kappaS * SMat
        kappaU<- theta[3]
        Q_u<- kappaU * R
        r<- theta[3+(1:time)]
        s<- theta[3+time+(1:12)]
        u<- theta[3+time+12+(1:ndept)]
        B<- 0
        intercepts<- theta[3+time+12+ndept+(1:nstrain)]
        logDensPosterior[j-post.s]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                                        sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                                        randomwalk2_cpp(r, kappaR) +
                                        seasonalComp2_cpp(s, kappaS, SMat) +
                                        logIGMRF1_cpp(u, kappaU, R, rankdef) +
                                      - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE))
      }else{
        kappaR<- theta[num_Gammas+1]
        Q_r<- kappaR * RW2PrecMat
        kappaS<- theta[num_Gammas+2]
        Q_s <- kappaS * SMat
        kappaU<- theta[num_Gammas+3]
        Q_u<- kappaU * R
        Gammas<- theta[1:num_Gammas]
        r<- theta[num_Gammas+3+(1:time)]
        s<- theta[num_Gammas+3+time+(1:12)]
        u<- theta[num_Gammas+3+time+12+(1:ndept)]
        B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
        intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

        if(Modeltype %in% c(1,2,5,6,7)){
          JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, nstrain, theta[length(theta)], Modeltype)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
        }else if(Modeltype %in% c(3,4)){
          JointTPM<- Multipurpose_JointTransitionMatrix2(Gammas, nstrain, theta[num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)], Modeltype, gh)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
        }
        logDensPosterior[j-post.s]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                                        sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                                        sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                                        randomwalk2_cpp(r, kappaR) +
                                        seasonalComp2_cpp(s, kappaS, SMat) +
                                        logIGMRF1_cpp(u, kappaU, R, rankdef) +
                                        sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
                                      - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE))
      }
    }
  }else{
    for(i in 1:num_samples){
      theta<- mvtnorm::rmvt(n=1, delta = mu, sigma = varcov, df = 3)
      if(Model==0){
        if(any(theta[1:3]<=0)){
          logDensProposal[i]<- NA
        }else{
          kappaR<- theta[1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[2]
          Q_s <- kappaS * SMat
          kappaU<- theta[3]
          Q_u<- kappaU * R
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          B<- 0
          intercepts<- theta[3+time+12+ndept+(1:nstrain)]
          logDensProposal[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                                  sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                                  randomwalk2_cpp(r, kappaR) +
                                  seasonalComp2_cpp(s, kappaS, SMat) +
                                  logIGMRF1_cpp(u, kappaU, R, rankdef) -
                                  mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(5,6) && nstrain>2){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0) || theta[length(theta)]<0){
          logDensProposal[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, nstrain, theta[length(theta)], Modeltype)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          logDensProposal[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                                  sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                                  sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                                  randomwalk2_cpp(r, kappaR) +
                                  seasonalComp2_cpp(s, kappaS, SMat) +
                                  logIGMRF1_cpp(u, kappaU, R, rankdef) +
                                  sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
                                - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(3,4)){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0) || any(abs(theta[num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)])>1)){
          logDensProposal[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix2(Gammas, nstrain, theta[num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)], Modeltype, gh)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          logDensProposal[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                                  sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                                  sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                                  randomwalk2_cpp(r, kappaR) +
                                  seasonalComp2_cpp(s, kappaS, SMat) +
                                  logIGMRF1_cpp(u, kappaU, R, rankdef) +
                                  sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
                                - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }else if(Modeltype %in% c(1,2,5,6,7)){
        if(any(theta[1:num_Gammas] <= 0) || any(theta[1:num_Gammas] >= 1) || any(theta[num_Gammas+(1:3)] <= 0)){
          logDensProposal[i]<- NA
        }else{
          kappaR<- theta[num_Gammas+1]
          Q_r<- kappaR * RW2PrecMat
          kappaS<- theta[num_Gammas+2]
          Q_s <- kappaS * SMat
          kappaU<- theta[num_Gammas+3]
          Q_u<- kappaU * R
          Gammas<- theta[1:num_Gammas]
          r<- theta[num_Gammas+3+(1:time)]
          s<- theta[num_Gammas+3+time+(1:12)]
          u<- theta[num_Gammas+3+time+12+(1:ndept)]
          B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
          intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

          JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, nstrain, theta[length(theta)], Modeltype)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

          logDensProposal[i]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                                  sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                                  sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                                  randomwalk2_cpp(r, kappaR) +
                                  seasonalComp2_cpp(s, kappaS, SMat) +
                                  logIGMRF1_cpp(u, kappaU, R, rankdef) +
                                  sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
                                - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE))
        }
      }
    }
    for(j in (post.s+1):nrow(inf.object)){
      theta<- as.numeric(inf.object[j, ])
      if(Model==0){
        kappaR<- theta[1]
        Q_r<- kappaR * RW2PrecMat
        kappaS<- theta[2]
        Q_s <- kappaS * SMat
        kappaU<- theta[3]
        Q_u<- kappaU * R
        r<- theta[3+(1:time)]
        s<- theta[3+time+(1:12)]
        u<- theta[3+time+12+(1:ndept)]
        B<- 0
        intercepts<- theta[3+time+12+ndept+(1:nstrain)]
        logDensPosterior[j-post.s]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                                        sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                                        randomwalk2_cpp(r, kappaR) +
                                        seasonalComp2_cpp(s, kappaS, SMat) +
                                        logIGMRF1_cpp(u, kappaU, R, rankdef) +
                                        - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE))
      }else{
        kappaR<- theta[num_Gammas+1]
        Q_r<- kappaR * RW2PrecMat
        kappaS<- theta[num_Gammas+2]
        Q_s <- kappaS * SMat
        kappaU<- theta[num_Gammas+3]
        Q_u<- kappaU * R
        Gammas<- theta[1:num_Gammas]
        r<- theta[num_Gammas+3+(1:time)]
        s<- theta[num_Gammas+3+time+(1:12)]
        u<- theta[num_Gammas+3+time+12+(1:ndept)]
        B<- theta[num_Gammas+3+time+12+ndept+(1:nstrain)]
        intercepts<- theta[num_Gammas+3+time+12+ndept+nstrain+(1:nstrain)]

        if(Modeltype %in% c(1,2,5,6,7)){
          JointTPM<- Multipurpose_JointTransitionMatrix(Gammas, nstrain, theta[length(theta)], Modeltype)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
        }else if(Modeltype %in% c(3,4)){
          JointTPM<- Multipurpose_JointTransitionMatrix2(Gammas, nstrain, theta[num_Gammas+3+time+12+ndept+nstrain+nstrain+(1:nstrain)], Modeltype, gh)
          JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
          JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
        }
        logDensPosterior[j-post.s]<- (SMOOTHINGgradmultstrainLoglikelihood_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=r, s=s, u=u, jointTPM=JointTPM, B=B, Bits=Bits, a_k=intercepts, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u,gradients=0,Qstz_r=Qstz_r, Qstz_s=Qstz_s, Qstz_u=Qstz_u,y_total=y_total)$loglike +
                                        sum(dbeta(Gammas, shape1 = shape1params, shape2 = shape2params, log = TRUE)) +
                                        sum(dgamma(c(kappaR, kappaS, kappaU), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                                        randomwalk2_cpp(r, kappaR) +
                                        seasonalComp2_cpp(s, kappaS, SMat) +
                                        logIGMRF1_cpp(u, kappaU, R, rankdef) +
                                        sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
                                      - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE))
      }
    }
  }

  log_w<- logDensProposal[!is.na(logDensProposal)] - logSumExp_cpp2(logDensProposal[!is.na(logDensProposal)])
  propESS<- exp(-logSumExp_cpp2(2 * log_w))
  postESS<- exp((2*logSumExp_cpp2(logDensPosterior))-(logSumExp_cpp2(2*logDensPosterior)))
#  print(paste("ESS from proposal samples is ", propESS))
#  print(paste("ESS from posterior samples is ", postESS))

  numerator<- numeric(num_samples)
  denominator<- numeric(post.s)
  delta<- 0.03

  while(delta > 0.01){
    for(i in 1:num_samples){
      numerator[i]<- logDensProposal[i] - logSumExp_cpp2(c(log(s1)+logDensProposal[i], log(s2) + logPy))
    }
    for(j in (post.s+1):nrow(inf.object)){
      denominator[j-post.s]<- log(1) - logSumExp_cpp2(c(log(s1)+logDensPosterior[j-post.s], log(s2) + logPy))
    }
    numerator<- numerator[!is.na(numerator)]
    denominator<- denominator[!is.na(denominator)]
    currentlogPy<- (-log(length(numerator))+logSumExp_cpp2(numerator))-(-log(post.s)+logSumExp_cpp2(denominator))
    delta<- abs(logPy-currentlogPy)
    print(currentlogPy)
    logPy<- currentlogPy
  }
  return(currentlogPy)
}
