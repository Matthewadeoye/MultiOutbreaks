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
  shape1params<- rep(c(2, 2), num_Gammas)
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
                    sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                    randomwalk2_cpp(r, kappaR) +
                    seasonalComp2_cpp(s, kappaS, SMat) +
                    logIGMRF1_cpp(u, kappaU, R, rankdef) +
                    sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
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
            sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                    sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
                    randomwalk2_cpp(r, kappaR) +
                    seasonalComp2_cpp(s, kappaS, SMat) +
                    logIGMRF1_cpp(u, kappaU, R, rankdef) +
                    sum(dgamma(B, shape = rep(2, length(B)), rate = rep(2, length(B)), log = TRUE))
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
                    sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                    sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                    sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
  shape1params<- rep(c(2, 2), num_Gammas)
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
                    sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                                  sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                                  sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                                        sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                                        sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                                  sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                                  sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                                  sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                                        sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
                                        sum(dgamma(c(kappaR, kappaS), shape = c(1, 1, 1), rate = c(0.0001, 0.001, 0.01), log = TRUE)) +
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
