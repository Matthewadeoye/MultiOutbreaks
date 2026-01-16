encodeBits<- function(K){
  S <- 2^K
  bits<- matrix(0, nrow = S, ncol = K)
  for (i in 0:(S - 1)) {
    for (k in 1:K) {
      bits[i + 1, k] <- (i %/% 2^(k - 1)) %% 2
    }
  }
  bits <- bits[, K:1]
  return(bits)
}

G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

Multipurpose_JointTransitionMatrix<- function(gammas, K, copParams, Modeltype){
  if(Modeltype==1){
    JointTPM<- JointTransitionMatrix_arma_cpp(G(gammas[1], gammas[2]), K)
  }else if(Modeltype==2){
    Glist<- BuildGamma_list_cpp(gammas)
    JointTPM<- JointTransitionMatrix_per_strain_cpp2(Glist, K)
  }else if(Modeltype == 3){
    JointTPM<- JointTransitionMatrix_Gaussiancopula_cpp(G(gammas[1], gammas[2]), K, copParams)
  }else if(Modeltype == 4){
    Glist<- BuildGamma_list_cpp(gammas)
    JointTPM<- JointTransitionMatrix_Gaussiancopula_perstrain_cpp(Glist, K, copParams)
  }else if(Modeltype==5){
    JointTPM<- JointTransitionMatrix_Frankcopula_cpp(G(gammas[1], gammas[2]), K, copParams)
  }else if(Modeltype==6){
    Glist<- BuildGamma_list_cpp(gammas)
    JointTPM<- JointTransitionMatrix_Frankcopula_perstrain_cpp(Glist, K, copParams)
  }else if(Modeltype == 7){
    nstate<- 2^K
    JointTPM<- matrix(gammas, nrow = nstate, ncol = nstate, byrow = TRUE)
  }
  return(JointTPM)
}

Multipurpose_JointTransitionMatrix2<- function(gammas, K, Lambdas, Modeltype, gh){
  nstate<- 2^K

  if(Modeltype==3){
    JointTPM<- JointTransitionMatrix_1FactorGaussiancopula_cpp(G(gammas[1], gammas[2]), K, Lambdas, gh$nodes, gh$weights)
  }else if(Modeltype==4){
    Glist<- BuildGamma_list_cpp(gammas)
    JointTPM<- JointTransitionMatrix_1FactorGaussiancopula_per_strain_cpp(Glist, K, Lambdas, gh$nodes, gh$weights)
  }
  return(JointTPM)
}

Posteriormultstrain.Decoding<- function(y, e_it, inf.object, Modeltype, thinningL=1000, burn.in=1000){
  ndept<- dim(y)[1]
  time<- dim(y)[2]
  nstrain<- dim(y)[3]
  nstate<- 2^nstrain
  Bits<- encodeBits(K=nstrain)

  fullG.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "G")]
  fullr.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "r")]
  fulls.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "s")]
  fullu.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "u")]
  fullB.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "B")]
  fulla_k.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "a")]
  if(Modeltype %in% c(3,4,5,6)) fullcop.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "c")]

  thinning<- numeric(floor(nrow(fullr.draws)/thinningL))
  thinning[1]<- thinningL
  for(i in 2:length(thinning)){
    thinning[i]<- thinning[i-1] + thinningL
  }

  G.draws<- fullG.draws[thinning, ]
  r.draws<- fullr.draws[thinning, ]
  s.draws<- fulls.draws[thinning, ]
  u.draws<- fullu.draws[thinning, ]
  B.draws<- fullB.draws[thinning, ]
  a_k.draws<- fulla_k.draws[thinning, ]
  if(Modeltype %in% c(3,4)){
    cop.draws<- fullcop.draws[thinning, ]
  }else if(Modeltype %in% c(5,6)){
    cop.draws<- fullcop.draws[thinning]
  }

  sum_P_itn<- array(0, dim = c(ndept, time, nstate))
  for(index in 1:length(thinning)){
    Gs<- as.numeric(G.draws[index,])
    r<- as.numeric(r.draws[index,])
    s<- as.numeric(s.draws[index,])
    u<- as.numeric(u.draws[index,])
    B<- as.numeric(B.draws[index,])
    a_k<- as.numeric(a_k.draws[index,])
    if(Modeltype %in% c(3,4)){
      cop<- as.numeric(cop.draws[index, ])
    }else if(Modeltype %in% c(5,6)){
      cop<- cop.draws[index]
    }else{
      cop<- 0
    }

    if(Modeltype == 1){
      JointTPM<- Multipurpose_JointTransitionMatrix(Gs, nstrain, cop, Modeltype)
    }else if(Modeltype == 2){
      JointTPM<- Multipurpose_JointTransitionMatrix(Gs, nstrain, cop, Modeltype)
    }else if(Modeltype == 3){
      JointTPM<- Multipurpose_JointTransitionMatrix(Gs, nstrain, cop, Modeltype)
      JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
      JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
    }else if(Modeltype == 4){
      JointTPM<- Multipurpose_JointTransitionMatrix(Gs, nstrain, cop, Modeltype)
      JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
      JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
    }else if(Modeltype == 5){
      JointTPM<- Multipurpose_JointTransitionMatrix(Gs, nstrain, cop, Modeltype)
      JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
      JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
    }else if(Modeltype == 6){
      JointTPM<- Multipurpose_JointTransitionMatrix(Gs, nstrain, cop, Modeltype)
      JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
      JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
    }else if(Modeltype == 7){
      JointTPM<- Multipurpose_JointTransitionMatrix(Gs, nstrain, cop, Modeltype)
    }

    P_itn <- PostOutbreakProbs_cpp(y = y, e_it = e_it, nstrain=nstrain, r = r, s = s, u = u, jointTPM = JointTPM, B = B, Bits = Bits, a_k = a_k)
    sum_P_itn<- sum_P_itn + P_itn
  }
  sum_P_itn<- sum_P_itn/length(thinning)

  perStrainProbs <- array(0, dim = c(ndept, time, nstrain))
  for(n in 1:nstate){
    for(k in 1:nstrain){
      perStrainProbs[,,k] <- perStrainProbs[,,k] + Bits[n, k] * sum_P_itn[,,n]
    }
  }
  return(perStrainProbs)
}


#Crude estimates
crudeEst<- function(y, e_it){
  ydot<- colSums(y, na.rm = T)
  edot<- colSums(e_it, na.rm = T)
  logydot<- log(ydot)
  logedot<- log(edot)
  lambdadot<- logydot-logedot
  nlambdadot<- lambdadot[lambdadot != -Inf]
  x<- 1:ncol(y)
  lambdadot<- ifelse(lambdadot== -Inf, mean(nlambdadot), lambdadot)
  success <- tryCatch({
    loess_fit <- loess(lambdadot ~ x, span = 0.5)
    TRUE
  }, error = function(e) {
    FALSE
  })

  if(success){
    loess_fit <- loess(lambdadot ~ x, span = 0.5)
    smoothed <- predict(loess_fit)
    crudeS<- lambdadot - smoothed
    crudeR<- smoothed
    crudeU<- log(rowSums(y/e_it, na.rm = T)/sum(exp(crudeR+crudeS)))
    crudeU[crudeU==-Inf]<- mean(crudeU[is.finite(crudeU)])
    crudeU<- crudeU-mean(crudeU[is.finite(crudeU)])
    #crudeU<- log(rowSums(y/e_it, na.rm = T)/sum(exp(crudeR+crudeS)))-mean(log(rowSums(y/e_it, na.rm = T)/sum(exp(crudeR+crudeS))))
    #crudeU<- rep(0, nrow(y))
  }else{
    crudeR<- rep(mean(lambdadot), ncol(y))
    crudeS<- lambdadot-mean(lambdadot)
    crudeU<- rep(0, nrow(y))
  }
  return(list(crudeR, crudeS, crudeU))
}
