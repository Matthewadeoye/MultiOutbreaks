sim_adjmat <- matrix(0, nrow = 9, ncol = 9)
uppertriang <- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
gdata::upperTriangle(sim_adjmat, byrow=TRUE) <- uppertriang
gdata::lowerTriangle(sim_adjmat, byrow=FALSE) <- uppertriang

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

Posteriormultstrain.Decoding<- function(y, e_it, inf.object, Modeltype, y_total=NULL,
                                        thinningL=1000, burn.in=1000){
  ndept<- dim(y)[1]
  time<- dim(y)[2]
  nstrain<- dim(y)[3]
  nstate<- 2^nstrain
  Bits<- encodeBits(K=nstrain)

  if(is.null(y_total)){
    sumY1<- y[,,1]
    for(k in 2:nstrain){
      sumY1<- sumY1 + y[,,1]
    }
    y_total<- sumY1
  }else{
    y_total<- y_total
  }

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

    P_itn <- PostOutbreakProbs_cpp(y = y, e_it = e_it, nstrain=nstrain, r = r, s = s, u = u, jointTPM = JointTPM, B = B, Bits = Bits, a_k = a_k, y_total = y_total)
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


Decode_to_Perstrain <- function(StatesMatrix, nstrain){
  Bits <- encodeBits(K = nstrain)
  ndept <- nrow(StatesMatrix)
  time <- ncol(StatesMatrix)

  Outbreaks <- array(NA, dim = c(ndept, time, nstrain))

  for (i in 1:ndept) {
    for (t in 1:time) {
      s <- StatesMatrix[i, t]   #joint state index
      Outbreaks[i, t, ] <- Bits[s + 1, ]  #per-strain indicators
    }
  }
  return(Outbreaks)
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

#Extract posterior credible interval
posterior_interval_custom <- function(posterior_samples, prob = 0.95) {

  lower <- (1 - prob) / 2
  upper <- 1 - lower

  # Apply quantile function across columns
  credible_intervals <- apply(posterior_samples, 2, quantile, probs = c(lower, upper))

  credible_intervals_df <- as.data.frame(t(credible_intervals))
  colnames(credible_intervals_df) <- c("2.5%", "97.5%")

  return(credible_intervals_df)
}

#custom legend
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

#Cholesky confirmation for using mvnfast
check_cholesky <- function(matrix) {
  result <- tryCatch({
    chol(matrix)
    TRUE
  }, error = function(e) {
    FALSE
  })
  return(result)
}

#QR-decomposition for sum-to-zero basis
QRstz_basis<- function(dimension){
  A<- diag(dimension)
  A[, 1]<- 1

  qrDecomposition<- qr(A)
  Q<- qr.Q(qrDecomposition)
  return(Q[, 2:dimension])
}

#log Generalized Double Pareto Density
log_GDP<- function(Psi, Alpha, Beta){
  res1<- log(Alpha) - log(2 * Beta)
  res2<- -(Alpha + 1) * log(1 + abs(Psi) / Beta)
  res<- res1 + res2
  return(res)
}

#log Dirichlet Density
log_Dirichlet<- function(Xvector, Alpha){
  res1<- sum((Alpha-1)*log(Xvector))
  res2<- sum(lgamma(Alpha)) - lgamma(sum(Alpha))
  res<- res1 - res2
  return(res)
}

#Datasets for application (Meningococcal)
#popn2010<- read.csv("C:/Users/Matthew Adeoye/Downloads/popn.csv")
#names(popn2010)<- NULL
#popn2010<- popn2010[ , -1]
#popn2010<- as.matrix(popn2010)
#popn2010<- DetectOutbreaks:::LinearInterp(popn2010)
#dim(popn2010)
#cleanedMultitypeData<- read.csv("C:/Users/Matthew Adeoye/Downloads/newdata5strains.csv")
#cleanedMultitypeData$Value<- as.numeric(cleanedMultitypeData$Value)
#ApplicationCounts<- read.csv("C:/Users/Matthew Adeoye/Downloads/updatedcases2.csv")
#ApplicationCounts$NumValue<- as.numeric(ApplicationCounts$NumValue)
#ApplicationCounts<- matrix(ApplicationCounts$NumValue, nrow = 28, byrow = T)
#MultitypeData<- array(NA, dim=c(28,120,4))
#dates= cleanedMultitypeData$Time
#countries <- sort(unique(cleanedMultitypeData$RegionName), decreasing = FALSE)
#dates=unique(dates)
#cleanedMultitypeData<- read.csv("C:/Users/Matthew Adeoye/Downloads/serogroups.csv")
#cleanedMultitypeData$Value<- as.numeric(cleanedMultitypeData$Value)
#strains<- c("NEIMENI_B", "NEIMENI_W", "NEIMENI_Y", "NEIMENI_C")

#for(t in 1:120){
#  for(i in 1:28){
#    for(k in 1:4){
#      w=which(cleanedMultitypeData$Time==dates[t] & cleanedMultitypeData$RegionName==countries[i] & cleanedMultitypeData$Category==strains[k])
#      if(length(w)==1) MultitypeData[i, t, k]<- cleanedMultitypeData$Value[w]
#    }
#  }
#}
#MultitypeData<- MultitypeData/100

#for(i in 1:28){
#  for(t in 1:120){
#    if(is.finite(ApplicationCounts[i, t]) && ApplicationCounts[i, t] == 0){
#      MultitypeData[i,t,] <- 0
#    }
#    for(k in 1:4){
#      MultitypeData[i,t,k]<- MultitypeData[i,t,k] * ApplicationCounts[i, t]
#    }
#  }
#}
#MultitypeData<- round(MultitypeData)


# FactorJointTransitionMatrix_copulaR<- function(gamma, K, copulaParams){
#
#   S <- 2^K
#   gamma[1, ] <- gamma[1, ][2:1]
#
#   Gamma <- matrix(0, S, S)
#
#   for(a in 0:(S-1)) {
#     for(b in 0:(S-1)) {
#       Indices <- integer(0)
#       IndicesComplement <- integer(0)
#       prob <- numeric(K)
#
#       for(k in 1:K){
#         from_k <- (a %/% 2^(k-1)) %% 2
#         to_k   <- (b %/% 2^(k-1)) %% 2
#
#         if(from_k == 1) Indices <- c(Indices,k)
#         else IndicesComplement <- c(IndicesComplement,k)
#
#         prob[k] <- gamma[from_k+1, to_k+1]
#       }
#
#       subsets <- sets::set_power(sets::as.set(IndicesComplement))
#       subsets <- lapply(subsets, function(g) unlist(as.vector(g)))
#
#       total <- 0
#       for(Tset in subsets){
#         sign <- (-1)^length(Tset)
#         idx <- c(Indices, Tset)
#         u <- rep(1, K)
#         if(length(idx)>0) u[idx] <- prob[idx]
#         total <- total + sign * one_factor_copula_cdf2(u, copulaParams)
#       }
#       Gamma[a+1, b+1] <- total
#     }
#   }
#   Gamma <- Gamma / rowSums(Gamma)
#   Gamma
# }
#

# matrix FactorJointTransitionMatrix_copula(matrix gamma, int K, vector lambda, vector gh_x, vector gh_w,
#                                           array[] int num_subsets, array[,] int subset_sizes, array[,] int subset_indices){
#
#   int S = intPower(2,K);
#   matrix[S, S] Gamma;
#   matrix[2,2] gamma2 = gamma;
#   gamma2[1,1] = gamma[1,2];
#   gamma2[1,2] = gamma[1,1];
#
#   for (a in 0:(S-1)) {
#     for (b in 0:(S-1)) {
#       vector[K] prob;
#       array[K] int ones;
#       array[K] int zeros;
#       int n1 = 0;
#       int n0 = 0;
#
#       for (k in 1:K) {
#         int from_k = get_bit(a, k);
#         int to_k   = get_bit(b, k);
#         if (from_k == 1) {
#           n1 += 1;
#           ones[n1] = k;
#         } else {
#           n0 += 1;
#           zeros[n0] = k;
#         }
#         prob[k] = gamma2[from_k + 1, to_k + 1];
#       }
#       real total = 0;
#
#       int offset = 0;
#
#       for (s in 1:num_subsets[a+1]) {
#         int size_s = subset_sizes[a+1, s];
#         int sign = (size_s % 2 == 0) ? 1 : -1;
#
#         vector[K] u = rep_vector(1.0, K);
#
#         //Indices
#         for (i in 1:n1){
#           u[ones[i]] = prob[ones[i]];
#         }
#
#         //complement elements
#         for (j in 1:size_s) {
#           int idx = subset_indices[a+1, offset + j];
#           if (idx > 0)
#             u[idx] = prob[idx];
#         }
#
#         offset += size_s;
#
#         total += sign * one_factor_copula_cdf(u, lambda, gh_x, gh_w);
#       }
#       Gamma[a+1, b+1] = fmax(total, 1e-300);;
#     }
#   }
#   for (i in 1:S) {
#     real row_sum = sum(Gamma[i]);
#     Gamma[i] /= row_sum;
#   }
#   return Gamma;
# }
#

PrecomputeCopulaInclusionExclusion<- function(K){
S<- 2^K
subset_indices<- matrix(NA, nrow = S, ncol = 101)
subset_sizes<- matrix(NA, nrow=S, ncol = 101)
num_subsets<- numeric(S)
for(a in 1:S) {

  Indices <- integer(0)
  IndicesComplement <- integer(0)

  for(k in 1:K){
    from_k <- ((a-1) %/% 2^(k-1)) %% 2

    if(from_k == 1) Indices <- c(Indices,k)
    else IndicesComplement <- c(IndicesComplement,k)
  }

  subsets <- sets::set_power(sets::as.set(IndicesComplement))
  subsets <- lapply(subsets, function(g) unlist(as.vector(g)))

  num_subsets[a] <- length(subsets)

  subset_lengths <- sapply(subsets, length)

  flattened <- unlist(lapply(subsets, function(sub) {
    if(length(sub) == 0) 0 else sub
  }))

  subset_indices[a, 1:length(flattened)] <- flattened
  subset_sizes[a, 1:length(subset_lengths)] <- subset_lengths
}
maxLength<- max(num_subsets)
subset_sizes<- subset_sizes[,1:maxLength]
subset_indices<-subset_indices[,-1]
subset_sizes<- ifelse(is.na(subset_sizes),-1,subset_sizes)
subset_indices<- ifelse(is.na(subset_indices),-1,subset_indices)
 allobjects<- list("n_subsets"=maxLength, "num_subsets"=num_subsets, "subset_sizes"=subset_sizes, "subset_indices"=subset_indices)
 return(allobjects)
}

# subset_indices_list <- lapply(
#   seq_len(nrow(subset_indices)),
#   function(i) as.integer(subset_indices[i, ])
# )
# subset_sizes_list <- lapply(
#   seq_len(nrow(subset_sizes)),
#   function(i) as.integer(subset_sizes[i, ])
# )
#
# FactorJointTransitionMatrix_copula(
#   G(0.1,0.2),
#   5,
#   c(0.9, 0.4, 0.7,0.1,0.2),
#   gh$nodes,
#   gh$weights,
#   num_subsets,
#   subset_sizes_list,
#   subset_indices_list
# )
#
# FactorJointTransitionMatrix_copulaR(G(0.1,0.2),5,c(0.9,0.4,0.7,0.1,0.2))
