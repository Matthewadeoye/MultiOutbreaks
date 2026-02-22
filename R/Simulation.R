sim.RW2mean0<- function(time, sd=0.0001, init.r1 = 0, init.r2 = 0){
  r <- numeric(time)
  r[1] <- init.r1
  r[2] <- init.r2
  for(t in 3:time){
    epsilon <- rnorm(1, mean = 0, sd = sd)
    r[t] <- 2*(r[t - 1]) - r[t - 2] + epsilon
  }
  r <- r - mean(r)
  return(r)
}

sim.Seasonals<- function(Amplitude, Cycle = 1:12){
  frequency <- 1/length(Cycle)
  s <- sin(2 * pi * frequency * Cycle)
  s <- Amplitude * s
  return(s)
}

sim.GMRF <- function(n, Q, tol=1e-12){
  Es <- eigen(Q)
  LAMBDAS <- rev(Es$values)
  ai <- LAMBDAS > tol
  Vtilde <- Es$vectors[,sum(ai):1]

  yij <- matrix(
    stats::rnorm(sum(ai) * n, 0, sqrt(rep(LAMBDAS[ai], n)^-1)),
    nrow=sum(ai),
    ncol=n)

  t(Vtilde %*% yij)
}

sim.Spatials<- function(adj.matrix, sd = 0.2){
  Q<- -1 * adj.matrix
  diag(Q)<- -rowSums(Q, na.rm = T)
  Q<- (1/sd^2)*Q
  u<- c(sim.GMRF(1,Q))
  return(u)
}

simulateMarkovChain<- function(nstep, TPM){
  state.space<- 0:(nrow(TPM)-1)
  initial_state <- 0
  states <- numeric(nstep)
  states[1] <- initial_state
  for(i in 2:nstep) {
    current_state <- states[i - 1]
    next_state <- sample(state.space, size = 1, prob = TPM[current_state + 1, ])
    states[i] <- next_state
  }
  return(states)
}

simulateMultiModel<- function(Modeltype, time, nstrain, adj.matrix, copulaParam,
                              e_it=matrix(c(rep(c(rpois(time, 500000), rpois(time, 1000000)), 4), rpois(time, 500000)),
                                          byrow = T, ncol = time),
                              B = runif(nstrain), T.prob = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = T),
                              r = sim.RW2mean0(time, sd=0.009), s = sim.Seasonals(Amplitude = 1.4),
                              u = sim.Spatials(adj.matrix)){
  ndept<- nrow(adj.matrix)
  y_itk<- array(NA, dim=c(ndept, time, nstrain))
  EpidemicIndicator<- matrix(NA, ndept, time)
  Jointstates<- 2^nstrain
  Bits<- encodeBits(K=nstrain)
  aVec<- numeric(nstrain)
  for(k in 1:nstrain){
    a_k<- runif(1, min = -14, max = -12)
    aVec[k]<- a_k
  }
  #Due to lazyloading
  r<- r;  s<- s;  u<- u;  e_it<- e_it;  B<- B

  Model<- ifelse(Modeltype>0,1,0)

  if(Modeltype == 1){
    copulaParam<- 0
    T.prob<- c(T.prob[1,2],T.prob[2,1])
    JointTPM<- Multipurpose_JointTransitionMatrix(T.prob, nstrain, copulaParam, Modeltype)
  }else if(Modeltype == 2){
    copulaParam<- 0
    T.prob<- runif(2*nstrain, min = 0.1, max = 0.2)
    JointTPM<- Multipurpose_JointTransitionMatrix(T.prob, nstrain, copulaParam, Modeltype)
  }else if(Modeltype == 3){
    T.prob<- c(T.prob[1,2],T.prob[2,1])
    JointTPM<- Multipurpose_JointTransitionMatrix(T.prob, nstrain, copulaParam, Modeltype)
  }else if(Modeltype == 4){
    T.prob<- runif(2*nstrain, min = 0.1, max = 0.2)
    JointTPM<- Multipurpose_JointTransitionMatrix(T.prob, nstrain, copulaParam, Modeltype)
  }else if(Modeltype == 5){
    T.prob<- c(T.prob[1,2],T.prob[2,1])
    JointTPM<- Multipurpose_JointTransitionMatrix(T.prob, nstrain, copulaParam, Modeltype)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
  }else if(Modeltype == 6){
    T.prob<- runif(2*nstrain, min = 0.1, max = 0.2)
    JointTPM<- Multipurpose_JointTransitionMatrix(T.prob, nstrain, copulaParam, Modeltype)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
  }else if(Modeltype == 7){
    copulaParam<- 0
    JointTPM<- matrix(NA, nrow = Jointstates, ncol = Jointstates)
    T.prob<- 0
    JointTPM<- gtools::rdirichlet(Jointstates, sample(2:7, size = Jointstates, replace = T))
  }

  if(Model == 0){
    for(i in 1:ndept){
      for(t in 1:time){
        m<- (t - 1) %% 12 + 1
        for(k in 1:nstrain){
          lograte <- aVec[k] + r[t] + s[m] + u[i]
          y_itk[i, t, k]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
        }
      }
    }
    return(list("y" =y_itk, "e_it"=e_it, "r"=r, "s"=s, "u"=u, "states"=EpidemicIndicator, "a_k"=aVec))
  }
  else{
    for(i in 1:ndept){
      for(k in 1:nstrain){
        lograte<- aVec[k] + r[1] + s[1] + u[i]
        y_itk[i, 1, k]<- rpois(1, lambda = e_it[i, 1] * exp(lograte))
      }
      EpidemicIndicator[i, ]<- simulateMarkovChain(nstep=time, JointTPM)
    }

    for(t in 2:time){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        for(k in 1:nstrain){
          newB<- rep(0, nstrain)
          newB[k]<- B[k]
          lograte<- aVec[k] + r[t] + s[m] + u[i] + (newB %*% Bits[EpidemicIndicator[i, t]+1, ])
          y_itk[i, t, k]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
        }
      }
    }
    return(list("y" =y_itk, "e_it"=e_it, "r"=r, "s"=s, "u"=u, "states"=EpidemicIndicator, "B"=B, "a_k"=aVec, "T.prob"=T.prob, "JointTPM"=JointTPM, "copulaParam"=copulaParam))
  }
}

