functions{
  int intPower(int a, int b){
    int res = 1;
    for (i in 1:b) {
    res = a * res;
    }
    return(res);
}

    int get_bit(int x, int k) {
     return (x / (intPower(2,k-1))) % 2;
   }

  vector logVecMatMult(vector logV, matrix logM){
    int S = dims(logV)[1];
    vector[S] res;
    for(s in 1:S){
      res[s] = log_sum_exp(logV + logM[, s]);
    }
    return(res);
  }

  matrix JointTransitionMatrix_Independent(array[,,] real gamma, int K){
    int S = intPower(2, K);
    matrix[S, S] Gamma;
      for(a in 0:(S - 1)){
      for(b in 0:(S - 1)){
      real prob = 1;
      for(k in 1:K){
       int from_k = (a %/% intPower(2, (k - 1))) % 2;
       int to_k = (b %/% intPower(2, (k - 1))) % 2;
        prob = prob * gamma[from_k + 1, to_k + 1, k];
      }
      Gamma[a + 1, b + 1] = prob;
    }
  }
  return(Gamma);
  }

    //Frank copula CDF
real FrankCDF(real Psi, vector u){
  int d = num_elements(u);
  real num = 1.0;
  for(i in 1:d){
    num *= (exp(-Psi * u[i]) - 1);
  }
  real denom = pow(exp(-Psi) - 1, d - 1);
  real result = -1.0 / Psi * log1p(num / denom);

  return result;
}

   real one_factor_copula(vector u, vector lambda, vector gh_x, vector gh_w) {
    int K = num_elements(u);
    int M = num_elements(gh_x);

    vector[K] q;
    real result = 0;
    for (k in 1:K) {
      real uk = fmin(fmax(u[k], 1e-12), 1 - 1e-12);
      q[k] = inv_Phi(uk);
    }

    // Gauss–Hermite quadrature
    for (j in 1:M) {
      real z = sqrt(2.0) * gh_x[j];
      real prod_term = 1.0;

      for (k in 1:K) {
        real denom = sqrt(1 - square(lambda[k]));
        real arg = (q[k] - lambda[k] * z) / denom;
        prod_term *= Phi(arg);
      }
      result += gh_w[j] * prod_term;
    }

    return result / sqrt(pi());
  }

  matrix JointTransitionMatrix_Copula(array[,,] real gamma, int K, vector lambda, vector gh_x, vector gh_w,
                                    array[] int num_subsets, array[,] int subset_sizes, array[,] int subset_indices, int Modeltype, real frankparam){

    int S = intPower(2,K);
    matrix[S, S] Gamma;
    array[2,2,K] real gamma2 = gamma;
    for(m in 1:K){
      gamma2[1,1,m] = gamma[1,2,m];
      gamma2[1,2,m] = gamma[1,1,m];
    }

    for (a in 0:(S-1)) {
      for (b in 0:(S-1)) {
        vector[K] prob;
        array[K] int ones;
        array[K] int zeros;
        int n1 = 0;
        int n0 = 0;

        for (k in 1:K) {
          int from_k = get_bit(a, k);
          int to_k   = get_bit(b, k);
          if (from_k == 1) {
            n1 += 1;
            ones[n1] = k;
          } else {
            n0 += 1;
            zeros[n0] = k;
          }
          prob[k] = gamma2[from_k + 1, to_k + 1, k];
        }
        real total = 0;

int TRACKER = 0;

for (s in 1:num_subsets[a+1]) {
  int size_s = subset_sizes[a+1, s];
  int sign = (size_s % 2 == 0) ? 1 : -1;

  vector[K] u = rep_vector(1.0, K);

  //Indices
  for (i in 1:n1){
    u[ones[i]] = prob[ones[i]];
  }

  //complement elements
  for (j in 1:size_s) {
    int idx = subset_indices[a+1, TRACKER + j];
    if (idx > 0)
      u[idx] = prob[idx];
  }

  TRACKER += size_s;
if(Modeltype == 3 || Modeltype == 4){
  total += sign * one_factor_copula(u, lambda, gh_x, gh_w);
}else if(Modeltype == 6 || Modeltype == 6){
  total += sign * FrankCDF(frankparam, u);
}
}
        Gamma[a+1, b+1] = fmax(total, 1e-300);;
      }
    }
    for (i in 1:S) {
      real row_sum = sum(Gamma[i]);
      Gamma[i] /= row_sum;
    }
    return Gamma;
  }

  //Stationary distribution
  vector stationarydist(matrix Gamma){
    int ncol = dims(Gamma)[2];

    matrix[ncol, ncol] mT = transpose(Gamma);

    complex_vector[ncol] E_values = eigenvalues(mT); vector[ncol] NE_values = get_real(E_values);
    complex_matrix[ncol, ncol] E_vectors = eigenvectors(mT);

    real min_dist = positive_infinity();
    int index = 1;
    for(i in 1:ncol){
    real dist = abs(NE_values[i] - 1);
    if(dist < min_dist){
    min_dist = dist;
    index = i;
    }
  }
    complex_vector[ncol] stationary_distribution = E_vectors[, index];
    vector[ncol] Nstationary_distribution = get_real(stationary_distribution);
    Nstationary_distribution /= sum(Nstationary_distribution);
    return(Nstationary_distribution);
}

  //colsums function
  array[] int colSums(array[,] int M){
    int ncol = dims(M)[2];
    array[ncol] int sums;

    for(i in 1:ncol){
      sums[i] = sum(M[,i]);
    }
    return(sums);
  }

  // Intrinsic GMRF density
  real IGMRF1(vector uconstrained, real kappa_u, matrix R, int rankdef) {
    int n = rows(R);
    return (((n - rankdef) / 2.0) * (log(kappa_u) - log(2.0 * pi())) - (kappa_u / 2.0) * quad_form(R, uconstrained));
  }

  // Random walk density
  real randomwalk2(vector r, real kappa_r) {
    int time = dims(r)[1];
    real res = 0;
    for (i in 3:time) {
      res += (r[i-2] - (2 * r[i-1]) + r[i])^2;
    }
    return (((time - 2) / 2.0) * log(kappa_r) - (kappa_r / 2.0) * res);
  }

  // Seasonal components' density
    real seasonalComp(vector s, real kappa_s, matrix SMat) {
    int n = rows(SMat);
    return (((n - 1) / 2.0) * (log(kappa_s) - log(2.0 * pi())) - (kappa_s / 2.0) * quad_form(SMat, s));
}

//log Generalized Double Pareto Density
  real logGDP(vector Psi){
  int n = num_elements(Psi);
  real res1 = log(3.0) - log(2 * 1.0);
  real res2 = 0;
  for(i in 1:n){
   res2 += -(3.0 + 1.0) * log(1 + abs(Psi[i]) / 1.0);
  }
  real res = res1 + res2;
  return res;
}

  //Transition matrix function
    array[,,] real MakeTransMatrix(vector transprob, int nstrain) {
    int n_transprob = num_elements(transprob);
    array[2,2,nstrain] real res;
      res[1,1,1] = 1-transprob[1];
      res[1,2,1] = transprob[1];
      res[2,1,1] = transprob[2];
      res[2,2,1] = 1-transprob[2];
      if(n_transprob == 2){
          for(i in 2:nstrain){
            res[1,1,i] = 1-transprob[1];
            res[1,2,i] = transprob[1];
            res[2,1,i] = transprob[2];
            res[2,2,i] = 1-transprob[2];
          }
        }else if(n_transprob>2){
          int index = 1;
            for(i in 2:nstrain){
              res[1,1,i] = 1-transprob[i+index];
              res[1,2,i] = transprob[i+index];
              res[2,1,i] = transprob[i+index+1];
              res[2,2,i] = 1-transprob[i+index+1];
              index += 1;
        }
      }
    return res;
  }

  //loglikelihood via forward filtering
  real Stan_Loglikelihood(array[,,] int y, vector a_k, vector r, vector s, vector u, array[,,] real gamma, matrix e_it, vector B, int Modeltype, matrix Bits,
                          vector lambda, vector gh_x, vector gh_w, array[] int num_subsets, array[,] int subset_sizes, array[,] int subset_indices, array[,] int y_total, array[] vector mod7TransitionArray, real frankparam){
  int ndept = dims(y)[1];
  int time = dims(y)[2];
  int nstrain = dims(a_k)[1];
  int nstate = intPower(2, nstrain);

  if(Modeltype == 0){
    real allLoglikelihood = 0;
  for(i in 1:ndept) {
    for (t in 1:time) {
    int month_index = (t - 1) % 12 + 1;
    vector[nstrain] PoisMean = e_it[i, t] * exp(a_k + r[t] + s[month_index] + u[i]);
    vector[nstrain] safePoisMean = fmax(PoisMean, 1e-12);

    //y slice
    array[nstrain] int y_obs = y[i, t];

    // mask: 1 if observed, 0 if missing
    vector[nstrain] mask;

    for(k in 1:nstrain){
      if(y_obs[k] == -1){
        mask[k] = 0;
        y_obs[k] = 0;
      } else{
        mask[k] = 1;
      }
    }
    vector[nstrain] loglike = to_vector(y_obs) .* log(safePoisMean) - safePoisMean;
    allLoglikelihood += dot_product(mask, loglike);

    if(sum(mask) == 0 && y_total[i, t] != -1){
      allLoglikelihood += y_total[i, t] * sum(log(safePoisMean)) - sum(safePoisMean);
    }
  }
}
      return allLoglikelihood;
 }else{
   matrix[nstate, nstate] jointTPM;
   if(Modeltype==1){
      jointTPM = JointTransitionMatrix_Independent(gamma, nstrain);
   }else if(Modeltype==2){
      jointTPM = JointTransitionMatrix_Independent(gamma, nstrain);
   }else if(Modeltype==3){
      jointTPM = JointTransitionMatrix_Copula(gamma, nstrain, lambda, gh_x, gh_w, num_subsets, subset_sizes, subset_indices, Modeltype, frankparam);
   }else if(Modeltype==4){
      jointTPM = JointTransitionMatrix_Copula(gamma, nstrain, lambda, gh_x, gh_w, num_subsets, subset_sizes, subset_indices, Modeltype, frankparam);
   }else if(Modeltype==5){
      jointTPM = JointTransitionMatrix_Copula(gamma, nstrain, lambda, gh_x, gh_w, num_subsets, subset_sizes, subset_indices, Modeltype, frankparam);
   }else if(Modeltype==6){
      jointTPM = JointTransitionMatrix_Copula(gamma, nstrain, lambda, gh_x, gh_w, num_subsets, subset_sizes, subset_indices, Modeltype, frankparam);
   }else if(Modeltype==7){
     for(m in 1:nstate){
       jointTPM[m] = transpose(mod7TransitionArray[m]);
     }
   }
  vector[nstate] alpha_prev;
  vector[nstate] alpha_curr;
  vector[nstate] log_init_density = log(stationarydist(jointTPM));
  jointTPM = log(jointTPM);
  vector[ndept] log_forwards;

  vector[nstate] BdotBits;
  for (n in 1:nstate) {
    BdotBits[n] = dot_product(B, Bits[n, ]);
  }

  for (i in 1:ndept){
    vector[nstate] prodEmission = rep_vector(0, nstate);
  for(n in 1:nstate){
  for(k in 1:nstrain){
  if(y[i, 1, k] != -1){
    prodEmission[n] += poisson_lpmf(y[i, 1, k] | e_it[i, 1] * exp(a_k[k] + r[1] + s[1] + u[i] + BdotBits[n]));
    }
  }
   if(sum(y[i, 1]) == -nstrain && y_total[i, 1] != -1){
       prodEmission[n] += y_total[i, 1] * sum(log(e_it[i, 1] * exp(a_k + r[1] + s[1] + u[i] + BdotBits[n]))) - sum(e_it[i, 1] * exp(a_k + r[1] + s[1] + u[i] + BdotBits[n])) - lgamma(y_total[i, 1] + 1);
    }
}
   alpha_curr = log_init_density + prodEmission;

      for(t in 2:time) {
        alpha_prev = alpha_curr;
        prodEmission = rep_vector(0, nstate);
        int month_index = (t - 1) % 12 + 1;
        real PoisMean_it = e_it[i, t] * exp(r[t] + s[month_index] + u[i]);
        for(n in 1:nstate){
        for(k in 1:nstrain){
        if(y[i, t, k] != -1){
          prodEmission[n] += poisson_lpmf(y[i, t, k] | PoisMean_it * exp(a_k[k] + BdotBits[n]));
         }
        }
        if(sum(y[i, t]) == -nstrain && y_total[i, t] != -1){
           prodEmission[n] += y_total[i, t] * sum(log(PoisMean_it .* exp(a_k + BdotBits[n]))) - sum(PoisMean_it .* exp(a_k + BdotBits[n])) - lgamma(y_total[i, t] + 1);
        }
      }
        alpha_curr = logVecMatMult(alpha_prev, jointTPM) + prodEmission;
    }
        log_forwards[i] = log_sum_exp(alpha_curr);
  }
      real fullLogLikelihood = sum(log_forwards);
      return fullLogLikelihood;
    }
   return 0;
    }
}

data {
  int<lower=1> ndept;                 // Number of departments
  int<lower=1> time;                  // Time
  int<lower=1> nstate;                // Number of states
  int<lower=1> rankdef;               // Rank deficiency of structure matrix (R)
  int<lower=1> nstrain;               // Number of strains
  int<lower=0> npar;
  array[ndept, time, nstrain] int y;  // data matrix
  matrix[ndept, time] e_it;           // initial Susceptibles
  matrix[ndept, ndept] R;             // Structure matrix (IGMRF1)
  int<lower=0> Modeltype;             // Model's functional form
  matrix[12, 12] SMat;                //Structure matrix (seasonal_comp)
  matrix[nstate, nstrain] Bits;       //Bits matrix
  vector[30] gh_x;
  vector[30] gh_w;
  int<lower=1> n_subsets;
  array[n_subsets] int num_subsets;
  array[nstate,n_subsets] int subset_sizes;
  array[nstate, 100] int subset_indices;
  array[ndept, time] int y_total;
  matrix[nstate, nstate] DirichletPrior;
  int mod7nstate;
}

parameters {
  vector<lower=0,upper=1>[2*npar] transitionParams;            // transition to hyperendemic
  real<lower=0> kappa_u;                // spatial precision parameter
  real<lower=0> kappa_r;                 // trend precision parameter
  real<lower=0> kappa_s;                 // seasonal precision parameter
  vector[ndept-1] u;                     // Spatial components
  vector[time] rraw;                     // Trend components
  vector[11] sraw;                       // cyclic seasonal components
  vector<lower=0>[npar] B;               // autoregressive parameters
  vector[nstrain] a_k;                   // background intercepts
  vector<lower=-4,upper=4>[nstrain-1] lambda_free;
  real frankparam;
  array[mod7nstate] simplex[mod7nstate] mod7TransitionArray;
//  row_stochastic_matrix[mod7nstate, mod7nstate] mod7TransitionMatrix;
}

transformed parameters {
  real sumC = sum(u[1:(ndept-1)]);
  vector[ndept] uconstrained;
  uconstrained = append_row(u[1:(ndept-1)], -sumC);
  real sumS = sum(sraw[1:11]);
  vector[12] s;
  s = append_row(sraw[1:11], -sumS);

  vector[time] r;
  real MeanR = mean(rraw);
  for(i in 1:time){
  r[i] = rraw[i] - MeanR;
  }
  vector[nstrain] phi_k;
  phi_k = exp(a_k[1:nstrain]);
  vector<lower=-0.9999,upper=0.9999>[nstrain] lambda;
  lambda = append_row(0.9999, 0.99*tanh(lambda_free));
}

model {
  target += gamma_lpdf(kappa_u | 1, 0.01);
  target += gamma_lpdf(kappa_r | 1, 0.0001);
  target += gamma_lpdf(kappa_s | 1, 0.001);
  target += gamma_lpdf(phi_k | 0.1, 0.01/exp(-15));
  if(Modeltype>0 && Modeltype !=7){
    if(Modeltype == 3){
        target += beta_lpdf(transitionParams[1] | 1, 11);
        target += beta_lpdf(transitionParams[2] | 6, 6);
        target += logGDP(lambda_free);
    }else if(Modeltype == 4){
      int index = 0;
      for(i in 1:nstrain){
        target += beta_lpdf(transitionParams[i+index] | 1, 11);
        target += beta_lpdf(transitionParams[i+index+1] | 6, 6);
        index += 1;
        }
        target += logGDP(lambda_free);
    }
      target += gamma_lpdf(B | 2, 2);
      if(nstrain==2){
      target += normal_lpdf(frankparam | 0, 100);
      }else{
        target += exponential_lpdf(frankparam | 0.5);
      }
  }
  if(Modeltype==7){
    for(n in 1:nstate){
      target += dirichlet_lpdf(mod7TransitionArray[n]|DirichletPrior[n,]);
    }
  }
  target += IGMRF1(uconstrained, kappa_u, R, rankdef);
  target += randomwalk2(r, kappa_r);
  target += seasonalComp(s, kappa_s, SMat);

  // Likelihood
  target += Stan_Loglikelihood(y, a_k, r, s, uconstrained, MakeTransMatrix(transitionParams,nstrain), e_it, B, Modeltype, Bits, lambda, gh_x, gh_w, num_subsets, subset_sizes, subset_indices, y_total, mod7TransitionArray, frankparam);
}

generated quantities{
}
