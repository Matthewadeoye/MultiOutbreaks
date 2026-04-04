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

  matrix JointTransitionMatrix_sameTPM(matrix gamma, int K){
    int S = intPower(2, K);
    matrix[S, S] Gamma;
      for(a in 0:(S - 1)){
      for(b in 0:(S - 1)){
      real prob = 1;
      for(k in 1:K){
       int from_k = (a %/% intPower(2, (k - 1))) % 2;
       int to_k = (b %/% intPower(2, (k - 1))) % 2;
        prob = prob * gamma[from_k + 1, to_k + 1];
      }
      Gamma[a + 1, b + 1] = prob;
    }
  }
  return(Gamma);
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

  matrix FactorJointTransitionMatrix_copula(matrix gamma, int K, vector lambda, vector gh_x, vector gh_w,
                                    array[] int num_subsets, array[,] int subset_sizes, array[,] int subset_indices){

    int S = intPower(2,K);
    matrix[S, S] Gamma;
    matrix[2,2] gamma2 = gamma;
    gamma2[1,1] = gamma[1,2];
    gamma2[1,2] = gamma[1,1];

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
          prob[k] = gamma2[from_k + 1, to_k + 1];
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

  total += sign * one_factor_copula(u, lambda, gh_x, gh_w);
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

  //Gamma matrix function
    matrix G(real G12, real G21) {
    matrix[2, 2] m;
    m[1, 1] = 1 - G12;
    m[1, 2] = G12;
    m[2, 1] = G21;
    m[2, 2] = 1 - G21;
    return m;
  }

  //loglikelihood via forward filtering
  real Stan_Loglikelihood(array[,,] int y, vector a_k, vector r, vector s, vector u, matrix gamma, matrix e_it, vector B, int Modeltype, matrix Bits,
                          vector lambda, vector gh_x, vector gh_w, array[] int num_subsets, array[,] int subset_sizes, array[,] int subset_indices){
  int ndept = dims(y)[1];
  int time = dims(y)[2];
  int nstrain = dims(a_k)[1];
    //Model0
  if(Modeltype == 0){
    real allLoglikelihood = 0;
  for(i in 1:ndept) {
    for(t in 1:time) {
      for(k in 1:nstrain){
      if(y[i, t, k] == -1){
        allLoglikelihood += 0;
      }else{
      int month_index = (t - 1) % 12 + 1;
        allLoglikelihood += poisson_lpmf(y[i, t, k] | e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i]));
        }
      }
    }
  }
      return allLoglikelihood;
}else if(Modeltype == 1){
  int nstate = intPower(2, nstrain);
  matrix[nstate, nstate] jointTPM = JointTransitionMatrix_sameTPM(gamma, nstrain);
  matrix[time, nstate] Alpha;
  vector[nstate] init_density = stationarydist(jointTPM);
  vector[ndept] log_forwards;

  for (i in 1:ndept){
    vector[nstate] prodEmission;
    for(m in 1:nstate){
      prodEmission[m] = 0;
    }
  // Initialization of the first time step for each department
  for(n in 1:nstate){
  for(k in 1:nstrain){
  if(y[i, 1, k] == -1){
    prodEmission[n] += 0;
    }else{
    prodEmission[n] += poisson_lpmf(y[i, 1, k] | e_it[i, 1] * exp(a_k[k] + r[1] + s[1] + u[i] + dot_product(B, Bits[n, ])));
    }
  }
}
  Alpha[1] = transpose(log(init_density) + prodEmission);

  // Dynamic programming loop for the remaining time steps
      for(t in 2:time) {
        int month_index = (t - 1) % 12 + 1;
        for(n in 1:nstate){
        for(k in 1:nstrain){
        if(y[i, t, k] == -1){
          prodEmission[n] += 0;
          }else{
         prodEmission[n] += poisson_lpmf(y[i, t, k] | e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + dot_product(B, Bits[n, ])));
         }
        }
      }
        Alpha[t] = transpose(logVecMatMult(transpose(Alpha[t-1, ]), log(jointTPM)) + prodEmission);
    }
        log_forwards[i] = log_sum_exp(Alpha[time, ]);
  }
      real fullLogLikelihood = sum(log_forwards);
      return fullLogLikelihood;
    }else if(Modeltype == 2){
  int nstate = intPower(2, nstrain);

  matrix[nstate, nstate] jointTPM = FactorJointTransitionMatrix_copula(gamma, nstrain, lambda, gh_x, gh_w, num_subsets, subset_sizes, subset_indices);
  matrix[time, nstate] Alpha;
  vector[nstate] init_density = stationarydist(jointTPM);
  vector[ndept] log_forwards;

  for (i in 1:ndept){
    vector[nstate] prodEmission;
    for(m in 1:nstate){
      prodEmission[m] = 0;
    }
  for(n in 1:nstate){
  for(k in 1:nstrain){
  if(y[i, 1, k] == -1){
    prodEmission[n] += 0;
    }else{
    prodEmission[n] += poisson_lpmf(y[i, 1, k] | e_it[i, 1] * exp(a_k[k] + r[1] + s[1] + u[i] + dot_product(B, Bits[n, ])));
    }
  }
}
  Alpha[1] = transpose(log(init_density) + prodEmission);

      for(t in 2:time) {
        int month_index = (t - 1) % 12 + 1;
        for(n in 1:nstate){
        for(k in 1:nstrain){
        if(y[i, t, k] == -1){
          prodEmission[n] += 0;
          }else{
         prodEmission[n] += poisson_lpmf(y[i, t, k] | e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + dot_product(B, Bits[n, ])));
         }
        }
      }
        Alpha[t] = transpose(logVecMatMult(transpose(Alpha[t-1, ]), log(jointTPM)) + prodEmission);
    }
        log_forwards[i] = log_sum_exp(Alpha[time, ]);
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
  int<lower=0> Modeltype;        // Model's functional form
  matrix[12, 12] SMat;                //Structure matrix (seasonal_comp)
  matrix[nstate, nstrain] Bits;       //Bits matrix
  vector[30] gh_x;
  vector[30] gh_w;
  int<lower=1> n_subsets;
  array[n_subsets] int num_subsets;
  array[nstate,n_subsets] int subset_sizes;
  array[nstate, 100] int subset_indices;
}

parameters {
  real<lower=0, upper=1> G12;            // transition to hyperendemic
  real<lower=0, upper=1> G21;            // transition to endemic
  real<lower=0> kappa_u;                 // spatial precision parameter
  real<lower=0> kappa_r;                 // trend precision parameter
  real<lower=0> kappa_s;                 // seasonal precision parameter
  vector[ndept-1] u;                     // Spatial components
  vector[time] rraw;                     // Trend components
  vector[11] sraw;                       // cyclic seasonal components
  vector<lower=0>[npar] B;            // autoregressive parameters
  vector[nstrain] a_k;                   // background intercepts
  vector<lower=-4,upper=4>[nstrain-1] lambda_free;
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
  lambda = append_row(0.99, 0.99*tanh(lambda_free));
}

model {
  target += beta_lpdf(G12 | 2, 2);
  target += beta_lpdf(G12 | 2, 2);
  target += gamma_lpdf(kappa_u | 1, 0.01);
  target += gamma_lpdf(kappa_r | 1, 0.0001);
  target += gamma_lpdf(kappa_s | 1, 0.001);
  target += gamma_lpdf(phi_k | 0.1, 0.01/exp(-15));
  if(Modeltype>0){
  target += gamma_lpdf(B | 2, 2);
  target += logGDP(lambda_free);
  }
  target += IGMRF1(uconstrained, kappa_u, R, rankdef);
  target += randomwalk2(r, kappa_r);
  target += seasonalComp(s, kappa_s, SMat);

  // Likelihood
  target += Stan_Loglikelihood(y, a_k, r, s, uconstrained, G(G12, G21), e_it, B, Modeltype, Bits, lambda, gh_x, gh_w, num_subsets, subset_sizes, subset_indices);
}

generated quantities{
}
