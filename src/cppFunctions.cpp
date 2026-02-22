#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
#include <omp.h>
#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
int intPower(int a, int b){
  int res = 1;
  for (int i = 0; i < b; ++i) {
    res = a * res;
  }
  return res;
}

// [[Rcpp::export]]
double logSumExp_cpp2(arma::vec x) {
  double max_val = max(x);
  return log(sum(exp(x - max_val))) + max_val;
}

// [[Rcpp::export]]
arma::vec logVecMatMult2(arma::vec logV, arma::mat logM) {
  int S = logV.size();
  arma::vec res(S);
  for (int p = 0; p < S; ++p) {
    arma::vec temp(S);
    for (int s = 0; s < S; ++s) {
      temp[s] = logV[s] + logM(s, p);
    }
    res[p] = logSumExp_cpp2(temp);
  }
  return res;
}

// [[Rcpp::export]]
arma::vec stationarydistArma_cpp(arma::mat Gamma) {
  int n = Gamma.n_cols;
  Eigen::MatrixXd m(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      m(i, j) = Gamma(i, j);
    }
  }

  Eigen::MatrixXd mT = m.transpose();

  // Eigen decomposition
  Eigen::EigenSolver<Eigen::MatrixXd> es(mT);
  Eigen::VectorXd eigenvalues = es.eigenvalues().real();
  Eigen::MatrixXd eigenvectors = es.eigenvectors().real();

  // Index of eigenvalue close to 1
  int index = 0;
  double min_diff = std::abs(eigenvalues(0) - 1.0);
  for (int i = 1; i < eigenvalues.size(); ++i) {
    double diff = std::abs(eigenvalues(i) - 1.0);
    if (diff < min_diff) {
      min_diff = diff;
      index = i;
    }
  }

  // corresponding eigenvector
  Eigen::VectorXd stationary_distribution = eigenvectors.col(index);

  // Normalize stationary distribution
  stationary_distribution = stationary_distribution / stationary_distribution.sum();

  // Convert result to NumericVector
  arma::vec result(stationary_distribution.data(),
                   stationary_distribution.size(),
                   /* copy_aux_mem = */ true);
  return result;
}

// [[Rcpp::export]]
double randomwalk2_cpp(arma::vec componentR, double PrecisionR){
  int time = componentR.size();
  double Sumres = 0;
  for(int i = 2; i < time; ++i) {
    double res = pow(componentR[i - 2] - (2 * componentR[i - 1]) + componentR[i], 2);
    Sumres += res;
  }
  return (time - 2) / 2.0 * log(PrecisionR) - PrecisionR / 2.0 * Sumres;
}

// [[Rcpp::export]]
double randomwalk2_sumzero_cpp(const arma::vec& x, double PrecisionR){
  arma::uword n = x.n_elem;
  double sumR = arma::accu(x);
  arma::vec x_new = arma::join_vert(x.head(n - 1), arma::vec(1).fill(-sumR));

  double Sumres = 0.0;
  for(int i = 2; i < n; ++i) {
    double res = pow(x_new[i - 2] - (2 * x_new[i - 1]) + x_new[i], 2);
    Sumres += res;
  }

  return (n - 2) / 2.0 * log(PrecisionR) - PrecisionR / 2.0 * Sumres;
}

// [[Rcpp::export]]
double seasonalComp2_cpp(const arma::vec& x, double y, const arma::mat& z) {
  arma::uword n = x.n_elem;
  double sumC = arma::accu(x);
  arma::vec x_new = arma::join_vert(x.head(n-1), arma::vec(1).fill(-sumC));

  double result = ((n - 1) / 2.0 ) * (std::log(y) - std::log(2 * M_PI)) - 0.5 * y * arma::as_scalar(x_new.t() * z * x_new);

  return result;
}

// [[Rcpp::export]]
double logIGMRF1_cpp(const arma::vec& x, double y, const arma::mat& z, int rankdef) {
  arma::uword n = x.n_elem;
  double sumC = arma::accu(x);
  arma::vec x_new = arma::join_vert(x.head(n-1), arma::vec(1).fill(-sumC));

  double result = ((n - rankdef) / 2.0 ) * (std::log(y) - std::log(2 * M_PI)) - 0.5 * y * arma::as_scalar(x_new.t() * z * x_new);

  return result;
}

//[[Rcpp::export]]
arma::mat JointTransitionMatrix_arma_cpp(arma::mat gamma, int K) {
  int S = intPower(2, K);
  arma::mat jointGamma(S, S, arma::fill::zeros);
  for (int a = 0; a < S; ++a) {
    for (int b = 0; b < S; ++b) {
      double prob = 1.0;
      for (int k = 0; k < K; ++k) {
        int from_k = (a / intPower(2, k)) % 2;
        int to_k   = (b / intPower(2, k)) % 2;
        prob *= gamma(from_k, to_k);
      }
      jointGamma(a, b) = prob;
    }
  }
  return jointGamma;
}

// [[Rcpp::export]]
arma::mat JointTransitionMatrix_per_strain_cpp2(List gamma_list, int K) {
  int S = intPower(2, K);
  arma::mat jointGamma(S, S, arma::fill::zeros);
  for (int a = 0; a < S; ++a) {
    for (int b = 0; b < S; ++b) {
      double prob = 1.0;
      for (int k = 0; k < K; ++k) {
        int from_k = (a / intPower(2, k)) % 2;
        int to_k   = (b / intPower(2, k)) % 2;
        arma::mat currentGamma = gamma_list[k];
        prob *= currentGamma(from_k, to_k);
      }
      jointGamma(a, b) = prob;
    }
  }
  return jointGamma;
}

void generate_subsets(const std::vector<int>& elements,
                      int idx,
                      std::vector<int>& current,
                      std::vector< std::vector<int> >& out){
  if (idx == (int)elements.size()) {
    out.push_back(current);
    return;
  }

  // Not include element
  generate_subsets(elements, idx + 1, current, out);

  // Include element
  current.push_back(elements[idx]);
  generate_subsets(elements, idx + 1, current, out);
  current.pop_back();
}

// [[Rcpp::export]]
double one_factor_copula_cdf_rcpp(const arma::vec& u,
                                  const arma::vec& Lambdas,
                                  const arma::vec& gh_x,
                                  const arma::vec& gh_w) {

  int K = u.size();
  int M = gh_x.size();

  const double eps = 1e-12;
  arma::vec q(K, arma::fill::zeros);

  // Compute qnorm(u)
  for (int k = 0; k < K; ++k) {
    double uk = std::min(std::max(u[k], eps), 1.0 - eps);
    q[k] = R::qnorm(uk, 0.0, 1.0, 1, 0);
  }

  double result = 0.0;

  for (int j = 0; j < M; ++j) {
    double z = std::sqrt(2.0) * gh_x[j];
    double prod_term = 1.0;

    for (int k = 0; k < K; ++k) {
      double denom = std::sqrt(1.0 - Lambdas[k] * Lambdas[k]);
      double arg = (q[k] - Lambdas[k] * z) / denom;
      prod_term *= R::pnorm(arg, 0.0, 1.0, 1, 0);
    }

    result += gh_w[j] * prod_term;
  }

  return result / std::sqrt(M_PI);
}

// [[Rcpp::export]]
arma::mat JointTransitionMatrix_1FactorGaussiancopula_cpp(const arma::mat& gamma,
                                                           int K,
                                                           const arma::vec& Lambdas,
                                                           const arma::vec& gh_x,
                                                           const arma::vec& gh_w){
  int S = intPower(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);

  arma::mat gamma2 = gamma;
  gamma2(0,0) = gamma(0, 1);
  gamma2(0,1) = gamma(0, 0);

  // parallelize
#pragma omp parallel for schedule(dynamic)
  for (int a = 0; a < S; a++) {
    for (int b = 0; b < S; b++) {

      std::vector<int> Ones;
      std::vector<int> Zeros;
      arma::vec prob(K);

      for (int k = 0; k < K; k++) {
        int from_k = (a >> k) & 1;
        int to_k   = (b >> k) & 1;

        if (from_k == 1)
          Ones.push_back(k);
        else
          Zeros.push_back(k);

        prob(k) = gamma2(from_k, to_k);
      }

      // power set of zeros
      std::vector<std::vector<int>> subsets;
      std::vector<int> cur;
      generate_subsets(Zeros, 0, cur, subsets);

      double total = 0.0;

      for (auto& Tset : subsets) {
        int sign = (Tset.size() % 2 == 0 ? 1 : -1);

        arma::vec u(K, arma::fill::ones);
        for (int idx : Ones) u(idx) = prob(idx);
        for (int idx : Tset) u(idx) = prob(idx);

        total += sign * one_factor_copula_cdf_rcpp(u, Lambdas, gh_x, gh_w);
      }

      GammaMat(a, b) = total;
    }
  }

  for (int i = 0; i < S; i++) {
    double s = arma::accu(GammaMat.row(i));
    GammaMat.row(i) /= s;
  }
  return GammaMat;
}

// [[Rcpp::export]]
arma::mat JointTransitionMatrix_1FactorGaussiancopula_per_strain_cpp(List gamma_list,
                                                                      int K,
                                                                      const arma::vec& Lambdas,
                                                                      const arma::vec& gh_x,
                                                                      const arma::vec& gh_w){
  int S = intPower(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);

  std::vector<arma::mat> gamma_cpp(K);
  for (int k = 0; k < K; k++)
    gamma_cpp[k] = as<arma::mat>(gamma_list[k]);

  // parallelize
#pragma omp parallel for schedule(dynamic)
  for (int a = 0; a < S; a++) {
    for (int b = 0; b < S; b++) {

      std::vector<int> Ones;
      std::vector<int> Zeros;
      arma::vec prob(K);

      for (int k = 0; k < K; k++) {
        arma::mat gamma = gamma_cpp[k];
        arma::mat gamma2 = gamma;
        gamma2(0,0) = gamma(0, 1);
        gamma2(0,1) = gamma(0, 0);
        int from_k = (a >> k) & 1;
        int to_k   = (b >> k) & 1;

        if (from_k == 1)
          Ones.push_back(k);
        else
          Zeros.push_back(k);

        prob(k) = gamma2(from_k, to_k);
      }

      // power set of zeros
      std::vector<std::vector<int>> subsets;
      std::vector<int> cur;
      generate_subsets(Zeros, 0, cur, subsets);

      double total = 0.0;

      for (auto& Tset : subsets) {
        int sign = (Tset.size() % 2 == 0 ? 1 : -1);

        arma::vec u(K, arma::fill::ones);
        for (int idx : Ones) u(idx) = prob(idx);
        for (int idx : Tset) u(idx) = prob(idx);

        total += sign * one_factor_copula_cdf_rcpp(u, Lambdas, gh_x, gh_w);
      }

      GammaMat(a, b) = total;
    }
  }

  for (int i = 0; i < S; i++) {
    double s = arma::accu(GammaMat.row(i));
    GammaMat.row(i) /= s;
  }
  return GammaMat;
}

// [[Rcpp::export]]
double frank_cdf_cpp2(const arma::vec &u, double theta) {
  int d = u.n_elem;

  if (std::abs(theta) < 1e-12) {
    double prod = 1.0;
    for (int i = 0; i < d; i++) {
      double ui = u[i];
      if (ui <= 0.0) return 0.0;
      if (ui >= 1.0) ui = 1.0;
      prod *= ui;
    }
    return prod;
  }

  for (int i = 0; i < d; i++) {
    if (!std::isfinite(u[i])) return NA_REAL;
    if (u[i] < 0.0 || u[i] > 1.0) return NA_REAL;
  }

  double den = std::expm1(-theta);
  if (den == 0.0) return NA_REAL;

  long double prod_term = 1.0L;

  for (int i = 0; i < d; i++) {
    double ui = u[i];
    double term = std::expm1(-theta * ui);
    prod_term *= (long double)term;
  }

  long double ratio = prod_term / std::pow(den, d - 1);

  double result = -(1.0 / theta) * std::log1p((double)ratio);

  if (!std::isfinite(result)) return NA_REAL;
  if (result < 0.0) result = 0.0;
  if (result > 1.0) result = 1.0;

  return result;
}

// [[Rcpp::export]]
arma::mat JointTransitionMatrix_Frankcopula_cpp(const arma::mat& gamma,
                                                        int K,
                                                        double copParam){
  int S = intPower(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);

  arma::mat gamma2 = gamma;
  gamma2(0,0) = gamma(0, 1);
  gamma2(0,1) = gamma(0, 0);

  // parallelize
#pragma omp parallel for schedule(dynamic)
  for (int a = 0; a < S; a++) {
    for (int b = 0; b < S; b++) {

      std::vector<int> Ones;
      std::vector<int> Zeros;
      arma::vec prob(K);

      for (int k = 0; k < K; k++) {
        int from_k = (a >> k) & 1;
        int to_k   = (b >> k) & 1;

        if (from_k == 1)
          Ones.push_back(k);
        else
          Zeros.push_back(k);

        prob(k) = gamma2(from_k, to_k);
      }

      // power set of zeros
      std::vector<std::vector<int>> subsets;
      std::vector<int> cur;
      generate_subsets(Zeros, 0, cur, subsets);

      double total = 0.0;

      for (auto& Tset : subsets) {
        int sign = (Tset.size() % 2 == 0 ? 1 : -1);

        arma::vec u(K, arma::fill::ones);
        for (int idx : Ones) u(idx) = prob(idx);
        for (int idx : Tset) u(idx) = prob(idx);

        total += sign * frank_cdf_cpp2(u, copParam);
      }

      GammaMat(a, b) = total;
    }
  }

  for (int i = 0; i < S; i++) {
    double s = arma::accu(GammaMat.row(i));
    GammaMat.row(i) /= s;
  }
  return GammaMat;
}


// [[Rcpp::export]]
arma::mat JointTransitionMatrix_Frankcopula_perstrain_cpp(List gamma_list,
                                                                  int K,
                                                                  double copParam){
  int S = intPower(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);

  std::vector<arma::mat> gamma_cpp(K);
  for (int k = 0; k < K; k++)
    gamma_cpp[k] = as<arma::mat>(gamma_list[k]);

  // parallelize
#pragma omp parallel for schedule(dynamic)
  for (int a = 0; a < S; a++) {
    for (int b = 0; b < S; b++) {

      std::vector<int> Ones;
      std::vector<int> Zeros;
      arma::vec prob(K);

      for (int k = 0; k < K; k++) {
        arma::mat gamma = gamma_cpp[k];
        arma::mat gamma2 = gamma;
        gamma2(0,0) = gamma(0, 1);
        gamma2(0,1) = gamma(0, 0);
        int from_k = (a >> k) & 1;
        int to_k   = (b >> k) & 1;

        if (from_k == 1)
          Ones.push_back(k);
        else
          Zeros.push_back(k);

        prob(k) = gamma2(from_k, to_k);
      }

      // power set of zeros
      std::vector<std::vector<int>> subsets;
      std::vector<int> cur;
      generate_subsets(Zeros, 0, cur, subsets);

      double total = 0.0;

      for (auto& Tset : subsets) {
        int sign = (Tset.size() % 2 == 0 ? 1 : -1);

        arma::vec u(K, arma::fill::ones);
        for (int idx : Ones) u(idx) = prob(idx);
        for (int idx : Tset) u(idx) = prob(idx);

        total += sign * frank_cdf_cpp2(u, copParam);
      }

      GammaMat(a, b) = total;
    }
  }

  for (int i = 0; i < S; i++) {
    double s = arma::accu(GammaMat.row(i));
    GammaMat.row(i) /= s;
  }
  return GammaMat;
}

// [[Rcpp::export]]
arma::mat makematrix_arma_cpp2(double g12, double g21){
  arma::mat Gmat(2, 2, arma::fill::zeros);
  Gmat(0, 0) = 1 - g12;
  Gmat(0, 1) = g12;
  Gmat(1, 0) = g21;
  Gmat(1, 1) = 1 - g21;
  return Gmat;
}

// [[Rcpp::export]]
Rcpp::List BuildGamma_list_cpp(const arma::vec& Gs) {
  int n_mat = Gs.n_elem / 2;
  Rcpp::List out(n_mat);

  int index = 0;
  for (int i = 0; i < n_mat; i++) {
    double a = Gs[i + index];
    double b = Gs[i + index + 1];

    out[i] = makematrix_arma_cpp2(a, b);
    index += 1;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat build_corr_from_params_cpp(int d, const arma::vec& params) {
  arma::mat R = arma::eye(d, d);
  int idx = 0;

  for (int i = 0; i < d; i++) {
    for (int j = i + 1; j < d; j++) {
      R(i, j) = params(idx);
      R(j, i) = params(idx);
      idx++;
    }
  }
  return R;
}

// [[Rcpp::export]]
double gaussian_copula_cdf_cpp(const arma::vec& u,
                               const arma::mat& corrMat){
  int d = u.n_elem;

  arma::vec u_clip = u;
  double eps = 1e-12;
  for (int i = 0; i < d; i++) {
    if (u_clip(i) < eps)     u_clip(i) = eps;
    if (u_clip(i) > 1 - eps) u_clip(i) = 1 - eps;
  }

  // Probit transform
  NumericVector upper(d);
  for (int i = 0; i < d; i++) {
    upper[i] = R::qnorm(u_clip(i), 0.0, 1.0, true, false);
  }

  // Convert lower to numeric vector
  NumericVector lower(d, -std::numeric_limits<double>::infinity());

  NumericMatrix sigma = wrap(corrMat);
  NumericVector mean(d, 0.0);

  Environment mvtnorm = Environment::namespace_env("mvtnorm");
  Function pmvnorm = mvtnorm["pmvnorm"];

  SEXP result = pmvnorm(
    _["lower"] = lower,
    _["upper"] = upper,
    _["mean"]  = mean,
    _["sigma"] = sigma
  );

  return as<double>(result);
}

// [[Rcpp::export]]
arma::mat JointTransitionMatrix_Gaussiancopula_cpp(const arma::mat& gamma,
                                           int K,
                                           const arma::vec& copParams){
  arma::mat corrMat = build_corr_from_params_cpp(K, copParams);
  int S = intPower(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);

  arma::mat gamma2 = gamma;
  gamma2(0,0) = gamma(0, 1);
  gamma2(0,1) = gamma(0, 0);

  for (int a = 0; a < S; a++) {
    for (int b = 0; b < S; b++) {

      std::vector<int> Ones;
      std::vector<int> Zeros;
      arma::vec prob(K);

      for (int k = 0; k < K; k++) {
        int from_k = (a >> k) & 1;
        int to_k   = (b >> k) & 1;

        if (from_k == 1)
          Ones.push_back(k);
        else
          Zeros.push_back(k);

        prob(k) = gamma2(from_k, to_k);
      }

      // power set of zeros
      std::vector< std::vector<int> > subsets;
      std::vector<int> cur;
      generate_subsets(Zeros, 0, cur, subsets);

      double total = 0.0;

      for (auto& Tset : subsets) {
        int sign = (Tset.size() % 2 == 0 ? 1 : -1);

        arma::vec u(K, arma::fill::ones);
        for (int idx : Ones) u(idx) = prob(idx);
        for (int idx : Tset) u(idx) = prob(idx);

        total += sign *  gaussian_copula_cdf_cpp(u, corrMat);
      }
      GammaMat(a,b) = total;
    }
  }

  for (int i = 0; i < S; i++) {
    double s = arma::accu(GammaMat.row(i));
    GammaMat.row(i) /= s;
  }
  return GammaMat;
}

// [[Rcpp::export]]
arma::mat JointTransitionMatrix_Gaussiancopula_perstrain_cpp(List gamma_list,
                                                     int K,
                                                     const arma::vec& copParams){
  arma::mat corrMat = build_corr_from_params_cpp(K, copParams);
  int S = intPower(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);



  for (int a = 0; a < S; a++) {
    for (int b = 0; b < S; b++) {

      std::vector<int> Ones;
      std::vector<int> Zeros;
      arma::vec prob(K);

      for (int k = 0; k < K; k++) {
        arma::mat gamma = gamma_list[k];
        arma::mat gamma2 = gamma;
        gamma2(0,0) = gamma(0, 1);
        gamma2(0,1) = gamma(0, 0);
        int from_k = (a >> k) & 1;
        int to_k   = (b >> k) & 1;

        if (from_k == 1)
          Ones.push_back(k);
        else
          Zeros.push_back(k);

        prob(k) = gamma2(from_k, to_k);
      }

      // power set of zeros
      std::vector< std::vector<int> > subsets;
      std::vector<int> cur;
      generate_subsets(Zeros, 0, cur, subsets);

      double total = 0.0;

      for (auto& Tset : subsets) {
        int sign = (Tset.size() % 2 == 0 ? 1 : -1);

        arma::vec u(K, arma::fill::ones);
        for (int idx : Ones) u(idx) = prob(idx);
        for (int idx : Tset) u(idx) = prob(idx);

        total += sign *  gaussian_copula_cdf_cpp(u, corrMat);
      }
      GammaMat(a,b) = total;
    }
  }

  for (int i = 0; i < S; i++) {
    double s = arma::accu(GammaMat.row(i));
    GammaMat.row(i) /= s;
  }
  return GammaMat;
}

// [[Rcpp::export]]
arma::mat replace_naMat_with_zero(arma::mat X){
  X.elem(arma::find_nonfinite(X)).zeros();
  return X;
}

// [[Rcpp::export]]
arma::vec replace_naVec_with_zero(arma::vec X){
  X.elem(arma::find_nonfinite(X)).zeros();
  return X;
}

// [[Rcpp::export]]
double add_untypedPoissonLoglikelihood(arma::cube y, arma::cube allPoisMean, arma::mat y_total){
  int ndept = y_total.n_rows;
  int time = y_total.n_cols;
  double singlePoissonLoglikelihood = 0;
  arma::mat y_strain1 = y.slice(0);
  for(int i = 0; i < ndept; ++i){
    for (int t = 0; t < time; ++t){
      if(arma::is_finite(y_total(i,t)) && !arma::is_finite(y_strain1(i,t))){
        double sumRisks = arma::accu(allPoisMean.tube(i, t));
        sumRisks = std::max(sumRisks, 1e-12);
        singlePoissonLoglikelihood += y_total(i, t) * log(sumRisks) - sumRisks - lgamma(y_total(i, t) + 1);
      }
    }
  }
  return singlePoissonLoglikelihood;
}

// [[Rcpp::export]]
arma::mat add_untyped_delta(arma::cube y, arma::cube allPoisMean, arma::mat y_total, arma::mat delta){
  int ndept = y_total.n_rows;
  int time = y_total.n_cols;
  arma::mat newDelta = delta;
  arma::mat y_strain1 = y.slice(0);
  for(int i = 0; i < ndept; ++i){
    for (int t = 0; t < time; ++t){
      if(arma::is_finite(y_total(i,t)) && !arma::is_finite(y_strain1(i,t))){
        double sumRisks = arma::accu(allPoisMean.tube(i, t));
        sumRisks = std::max(sumRisks, 1e-12);
        newDelta(i, t) = y_total(i, t) - sumRisks;
      }
    }
  }
  return newDelta;
}

// [[Rcpp::export]]
double add_untyped_logemission(arma::vec y_vec, arma::vec lambda_vec, double y_totalscalar){
  double untypedlogEmission = 0;
  if(arma::is_finite(y_totalscalar) && !arma::is_finite(y_vec[0])){
    double sumRisks = arma::accu(lambda_vec);
    sumRisks = std::max(sumRisks, 1e-12);
    untypedlogEmission = y_totalscalar * log(sumRisks) - sumRisks - lgamma(y_totalscalar + 1);
  }
  return untypedlogEmission;
}

// [[Rcpp::export]]
List SMOOTHINGgradmultstrainLoglikelihood_cpp(arma::cube y, arma::mat e_it, int nstrain, arma::vec r, arma::vec s,
                                              arma::vec u, arma::mat jointTPM, arma::vec B, arma::mat Bits, arma::vec a_k,
                                              int Model, arma::mat Q_r, arma::mat Q_s, arma::mat Q_u, int gradients,
                                              arma::mat Qstz_r, arma::mat Qstz_s, arma::mat Qstz_u, arma::mat y_total){

  int ndept = e_it.n_rows;
  int time = e_it.n_cols;
  int nstate = intPower(2, nstrain);

  if(Model == 0){
    arma::uvec month_indexes(time);
    for (int t = 0; t < time; t++) {
      month_indexes(t) = (t % 12);
    }
    arma::mat r_mat = arma::repmat(r.t(), ndept, 1);

    arma::vec s_sub = s.elem(month_indexes);
    arma::mat s_mat = arma::repmat(s_sub.t(), ndept, 1);

    arma::mat u_mat = arma::repmat(u, 1, time);

    arma::mat log_risk = r_mat + s_mat + u_mat;

    arma::mat poisMean(ndept, time, arma::fill::zeros);
    arma::cube allPoisMean(ndept, time, nstrain, arma::fill::zeros);
    arma::mat delta(ndept, time, arma::fill::zeros);

    for (int k = 0; k < nstrain; ++k) {
      arma::mat lambda = e_it % arma::exp(log_risk + a_k[k]);
      arma::mat delta2 = y.slice(k) - lambda;
      delta2 = replace_naMat_with_zero(delta2);
      delta += delta2;
      poisMean += lambda;
      allPoisMean.slice(k) = lambda;
    }

    // compute log-likelihood
    double loglike = 0.0;
    for (int k = 0; k < nstrain; ++k){
      arma::mat Y = y.slice(k);
      arma::mat Lambda = allPoisMean.slice(k);
      arma::mat safeLambda = Lambda;
      safeLambda.transform( [](double val) { return (val <= 0) ? 1e-12 : val; } );
      arma::mat tempPoisDensity = Y % arma::log(safeLambda) - Lambda - lgamma(Y + 1);
      tempPoisDensity = replace_naMat_with_zero(tempPoisDensity);
      loglike += arma::accu(tempPoisDensity);
    }

    loglike += add_untypedPoissonLoglikelihood(y, allPoisMean, y_total);
    delta = add_untyped_delta(y, allPoisMean, y_total, delta);

    // Temporal trend r gradients
    arma::vec grad_r = arma::sum(delta, 0).t() - Q_r * r;
    grad_r = Qstz_r.t() * grad_r;
    arma::mat diag_pois_colsum = arma::diagmat(arma::sum(poisMean, 0));
    arma::mat cov_r = arma::inv_sympd(diag_pois_colsum + Q_r + arma::eye(time, time) * 1e-8);
    cov_r = Qstz_r.t() * cov_r * Qstz_r;

    // Seasonal s gradients
    arma::vec grad_s(12, arma::fill::zeros);
    arma::vec fishervec_s(12, arma::fill::zeros);

    for (int month_index = 0; month_index < 12; ++month_index) {
      for (int t = 0; t < time; ++t) {
        if ((t % 12) == month_index) {
          grad_s(month_index)     += arma::accu(delta.col(t));
          fishervec_s(month_index)+= arma::accu(poisMean.col(t));
        }
      }
    }
    grad_s -= Q_s * s;
    grad_s = Qstz_s.t() * grad_s;
    arma::mat cov_s = arma::inv_sympd(arma::diagmat(fishervec_s) + Q_s);
    cov_s = Qstz_s.t() * cov_s * Qstz_s;


    // Spatial u gradients
    arma::vec grad_u = arma::sum(delta, 1) - Q_u * u;
    grad_u = Qstz_u.t() * grad_u;
    arma::mat diag_pois_rowsum = arma::diagmat(arma::sum(poisMean, 1));
    arma::mat cov_u = arma::inv_sympd(diag_pois_rowsum + Q_u + arma::eye(ndept, ndept) * 1e-8);
    cov_u = Qstz_u.t() * cov_u * Qstz_u;

    double poisMean4GibbsUpdate = arma::accu(e_it % arma::exp(log_risk));

    return List::create(
      Named("loglike") = loglike,
      Named("grad_r") = grad_r,
      Named("grad_s") = grad_s,
      Named("grad_u") = grad_u,
      Named("cov_r") = cov_r,
      Named("cov_s") = cov_s,
      Named("cov_u") = cov_u,
      Named("poisMean4GibbsUpdate") = poisMean4GibbsUpdate
    );
  }else{

    double loglike_total = 0.0;

    arma::mat safeTPM = jointTPM;
    safeTPM.transform([](double val){ return (val <= 0) ? 1e-12 : val; });
    arma::mat logjointTPM = arma::log(safeTPM);
    arma::mat logjointTPM_t = logjointTPM.t();

    arma::vec init_density = stationarydistArma_cpp(jointTPM);
    arma::vec safeinitdensity = init_density;
    safeinitdensity.transform( [](double val) { return (val <= 0) ? 1e-12 : val; });
    arma::vec loginit_density = arma::log(safeinitdensity);

    arma::cube E_lambda_itk(ndept, time, nstrain, arma::fill::zeros);
    arma::cube E_lambda_itk2(ndept, time, nstrain, arma::fill::zeros);

    arma::vec grad_r(time, arma::fill::zeros);
    arma::vec grad_s(12, arma::fill::zeros);
    arma::vec grad_u(ndept, arma::fill::zeros);
    arma::mat cov_r(time, time, arma::fill::zeros);
    arma::mat cov_s(12, 12, arma::fill::zeros);
    arma::mat cov_u(ndept, ndept, arma::fill::zeros);
    arma::vec poisMean4GibbsUpdate(nstrain, arma::fill::zeros);

    if(gradients == 0){

      for(int i = 0; i < ndept; ++i){
        arma::mat logEmissions(time, nstate, arma::fill::zeros);
        arma::cube lambda_array(time, nstate, nstrain, arma::fill::zeros);
        for(int t = 0; t < time; ++t){
          int month_index = t % 12;
          for(int n = 0; n < nstate; ++n){
            for(int k = 0; k < nstrain; ++k){
              lambda_array(t, n, k) = e_it(i, t) * std::exp(a_k[k] + r[t] + s[month_index] + u[i] + Bits(n, k) * B[k]);
            }
            arma::vec y_vec = y.tube(i, t);
            arma::vec lambda_vec = lambda_array.tube(t, n);
            arma::vec safelambda_vec = lambda_vec;
            safelambda_vec.transform( [](double val) { return (val <= 0) ? 1e-12 : val; });
            arma::vec tempPoisDensity = y_vec % arma::log(safelambda_vec) - lambda_vec - lgamma(y_vec + 1);
            tempPoisDensity = replace_naVec_with_zero(tempPoisDensity);
            logEmissions(t, n) = arma::accu(tempPoisDensity);
            logEmissions(t, n) += add_untyped_logemission(y_vec, lambda_vec, y_total(i, t));
          }
        }
        //forward pass
        arma::mat logalpha(time, nstate, arma::fill::zeros);

        logalpha.row(0) = loginit_density.t() + logEmissions.row(0);
        for(int t = 1; t < time; ++t){
          logalpha.row(t) = (logVecMatMult2(logalpha.row(t-1).t(), logjointTPM) + logEmissions.row(t).t()).t();
        }

        double loglike_i = logSumExp_cpp2(logalpha.row(time-1).t());
        loglike_total += loglike_i;
      }
    }else{

      for(int i = 0; i < ndept; ++i){
        arma::mat logEmissions(time, nstate, arma::fill::zeros);
        arma::cube lambda_array(time, nstate, nstrain, arma::fill::zeros);
        arma::cube lambda_array2(time, nstate, nstrain, arma::fill::zeros);
        for(int t = 0; t < time; ++t){
          int month_index = t % 12;
          for(int n = 0; n < nstate; ++n){
            for(int k = 0; k < nstrain; ++k){
              lambda_array(t, n, k) = e_it(i, t) * std::exp(a_k[k] + r[t] + s[month_index] + u[i] + Bits(n, k) * B[k]);
              lambda_array2(t, n, k) = e_it(i, t) * std::exp(r[t] + s[month_index] + u[i] + Bits(n, k) * B[k]);
            }
            arma::vec y_vec = y.tube(i, t);
            arma::vec lambda_vec = lambda_array.tube(t, n);
            arma::vec safelambda_vec = lambda_vec;
            safelambda_vec.transform( [](double val) { return (val <= 0) ? 1e-12 : val; });
            arma::vec tempPoisDensity = y_vec % arma::log(safelambda_vec) - lambda_vec - lgamma(y_vec + 1);
            tempPoisDensity = replace_naVec_with_zero(tempPoisDensity);
            logEmissions(t, n) = arma::accu(tempPoisDensity);
            logEmissions(t, n) += add_untyped_logemission(y_vec, lambda_vec, y_total(i, t));
          }
        }
        //forward pass
        arma::mat logalpha(time, nstate, arma::fill::zeros);

        logalpha.row(0) = loginit_density.t() + logEmissions.row(0);
        for(int t = 1; t < time; ++t){
          logalpha.row(t) = (logVecMatMult2(logalpha.row(t-1).t(), logjointTPM) + logEmissions.row(t).t()).t();
        }

        double loglike_i = logSumExp_cpp2(logalpha.row(time-1).t());
        loglike_total += loglike_i;

        //backward pass
        arma::mat logbeta(time, nstate, arma::fill::zeros);

        for(int t = time - 2; t >= 0; --t){
          arma::vec vec = (logEmissions.row(t + 1) + logbeta.row(t + 1)).t();
          logbeta.row(t) = logVecMatMult2(vec, logjointTPM_t).t();
        }
        //Marginal posterior probabilities and smoothing
        arma::mat logP_s = (logalpha + logbeta) - loglike_i;
        arma::mat P_s = arma::exp(logP_s);

        for(int t = 0; t < time; ++t){
          for(int k = 0; k < nstrain; ++k){
            arma::vec newlambtube(nstate, arma::fill::zeros);
            arma::vec newlambtube2(nstate, arma::fill::zeros);
            for(int n = 0; n < nstate; ++n){
              newlambtube[n] = lambda_array(t, n, k);
              newlambtube2[n] = lambda_array2(t, n, k);
            }
            arma::vec probvec = P_s.row(t).t();
            E_lambda_itk(i, t, k) = arma::dot(probvec, newlambtube);
            E_lambda_itk2(i, t, k) = arma::dot(probvec, newlambtube2);
          }
        }
      }

      arma::mat poisMean(ndept, time, arma::fill::zeros);
      arma::mat delta(ndept, time, arma::fill::zeros);

      for (int k = 0; k < nstrain; ++k){
        arma::mat currentY = y.slice(k);
        arma::mat currentE_lambda = E_lambda_itk.slice(k);
        arma::mat currentE_lambda2 = E_lambda_itk2.slice(k);
        arma::mat delta2   = currentY - currentE_lambda;
        delta2 = replace_naMat_with_zero(delta2);
        delta += delta2;
        poisMean += currentE_lambda;
        poisMean4GibbsUpdate[k] = arma::accu(currentE_lambda2);
      }

      delta = add_untyped_delta(y, E_lambda_itk, y_total, delta);

      // Temporal trend r gradients
      grad_r = arma::sum(delta, 0).t() - Q_r * r;
      grad_r = Qstz_r.t() * grad_r;
      arma::mat diag_pois_colsum = arma::diagmat(arma::sum(poisMean, 0));
      cov_r = arma::inv_sympd(diag_pois_colsum + Q_r + arma::eye(time, time) * 1e-8);
      cov_r = Qstz_r.t() * cov_r * Qstz_r;

      // Seasonal s gradients
      arma::vec fishervec_s(12, arma::fill::zeros);

      for(int month_index = 0; month_index < 12; ++month_index){
        for(int t = 0; t < time; ++t){
          if((t % 12) == month_index){
            grad_s(month_index)     += arma::accu(delta.col(t));
            fishervec_s(month_index)+= arma::accu(poisMean.col(t));
          }
        }
      }
      grad_s -= Q_s * s;
      grad_s = Qstz_s.t() * grad_s;
      cov_s = arma::inv_sympd(arma::diagmat(fishervec_s) + Q_s);
      cov_s = Qstz_s.t() * cov_s * Qstz_s;


      // Spatial u gradients
      grad_u = arma::sum(delta, 1) - Q_u * u;
      grad_u = Qstz_u.t() * grad_u;
      arma::mat diag_pois_rowsum = arma::diagmat(arma::sum(poisMean, 1));
      cov_u = arma::inv_sympd(diag_pois_rowsum + Q_u + arma::eye(ndept, ndept) * 1e-8);
      cov_u = Qstz_u.t() * cov_u * Qstz_u;
    }

    return List::create(
      Named("loglike") = loglike_total,
      Named("grad_r") = grad_r,
      Named("grad_s") = grad_s,
      Named("grad_u") = grad_u,
      Named("cov_r") = cov_r,
      Named("cov_s") = cov_s,
      Named("cov_u") = cov_u,
      Named("poisMean4GibbsUpdate") = poisMean4GibbsUpdate
    );
  }
}

// [[Rcpp::export]]
List FFBSgradmultstrainLoglikelihood_cpp(arma::cube y, arma::mat e_it, int nstrain, arma::vec r, arma::vec s,
                                         arma::vec u, arma::mat jointTPM, arma::vec B, arma::mat Bits, arma::vec a_k,
                                         int Model, arma::mat Q_r, arma::mat Q_s, arma::mat Q_u, int gradients,
                                         arma::mat Qstz_r, arma::mat Qstz_s, arma::mat Qstz_u, arma::mat y_total){

  int ndept = e_it.n_rows;
  int time = e_it.n_cols;
  int nstate = intPower(2, nstrain);

  if(Model == 0){
    arma::uvec month_indexes(time);
    for (int t = 0; t < time; t++) {
      month_indexes(t) = (t % 12);
    }
    arma::mat r_mat = arma::repmat(r.t(), ndept, 1);

    arma::vec s_sub = s.elem(month_indexes);
    arma::mat s_mat = arma::repmat(s_sub.t(), ndept, 1);

    arma::mat u_mat = arma::repmat(u, 1, time);

    arma::mat log_risk = r_mat + s_mat + u_mat;

    arma::mat poisMean(ndept, time, arma::fill::zeros);
    arma::cube allPoisMean(ndept, time, nstrain, arma::fill::zeros);
    arma::mat delta(ndept, time, arma::fill::zeros);

    for (int k = 0; k < nstrain; ++k) {
      arma::mat lambda = e_it % arma::exp(log_risk + a_k[k]);
      arma::mat delta2 = y.slice(k) - lambda;
      delta2 = replace_naMat_with_zero(delta2);
      delta += delta2;

      poisMean += lambda;
      allPoisMean.slice(k) = lambda;
    }

    // compute log-likelihood
    double loglike = 0.0;
    for (int k = 0; k < nstrain; ++k){
      arma::mat Y = y.slice(k);
      arma::mat Lambda = allPoisMean.slice(k);
      arma::mat safeLambda = Lambda;
      safeLambda.transform( [](double val) { return (val <= 0) ? 1e-12 : val; } );
      arma::mat tempPoisDensity = Y % arma::log(safeLambda) - Lambda - lgamma(Y + 1);
      tempPoisDensity = replace_naMat_with_zero(tempPoisDensity);
      loglike += arma::accu(tempPoisDensity);
    }

    loglike += add_untypedPoissonLoglikelihood(y, allPoisMean, y_total);
    delta = add_untyped_delta(y, allPoisMean, y_total, delta);

    // Temporal trend r gradients
    arma::vec grad_r = arma::sum(delta, 0).t() - Q_r * r;
    grad_r = Qstz_r.t() * grad_r;
    arma::mat diag_pois_colsum = arma::diagmat(arma::sum(poisMean, 0));
    arma::mat cov_r = arma::inv_sympd(diag_pois_colsum + Q_r + arma::eye(time, time) * 1e-8);
    cov_r = Qstz_r.t() * cov_r * Qstz_r;

    // Seasonal s gradients
    arma::vec grad_s(12, arma::fill::zeros);
    arma::vec fishervec_s(12, arma::fill::zeros);

    for (int month_index = 0; month_index < 12; ++month_index) {
      for (int t = 0; t < time; ++t) {
        if ((t % 12) == month_index) {
          grad_s(month_index)     += arma::accu(delta.col(t));
          fishervec_s(month_index)+= arma::accu(poisMean.col(t));
        }
      }
    }
    grad_s -= Q_s * s;
    grad_s = Qstz_s.t() * grad_s;
    arma::mat cov_s = arma::inv_sympd(arma::diagmat(fishervec_s) + Q_s);
    cov_s = Qstz_s.t() * cov_s * Qstz_s;


    // Spatial u gradients
    arma::vec grad_u = arma::sum(delta, 1) - Q_u * u;
    grad_u = Qstz_u.t() * grad_u;
    arma::mat diag_pois_rowsum = arma::diagmat(arma::sum(poisMean, 1));
    arma::mat cov_u = arma::inv_sympd(diag_pois_rowsum + Q_u + arma::eye(ndept, ndept) * 1e-8);
    cov_u = Qstz_u.t() * cov_u * Qstz_u;

    double poisMean4GibbsUpdate = arma::accu(e_it % arma::exp(log_risk));

    return List::create(
      Named("loglike") = loglike,
      Named("grad_r") = grad_r,
      Named("grad_s") = grad_s,
      Named("grad_u") = grad_u,
      Named("cov_r") = cov_r,
      Named("cov_s") = cov_s,
      Named("cov_u") = cov_u,
      Named("poisMean4GibbsUpdate") = poisMean4GibbsUpdate
    );
  }else{

    double loglike_total = 0.0;

    arma::mat safeTPM = jointTPM;
    safeTPM.transform([](double val){ return (val <= 0) ? 1e-12 : val; });
    arma::mat logjointTPM = arma::log(safeTPM);
    arma::mat logjointTPM_t = logjointTPM.t();

    arma::vec init_density = stationarydistArma_cpp(jointTPM);
    arma::vec safeinitdensity = init_density;
    safeinitdensity.transform( [](double val) { return (val <= 0) ? 1e-12 : val; });
    arma::vec loginit_density = arma::log(safeinitdensity);

    arma::cube Actual_lambda_itk(ndept, time, nstrain, arma::fill::zeros);
    arma::cube Actual_lambda_itk2(ndept, time, nstrain, arma::fill::zeros);

    arma::vec grad_r(time, arma::fill::zeros);
    arma::vec grad_s(12, arma::fill::zeros);
    arma::vec grad_u(ndept, arma::fill::zeros);
    arma::mat cov_r(time, time, arma::fill::zeros);
    arma::mat cov_s(12, 12, arma::fill::zeros);
    arma::mat cov_u(ndept, ndept, arma::fill::zeros);
    arma::vec poisMean4GibbsUpdate(nstrain, arma::fill::zeros);

    if(gradients == 0){

      for(int i = 0; i < ndept; ++i){
        arma::mat logEmissions(time, nstate, arma::fill::zeros);
        arma::cube lambda_array(time, nstate, nstrain, arma::fill::zeros);
        for(int t = 0; t < time; ++t){
          int month_index = t % 12;
          for(int n = 0; n < nstate; ++n){
            for(int k = 0; k < nstrain; ++k){
              lambda_array(t, n, k) = e_it(i, t) * std::exp(a_k[k] + r[t] + s[month_index] + u[i] + Bits(n, k) * B[k]);
            }
            arma::vec y_vec = y.tube(i, t);
            arma::vec lambda_vec = lambda_array.tube(t, n);
            arma::vec safelambda_vec = lambda_vec;
            safelambda_vec.transform( [](double val) { return (val <= 0) ? 1e-12 : val; });
            arma::vec tempPoisDensity = y_vec % arma::log(safelambda_vec) - lambda_vec - lgamma(y_vec + 1);
            tempPoisDensity = replace_naVec_with_zero(tempPoisDensity);
            logEmissions(t, n) = arma::accu(tempPoisDensity);
            logEmissions(t, n) += add_untyped_logemission(y_vec, lambda_vec, y_total(i, t));
          }
        }

        //forward filtering
        arma::mat logalpha(time, nstate, arma::fill::zeros);

        logalpha.row(0) = loginit_density.t() + logEmissions.row(0);
        for(int t = 1; t < time; ++t){
          logalpha.row(t) = (logVecMatMult2(logalpha.row(t-1).t(), logjointTPM) + logEmissions.row(t).t()).t();
        }

        double loglike_i = logSumExp_cpp2(logalpha.row(time-1).t());
        loglike_total += loglike_i;
      }
    }else{

      for(int i = 0; i < ndept; ++i){
        arma::mat logEmissions(time, nstate, arma::fill::zeros);
        arma::cube lambda_array(time, nstate, nstrain, arma::fill::zeros);
        for(int t = 0; t < time; ++t){
          int month_index = t % 12;
          for(int n = 0; n < nstate; ++n){
            for(int k = 0; k < nstrain; ++k){
              lambda_array(t, n, k) = e_it(i, t) * std::exp(a_k[k] + r[t] + s[month_index] + u[i] + Bits(n, k) * B[k]);
            }
            arma::vec y_vec = y.tube(i, t);
            arma::vec lambda_vec = lambda_array.tube(t, n);
            arma::vec safelambda_vec = lambda_vec;
            safelambda_vec.transform( [](double val) { return (val <= 0) ? 1e-12 : val; });
            arma::vec tempPoisDensity = y_vec % arma::log(safelambda_vec) - lambda_vec - lgamma(y_vec + 1);
            tempPoisDensity = replace_naVec_with_zero(tempPoisDensity);
            logEmissions(t, n) = arma::accu(tempPoisDensity);
            logEmissions(t, n) += add_untyped_logemission(y_vec, lambda_vec, y_total(i, t));
          }
        }
        //forward filtering
        arma::mat logalpha(time, nstate, arma::fill::zeros);

        logalpha.row(0) = loginit_density.t() + logEmissions.row(0);
        for(int t = 1; t < time; ++t){
          logalpha.row(t) = (logVecMatMult2(logalpha.row(t-1).t(), logjointTPM) + logEmissions.row(t).t()).t();
        }

        double loglike_i = logSumExp_cpp2(logalpha.row(time-1).t());
        loglike_total += loglike_i;

        //backward sampling
        arma::uvec states(time, arma::fill::zeros);
        arma::vec statevec = arma::regspace(0, nstate-1);
        arma::vec logP_T = logalpha.row(time-1).t() - loglike_i;
        arma::vec P_T = arma::exp(logP_T - logSumExp_cpp2(logP_T));
        states[time-1] = Rcpp::RcppArmadillo::sample(statevec, 1, false, P_T)[0];

        arma::mat logbeta(time, nstate, arma::fill::zeros);

        for(int t = time - 2; t >= 0; --t){
          arma::vec logP_t = logalpha.row(t).t() + logjointTPM.col(states[t + 1]);
          arma::vec P_t = arma::exp(logP_t - logSumExp_cpp2(logP_t));
          states[t] = Rcpp::RcppArmadillo::sample(statevec, 1, false, P_t)[0];
        }

        //Actual Poisson mean
        for(int t = 0; t < time; ++t){
          int nmonth_index = t % 12;
          for(int k = 0; k < nstrain; ++k){
            Actual_lambda_itk(i,t,k) = e_it(i, t) * std::exp(a_k[k] + r[t] + s[nmonth_index] + u[i] + Bits(states[t], k) * B[k]);
            Actual_lambda_itk2(i,t,k) = e_it(i, t) * std::exp(r[t] + s[nmonth_index] + u[i] + Bits(states[t], k) * B[k]);
          }
        }
      }

      arma::mat poisMean(ndept, time, arma::fill::zeros);
      arma::mat delta(ndept, time, arma::fill::zeros);

      for (int k = 0; k < nstrain; ++k){
        arma::mat currentY = y.slice(k);
        arma::mat currentActual_lambda = Actual_lambda_itk.slice(k);
        arma::mat currentActual_lambda2 = Actual_lambda_itk2.slice(k);
        arma::mat delta2   = currentY - currentActual_lambda;
        delta2 = replace_naMat_with_zero(delta2);
        delta += delta2;
        poisMean += currentActual_lambda;
        poisMean4GibbsUpdate[k] = arma::accu(currentActual_lambda2);
      }

      delta = add_untyped_delta(y, Actual_lambda_itk, y_total, delta);

      // Temporal trend r gradients
      grad_r = arma::sum(delta, 0).t() - Q_r * r;
      grad_r = Qstz_r.t() * grad_r;
      arma::mat diag_pois_colsum = arma::diagmat(arma::sum(poisMean, 0));
      cov_r = arma::inv_sympd(diag_pois_colsum + Q_r + arma::eye(time, time) * 1e-8);
      cov_r = Qstz_r.t() * cov_r * Qstz_r;

      // Seasonal s gradients
      arma::vec fishervec_s(12, arma::fill::zeros);

      for(int month_index = 0; month_index < 12; ++month_index){
        for(int t = 0; t < time; ++t){
          if((t % 12) == month_index){
            grad_s(month_index)     += arma::accu(delta.col(t));
            fishervec_s(month_index)+= arma::accu(poisMean.col(t));
          }
        }
      }
      grad_s -= Q_s * s;
      grad_s = Qstz_s.t() * grad_s;
      cov_s = arma::inv_sympd(arma::diagmat(fishervec_s) + Q_s);
      cov_s = Qstz_s.t() * cov_s * Qstz_s;

      // Spatial u gradients
      grad_u = arma::sum(delta, 1) - Q_u * u;
      grad_u = Qstz_u.t() * grad_u;
      arma::mat diag_pois_rowsum = arma::diagmat(arma::sum(poisMean, 1));
      cov_u = arma::inv_sympd(diag_pois_rowsum + Q_u + arma::eye(ndept, ndept) * 1e-8);
      cov_u = Qstz_u.t() * cov_u * Qstz_u;
    }

    return List::create(
      Named("loglike") = loglike_total,
      Named("grad_r") = grad_r,
      Named("grad_s") = grad_s,
      Named("grad_u") = grad_u,
      Named("cov_r") = cov_r,
      Named("cov_s") = cov_s,
      Named("cov_u") = cov_u,
      Named("poisMean4GibbsUpdate") = poisMean4GibbsUpdate
    );
  }
}

// [[Rcpp::export]]
arma::cube PostOutbreakProbs_cpp(arma::cube y, arma::mat e_it, int nstrain, arma::vec r, arma::vec s,
                                              arma::vec u, arma::mat jointTPM, arma::vec B, arma::mat Bits,
                                              arma::vec a_k, arma::mat y_total){

  int ndept = e_it.n_rows;
  int time = e_it.n_cols;
  int nstate = intPower(2, nstrain);

  arma::cube PosteriorProbabilities_array(ndept, time, nstate, arma::fill::zeros);

    arma::mat safeTPM = jointTPM;
    safeTPM.transform([](double val){ return (val <= 0) ? 1e-12 : val; });
    arma::mat logjointTPM = arma::log(safeTPM);
    arma::mat logjointTPM_t = logjointTPM.t();

    arma::vec init_density = stationarydistArma_cpp(jointTPM);
    arma::vec safeinitdensity = init_density;
    safeinitdensity.transform( [](double val) { return (val <= 0) ? 1e-12 : val; });
    arma::vec loginit_density = arma::log(safeinitdensity);

      for(int i = 0; i < ndept; ++i){
        arma::mat logEmissions(time, nstate, arma::fill::zeros);
        arma::cube lambda_array(time, nstate, nstrain, arma::fill::zeros);
        for(int t = 0; t < time; ++t){
          int month_index = t % 12;
          for(int n = 0; n < nstate; ++n){
            for(int k = 0; k < nstrain; ++k){
              lambda_array(t, n, k) = e_it(i, t) * std::exp(a_k[k] + r[t] + s[month_index] + u[i] + Bits(n, k) * B[k]);
            }
            arma::vec y_vec = y.tube(i, t);
            arma::vec lambda_vec = lambda_array.tube(t, n);
            arma::vec safelambda_vec = lambda_vec;
            safelambda_vec.transform( [](double val) { return (val <= 0) ? 1e-12 : val; });
            arma::vec tempPoisDensity = y_vec % arma::log(safelambda_vec) - lambda_vec - lgamma(y_vec + 1);
            tempPoisDensity = replace_naVec_with_zero(tempPoisDensity);
            logEmissions(t, n) = arma::accu(tempPoisDensity);
            logEmissions(t, n) += add_untyped_logemission(y_vec, lambda_vec, y_total(i, t));
          }
        }
        //forward pass
        arma::mat logalpha(time, nstate, arma::fill::zeros);

        logalpha.row(0) = loginit_density.t() + logEmissions.row(0);
        for(int t = 1; t < time; ++t){
          logalpha.row(t) = (logVecMatMult2(logalpha.row(t-1).t(), logjointTPM) + logEmissions.row(t).t()).t();
        }

        double loglike_i = logSumExp_cpp2(logalpha.row(time-1).t());

        //backward pass
        arma::mat logbeta(time, nstate, arma::fill::zeros);

        for(int t = time - 2; t >= 0; --t){
          arma::vec vec = (logEmissions.row(t + 1) + logbeta.row(t + 1)).t();
          logbeta.row(t) = logVecMatMult2(vec, logjointTPM_t).t();
        }
        //Marginal posterior probabilities
        arma::mat logP_s = (logalpha + logbeta) - loglike_i;
        arma::mat P_s = arma::exp(logP_s);
        for (int t = 0; t < time; ++t) {
          for (int n = 0; n < nstate; ++n) {
            PosteriorProbabilities_array(i, t, n) = P_s(t, n);
          }
        }
      }

    return PosteriorProbabilities_array;
}
