// ConoDictor superfamily classifier: hierarchical conditional-logit model.
//
// Classes: K conopeptide (super)families + NONE (reference, eta = 0).
// Features are attributes of each candidate family: X[k][i, p] is sequence
// i's (standardised, out-of-fold) score against family k's profile p
// (p = region x method, e.g. SIG:hmm, MAT:pssm).
//
//   eta[i, k] = alpha[k] + sum_p X[k][i, p] * B[k, p] + sum_p M[k][i, p] * gamma[p]
//   B[k, p]   = beta[p] + tau[p] * delta_raw[k, p]          (partial pooling)
//
// beta: how much a score on a family's own profile counts as evidence for
// that family, shared across families; tau: how much families differ from
// the shared weights. Rare families shrink towards beta.
// M marks out-of-fold cells with no profile (holdout too small to build);
// gamma absorbs them and is not used at prediction time.
// w: case weights (1 / size of the entry group), so redundant groups do not
// dominate.
//
// Parallelism: the log-likelihood is a sum over sequences, computed with
// reduce_sum over slices of sequences, so each chain uses several threads
// (compile with STAN_THREADS; threads per chain set at run time). Within a
// slice, the linear predictor is K matrix-vector products.

functions {
  real partial_loglik(array[] int y_slice, int start, int end,
                      array[] matrix X, array[] matrix M, vector w,
                      vector alpha, matrix B, vector gamma) {
    int n = end - start + 1;
    int K = rows(alpha);
    matrix[n, K + 1] eta;
    real lp = 0;
    for (k in 1:K) {
      eta[:, k] = alpha[k] + X[k][start:end] * B[k]' + M[k][start:end] * gamma;
    }
    eta[:, K + 1] = rep_vector(0, n);
    for (i in 1:n) {
      lp += w[start + i - 1] * (eta[i, y_slice[i]] - log_sum_exp(eta[i]));
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> P;
  array[N] int<lower=1, upper=K + 1> y;
  array[K] matrix[N, P] X;
  array[K] matrix[N, P] M;
  vector<lower=0>[N] w;
  real<lower=0> beta_prior_sd;
  real<lower=0> tau_prior_sd;
  real<lower=0> alpha_prior_sd;
  int<lower=1> grainsize;
}

parameters {
  real mu_alpha;
  real<lower=0> sigma_alpha;
  vector[K] alpha_raw;
  vector[P] beta;
  vector<lower=0>[P] tau;
  matrix[K, P] delta_raw;
  vector[P] gamma;
}

transformed parameters {
  vector[K] alpha = mu_alpha + sigma_alpha * alpha_raw;
  matrix[K, P] B = rep_matrix(beta', K) + diag_post_multiply(delta_raw, tau);
}

model {
  mu_alpha ~ normal(0, alpha_prior_sd);
  sigma_alpha ~ normal(0, alpha_prior_sd);
  alpha_raw ~ std_normal();
  beta ~ normal(0, beta_prior_sd);
  tau ~ normal(0, tau_prior_sd);
  to_vector(delta_raw) ~ std_normal();
  gamma ~ normal(0, 1);

  target += reduce_sum(partial_loglik, y, grainsize, X, M, w, alpha, B, gamma);
}
