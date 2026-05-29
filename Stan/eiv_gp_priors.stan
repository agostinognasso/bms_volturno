// eiv_gp_priors.stan
// Same model as eiv_gp.stan, but exposes the scale parameters of
// the beta and sigma_u priors as data inputs (prior_s_beta and
// prior_s_sigma_u). Used for the prior-sensitivity analysis
// (Section "Prior sensitivity" of the revised manuscript).

functions {
  matrix cov_exp_quad_depth(vector d, real sigma_u, real ell) {
    int S = num_elements(d);
    matrix[S, S] K;
    for (i in 1:S) {
      for (j in i:S) {
        real val = square(sigma_u) * exp(-fabs(d[i] - d[j]) / ell);
        K[i, j] = val;
        K[j, i] = val;
      }
    }
    return K;
  }
}

data {
  int<lower=1> S;
  int<lower=1> K;
  vector[S] y_obs;
  matrix[S, K] X_obs;
  vector<lower=0>[S] tau_y;
  matrix<lower=0>[S, K] tau_x;
  vector[S] d;
  real<lower=0> prior_s_beta;
  real<lower=0> prior_s_sigma_u;
}

parameters {
  vector[K] beta;
  vector[K] mu_x;
  vector<lower=0>[K] sigma_x;
  matrix[S, K] Z_x;
  real<lower=0> sigma_u;
  real<lower=0> ell;
}

transformed parameters {
  matrix[S, K] X_true;
  vector[S] mu_y;
  matrix[S, S] K_theta;
  matrix[S, S] Sigma_y;
  matrix[S, S] L_Sigma_y;
  for (k in 1:K) X_true[, k] = mu_x[k] + sigma_x[k] * Z_x[, k];
  mu_y = X_true * beta;
  K_theta = cov_exp_quad_depth(d, sigma_u, ell);
  Sigma_y = K_theta;
  for (s in 1:S) Sigma_y[s, s] = Sigma_y[s, s] + square(tau_y[s]) + 1e-8;
  L_Sigma_y = cholesky_decompose(Sigma_y);
}

model {
  beta    ~ normal(0, prior_s_beta);
  mu_x    ~ normal(0, 1);
  sigma_x ~ normal(0, 1);
  to_vector(Z_x) ~ normal(0, 1);
  sigma_u ~ normal(0, prior_s_sigma_u);
  ell     ~ lognormal(log(0.2), 0.5);
  for (s in 1:S) for (k in 1:K)
    X_obs[s, k] ~ normal(X_true[s, k], tau_x[s, k]);
  y_obs ~ multi_normal_cholesky(mu_y, L_Sigma_y);
}
