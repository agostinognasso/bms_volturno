// M_noGP: EIV regression without Gaussian-process depth dependence.
// The latent X* layer is retained; sigma_u is removed (Sigma_y = D_y).
// Residuals are independent across depth bins (heteroskedastic normal).
// This model differs from M_full by setting K_theta = 0.

data {
  int<lower=1> S;
  int<lower=1> K;
  vector[S] y_obs;
  matrix[S, K] X_obs;
  vector<lower=0>[S] tau_y;
  matrix<lower=0>[S, K] tau_x;
  vector[S] d;                  // unused structurally; kept for interface parity
}

parameters {
  vector[K] beta;
  vector[K] mu_x;
  vector<lower=0>[K] sigma_x;
  matrix[S, K] Z_x;             // non-centred latent covariates
}

transformed parameters {
  matrix[S, K] X_true;
  vector[S] mu_y;

  for (k in 1:K)
    X_true[, k] = mu_x[k] + sigma_x[k] * Z_x[, k];

  mu_y = X_true * beta;
}

model {
  beta    ~ normal(0, 0.5);
  mu_x    ~ normal(0, 1);
  sigma_x ~ normal(0, 1);
  to_vector(Z_x) ~ normal(0, 1);

  // Measurement model for X
  for (s in 1:S)
    for (k in 1:K)
      X_obs[s, k] ~ normal(X_true[s, k], tau_x[s, k]);

  // Independent heteroskedastic likelihood (no GP covariance)
  y_obs ~ normal(mu_y, tau_y);
}

generated quantities {
  vector[S] y_rep;
  vector[S] log_lik;

  for (s in 1:S) {
    y_rep[s]    = normal_rng(mu_y[s], tau_y[s]);
    log_lik[s]  = normal_lpdf(y_obs[s] | mu_y[s], tau_y[s]);
  }
}
