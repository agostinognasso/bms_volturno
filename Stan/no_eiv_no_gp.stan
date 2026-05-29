// M_base: classical weighted least squares — no EIV, no GP.
// X_obs treated as error-free; residuals are independent and heteroskedastic
// with known variance tau_y^2.  This is the baseline model.

data {
  int<lower=1> S;
  int<lower=1> K;
  vector[S] y_obs;
  matrix[S, K] X_obs;
  vector<lower=0>[S] tau_y;
  vector[S] d;                  // unused; kept for interface parity
}

parameters {
  vector[K] beta;
}

transformed parameters {
  vector[S] mu_y = X_obs * beta;
}

model {
  beta ~ normal(0, 0.5);
  y_obs ~ normal(mu_y, tau_y);
}

generated quantities {
  vector[S] y_rep;
  vector[S] log_lik;

  for (s in 1:S) {
    y_rep[s]   = normal_rng(mu_y[s], tau_y[s]);
    log_lik[s] = normal_lpdf(y_obs[s] | mu_y[s], tau_y[s]);
  }
}
