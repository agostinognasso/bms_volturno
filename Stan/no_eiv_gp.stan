// M_noEIV: GP regression without errors-in-variables correction on X.
// X_obs is treated as error-free; the GP depth structure is retained.
// The measurement error on y (tau_y) still enters as heteroskedastic
// observational noise on the diagonal of Sigma_y.
// This model differs from M_full by removing the latent X* layer.

functions {
  matrix ou_kernel(vector d, real sigma_u, real ell) {
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
  matrix[S, K] X_obs;           // treated as exact (no latent layer)
  vector<lower=0>[S] tau_y;     // known std dev of measurement error for y
  vector[S] d;                  // standardised depth in [0, 1]
}

parameters {
  vector[K] beta;
  real<lower=0> sigma_u;
  real<lower=0> ell;
}

transformed parameters {
  vector[S] mu_y = X_obs * beta;
  matrix[S, S] K_theta = ou_kernel(d, sigma_u, ell);
  matrix[S, S] Sigma_y;
  matrix[S, S] L_Sigma_y;

  Sigma_y = K_theta;
  for (s in 1:S)
    Sigma_y[s, s] = Sigma_y[s, s] + square(tau_y[s]) + 1e-8;

  L_Sigma_y = cholesky_decompose(Sigma_y);
}

model {
  beta    ~ normal(0, 0.5);
  sigma_u ~ normal(0, 1);
  ell     ~ lognormal(log(0.2), 0.5);

  y_obs ~ multi_normal_cholesky(mu_y, L_Sigma_y);
}

generated quantities {
  vector[S] y_rep;
  vector[S] log_lik;

  y_rep = multi_normal_cholesky_rng(mu_y, L_Sigma_y);

  {
    vector[S] v       = mdivide_left_tri_low(L_Sigma_y, y_obs - mu_y);
    real log_det      = 2 * sum(log(diagonal(L_Sigma_y)));
    real quad         = dot_self(v);
    real full_ll      = -0.5 * (S * log(2 * pi()) + log_det + quad);
    for (s in 1:S) log_lik[s] = full_ll / S;
  }
}
