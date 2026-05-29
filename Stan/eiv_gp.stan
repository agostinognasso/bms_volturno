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
  int<lower=1> S;                       // number of depth bins
  int<lower=1> K;                       // number of covariates
  vector[S] y_obs;                      // standardised observed response
  matrix[S, K] X_obs;                   // standardised observed covariates
  vector<lower=0>[S] tau_y;             // std. dev. of measurement error for y
  matrix<lower=0>[S, K] tau_x;          // std. dev. of measurement error for X
  vector[S] d;                          // standardised depth locations in [0,1]
}

parameters {
  // beta0 is dropped: y and all X are standardised to mean 0 before fitting,
  // so the intercept is absorbed by construction and set to zero.
  vector[K] beta;

  vector[K] mu_x;
  vector<lower=0>[K] sigma_x;

  matrix[S, K] Z_x;                     // non-centred latent covariates

  real<lower=0> sigma_u;
  real<lower=0> ell;
}

transformed parameters {
  matrix[S, K] X_true;
  vector[S] mu_y;
  matrix[S, S] K_theta;
  matrix[S, S] Sigma_y;
  matrix[S, S] L_Sigma_y;

  // Non-centred parameterisation for latent X*
  for (k in 1:K) {
    X_true[, k] = mu_x[k] + sigma_x[k] * Z_x[, k];
  }

  // Structural mean (GP residual u is marginalised out into Sigma_y)
  mu_y = X_true * beta;

  // GP covariance kernel (Ornstein-Uhlenbeck / Matern-1/2)
  K_theta = cov_exp_quad_depth(d, sigma_u, ell);

  // Marginal covariance: GP structure + heteroskedastic measurement error
  // Sigma_y = K_theta + diag(tau_y^2)
  // The measurement error tau_y[s] enters additively on the diagonal,
  // separating observational noise from structural depth dependence.
  Sigma_y = K_theta;
  for (s in 1:S) {
    Sigma_y[s, s] = Sigma_y[s, s] + square(tau_y[s]) + 1e-8;
  }

  L_Sigma_y = cholesky_decompose(Sigma_y);
}

model {
  // Weakly informative priors on standardised scale
  beta    ~ normal(0, 0.5);
  mu_x    ~ normal(0, 1);
  sigma_x ~ normal(0, 1);             // half-normal by lower bound
  to_vector(Z_x) ~ normal(0, 1);

  sigma_u ~ normal(0, 1);             // half-normal by lower bound
  ell     ~ lognormal(log(0.2), 0.5); // depth standardised to [0,1]

  // Measurement model for X (errors-in-variables layer)
  for (s in 1:S) {
    for (k in 1:K) {
      X_obs[s, k] ~ normal(X_true[s, k], tau_x[s, k]);
    }
  }

  // Marginal GP likelihood for y (u already integrated out into Sigma_y)
  y_obs ~ multi_normal_cholesky(mu_y, L_Sigma_y);
}

generated quantities {
  vector[S] y_rep;
  vector[S] log_lik;
  corr_matrix[S] Corr_y;

  // Posterior predictive draw from the full multivariate distribution.
  // This correctly propagates both the GP covariance structure and the
  // heteroskedastic measurement error (both encoded in L_Sigma_y).
  y_rep = multi_normal_cholesky_rng(mu_y, L_Sigma_y);

  // Posterior correlation structure implied by Sigma_y.
  // Computed directly from Sigma_y — no matrix inversion required.
  for (s in 1:S) {
    Corr_y[s, s] = 1;
    for (t in (s + 1):S) {
      Corr_y[s, t] = Sigma_y[s, t] / sqrt(Sigma_y[s, s] * Sigma_y[t, t]);
      Corr_y[t, s] = Corr_y[s, t];
    }
  }

  // Log-likelihood: full joint log p(y_obs | X*, beta, sigma_u, ell),
  // split equally across bins as a convenience scalar.
  // NOTE: log_lik[s] = full_ll / S is NOT a proper pointwise LOO log-lik
  // for correlated data. True LOO requires conditional predictive densities
  // p(y_s | y_{-s}, ...). Use k-fold CV for proper out-of-sample evaluation.
  // The quad form is computed via Cholesky forward-solve (O(S^2))
  // instead of explicit matrix inversion (O(S^3)).
  {
    vector[S] v = mdivide_left_tri_low(L_Sigma_y, y_obs - mu_y);
    real log_det_Sigma = 2 * sum(log(diagonal(L_Sigma_y)));
    real quad    = dot_self(v);
    real full_ll = -0.5 * (S * log(2 * pi()) + log_det_Sigma + quad);
    for (s in 1:S) log_lik[s] = full_ll / S;
  }
}
