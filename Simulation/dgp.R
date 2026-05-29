## ============================================================
## sim/dgp.R — data-generating process for the simulation study
## Truth is calibrated on the real M_full posteriors (truth_params.rds).
##
## simulate_dataset(truth, scenario, seed) returns a Stan-ready data
## list with KNOWN beta/sigma_u/ell so estimators can be scored.
##
## Scenarios (each targets one referee objection):
##   S1_jan  Jan-like: sigma_u ~ 0.024, OU kernel, Gaussian error.
##           -> beta recovery under extreme collinearity (|r|~0.99).
##   S2_jul  Jul-like: sigma_u ~ 0.52, OU kernel, Gaussian error.
##           -> GP identifiability when residual depth structure is real.
##   S3_nogp Jan-like but sigma_u = 0 exactly (no GP residual).
##           -> model must NOT manufacture depth structure.
##   S4_miss Jul-like but misspecified: squared-exponential GP kernel
##           (not OU) AND Student-t(3) measurement error.
##           -> robustness to kernel/error misspecification.
## ============================================================

suppressMessages(library(mvtnorm))

## OU / exponential (Matern-1/2) kernel — matches the fitted model.
kern_ou <- function(d, sigma_u, ell) {
  S <- length(d)
  D <- abs(outer(d, d, "-"))
  sigma_u^2 * exp(-D / ell)
}

## Squared-exponential kernel — used only for the misspecified scenario.
kern_se <- function(d, sigma_u, ell) {
  S <- length(d)
  D <- outer(d, d, "-")
  sigma_u^2 * exp(-0.5 * (D / ell)^2)
}

## Draw a zero-mean GP residual with jitter for Cholesky stability.
draw_gp <- function(d, sigma_u, ell, kernel = "ou") {
  S <- length(d)
  if (sigma_u <= 0) return(rep(0, S))
  K <- if (kernel == "se") kern_se(d, sigma_u, ell) else kern_ou(d, sigma_u, ell)
  diag(K) <- diag(K) + 1e-8
  as.numeric(rmvnorm(1, mean = rep(0, S), sigma = K))
}

## Heteroskedastic noise: Gaussian, or scaled Student-t (same sd).
add_noise <- function(mu, tau, dist = "gauss", df = 3) {
  if (dist == "t") {
    scale <- tau / sqrt(df / (df - 2))          # match sd to tau
    mu + scale * rt(length(mu), df = df)
  } else {
    mu + rnorm(length(mu), 0, tau)
  }
}

SCENARIOS <- list(
  S1_jan  = list(base = "jan", sigma_u_override = NULL, kernel = "ou", err = "gauss",
                 sigma_x_perturb = 0),
  S1b_pert= list(base = "jan", sigma_u_override = NULL, kernel = "ou", err = "gauss",
                 sigma_x_perturb = 0.20),       # M16: perturb Sigma_X eigenvalues by +/-20%
  S2_jul  = list(base = "jul", sigma_u_override = NULL, kernel = "ou", err = "gauss",
                 sigma_x_perturb = 0),
  S3_nogp = list(base = "jan", sigma_u_override = 0,    kernel = "ou", err = "gauss",
                 sigma_x_perturb = 0),
  S4_miss = list(base = "jul", sigma_u_override = NULL, kernel = "se", err = "t",
                 sigma_x_perturb = 0)
)

## Perturb Sigma_X by random eigenvalue rescaling within +/- p (multiplicative).
## Preserves positive-definiteness and the eigenvectors (i.e. preserves the
## qualitative correlation structure of the empirical covariance) but moves
## the spectrum so the simulation no longer matches the model's independent-X
## prior even after the marginalisation already present in the harness.
perturb_sigma_x <- function(Sigma, p, seed = NULL) {
  if (p <= 0) return(Sigma)
  if (!is.null(seed)) set.seed(seed)
  e <- eigen(Sigma, symmetric = TRUE)
  scl <- 1 + runif(length(e$values), -p, p)
  Sp <- e$vectors %*% diag(pmax(e$values * scl, 1e-8)) %*% t(e$vectors)
  (Sp + t(Sp)) / 2
}

simulate_dataset <- function(truth, scenario, seed) {
  set.seed(seed)
  sc <- SCENARIOS[[scenario]]
  if (is.null(sc)) stop("Unknown scenario: ", scenario)
  tr <- truth[[sc$base]]

  S <- tr$S; K <- tr$K
  beta    <- tr$beta
  sigma_u <- if (is.null(sc$sigma_u_override)) tr$sigma_u else sc$sigma_u
  ell     <- tr$ell
  d       <- tr$d
  tau_y   <- tr$tau_y
  tau_x   <- tr$tau_x                                   # S x K (plug-in, known)

  ## 1. Latent covariates with the REAL empirical collinearity.
  Sig <- tr$Sigma_X
  if (!is.null(sc$sigma_x_perturb) && sc$sigma_x_perturb > 0)
    Sig <- perturb_sigma_x(Sig, sc$sigma_x_perturb, seed = seed + 7919L)
  diag(Sig) <- diag(Sig) + 1e-8
  X_true <- rmvnorm(S, mean = tr$Xbar, sigma = Sig)     # S x K

  ## 2. GP residual depth structure.
  u <- draw_gp(d, sigma_u, ell, kernel = sc$kernel)

  ## 3. Latent response.
  y_true <- as.numeric(X_true %*% beta + u)

  ## 4. Heteroskedastic observation layer (EIV).
  X_obs <- X_true
  for (k in seq_len(K))
    X_obs[, k] <- add_noise(X_true[, k], tau_x[, k], dist = sc$err)
  y_obs <- add_noise(y_true, tau_y, dist = sc$err)

  list(
    stan = list(S = S, K = K,
                y_obs = as.vector(y_obs), X_obs = X_obs,
                tau_y = as.vector(tau_y), tau_x = tau_x,
                d = as.vector(d)),
    truth = list(beta = beta, sigma_u = sigma_u, ell = ell,
                 scenario = scenario, base = sc$base,
                 cov_names = tr$cov_names)
  )
}
