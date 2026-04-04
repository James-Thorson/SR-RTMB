# =============================================================
# dail-madsen + spatial GMRF random effect on lambda via
# SPDE approach used in spatio-temporal models for ecologists
# =============================================================
#
# N_i1 ~ Poisson(lambda_i),  log(lambda_i) = mu_lambda + omega_i
# S_it ~ Binomial(N_it, omega)
# G_it ~ Poisson(gamma)
# N_i,t+1 = S_it + G_it
# y_it ~ Binomial(N_it, p)
# omega_s ~ GMRF(0, Q / tau^2)

library(RTMB)
library(fmesher)

set.seed(1234)

M <- 150
T <- 5
K <- 30

# true parameters
mu_lambda_true <- log(3.5)
omega_true <- 0.7
gamma_true <- 1.5
p_true <- 0.4
ln_tauO_true <- log(0.2)
ln_kappa_true <- log(sqrt(8) / 0.3)

# -------------------------------------------------------------
# 1. simulation
# -------------------------------------------------------------

coords <- matrix(runif(M * 2), ncol = 2, dimnames = list(NULL, c("x", "y")))

mesh <- fm_mesh_2d(coords, cutoff = 0.05, refine = list(max.edge = c(0.1, 0.3)))
spde <- fm_fem(mesh, order = 2)
A_is <- fm_evaluator(mesh, loc = coords)$proj$A

Q_sim <- exp(4 * ln_kappa_true) * spde$c0 +
  2 * exp(2 * ln_kappa_true) * spde$g1 +
  spde$g2

# draw spatial field on mesh vertices
L <- chol(as.matrix(Q_sim))
omega_v <- backsolve(L, rnorm(mesh$n)) / exp(ln_tauO_true)
omega_i <- as.numeric(A_is %*% omega_v)
lambda_i <- exp(mu_lambda_true + omega_i)

N <- matrix(NA, M, T)
N[, 1] <- rpois(M, lambda_i)
for (t in 1:(T - 1)) {
  N[, t + 1] <- rbinom(M, N[, t], omega_true) + rpois(M, gamma_true)
}
y <- matrix(rbinom(M * T, N, p_true), M, T)

# -------------------------------------------------------------
# 2. plot: mesh and site locations
# -------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)

png("figures/spde_mesh.png", width = 600, height = 600, res = 150)
par(mar = c(4, 4, 3, 1))
plot(mesh,
  main = "mesh and site locations",
  edge.color = "grey70", col = "grey90"
)
points(coords, pch = 19, cex = 0.7, col = "#2166ac")
dev.off()

# -------------------------------------------------------------
# 3. fitting
# -------------------------------------------------------------

dat <- list(
  y = y, M = M, T = T, K = K, ks = 0:K,
  A_is = A_is, M0 = spde$c0, M1 = spde$g1, M2 = spde$g2
)

par <- list(
  mu_lambda   = log(mean(y[, 1]) + 0.1),
  log_gamma   = log(1),
  logit_omega = 0,
  logit_p     = 0,
  ln_tauO     = log(1),
  ln_kappa    = log(1),
  omega_s     = rep(0, mesh$n)
)

f <- function(par) {
  getAll(dat, par, warn = FALSE)
  gamma <- exp(log_gamma)
  omega <- plogis(logit_omega)
  p <- plogis(logit_p)
  jnll <- 0
  Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
  jnll <- jnll - dgmrf(omega_s,
    mu = 0, Q = Q, log = TRUE,
    scale = 1 / exp(ln_tauO)
  )
  lambda_i <- exp(mu_lambda + A_is %*% omega_s)

  # transition matrix -- fixed, built once
  trans <- matrix(0, K + 1, K + 1)
  for (n in 0:K) {
    sp <- dbinom(0:n, n, omega)
    rp <- dpois(0:K, gamma)
    for (m in 0:K) {
      for (s in 0:min(n, m)) {
        trans[n + 1, m + 1] <- trans[n + 1, m + 1] + sp[s + 1] * rp[m - s + 1]
      }
    }
  }

  for (i in 1:M) {
    alpha <- dpois(ks, lambda_i[i]) * dbinom(y[i, 1], ks, p)
    for (t in 2:T) {
      alpha <- as.vector(t(trans) %*% alpha) * dbinom(y[i, t], ks, p)
    }
    jnll <- jnll - log(sum(alpha))
  }
  jnll
}

obj <- MakeADFun(f, par, random = "omega_s")
opt <- nlminb(obj$par, obj$fn, obj$gr,
  control = list(eval.max = 1e4, iter.max = 1e4)
)
sdr <- sdreport(obj)

# -------------------------------------------------------------
# 4. results
# -------------------------------------------------------------

est <- summary(sdr, "fixed")

print(data.frame(
  parameter = c("mu_lambda", "gamma", "omega", "p", "ln_tauO", "ln_kappa"),
  truth = c(
    mu_lambda_true, gamma_true, omega_true, p_true,
    ln_tauO_true, ln_kappa_true
  ),
  estimate = c(
    est["mu_lambda", "Estimate"],
    exp(est["log_gamma", "Estimate"]),
    plogis(est["logit_omega", "Estimate"]),
    plogis(est["logit_p", "Estimate"]),
    est["ln_tauO", "Estimate"],
    est["ln_kappa", "Estimate"]
  )
), digits = 3)

# -------------------------------------------------------------
# 5. plot: estimated vs true spatial field
# -------------------------------------------------------------

omega_est_sites <- as.numeric(A_is %*% obj$env$parList()$omega_s)

lim <- max(abs(c(omega_i, omega_est_sites)))
pal <- colorRampPalette(c("#2166ac", "#f7f7f7", "#d6604d"))(200)
col_fn <- function(x) {
  pal[findInterval(x, seq(-lim, lim, length.out = 201),
    all.inside = TRUE
  )]
}

png("figures/spde_field.png", width = 900, height = 450, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(coords,
  col = col_fn(omega_i), pch = 19, cex = 1, asp = 1,
  xlab = "x", ylab = "y", main = expression(lambda[i] ~ "true")
)

plot(coords,
  col = col_fn(omega_est_sites), pch = 19, cex = 1, asp = 1,
  xlab = "x", ylab = "y", main = expression(lambda[i] ~ "estimated")
)
dev.off()
