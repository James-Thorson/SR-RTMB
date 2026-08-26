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
#
# omega_s is continuous, integrated out via RTMB's Laplace approximation;
# N_it and S_it are discrete, marginalized jointly via Sequential Reduction
# (SR) over {0, ..., K} - the same redux state-space formulation used in
# models/dail_madsen.R.

# Make sure to use github version of RTMB
# pak::pak("kaskr/RTMB/RTMB")

library(RTMB)
library(fmesher)
library(ggplot2)
library(patchwork)

dir.create("plots", showWarnings = FALSE)

set.seed(1234)
M <- 200
T <- 3

# true parameters
mu_lambda_true <- log(3.5)
omega_true <- 0.7
gamma_true <- 1.5
p_true <- 0.4
omega_range <- 0.3
omega_SD <- 1
ln_kappa_true <- log(sqrt(8) / omega_range)
ln_tauO_true <- log(1 / (omega_SD * exp(ln_kappa_true) * sqrt(4*pi)))

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
K <- max(y) * 2
if(K > 30) stop("check")

# -------------------------------------------------------------
# 2. plot: mesh and site locations
# -------------------------------------------------------------

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
  y = y, M = M, T = T,
  A_is = A_is, M0 = spde$c0, M1 = spde$g1, M2 = spde$g2
)

par <- list(
  mu_lambda = log(mean(y[, 1]) + 0.1),
  log_gamma = log(1),
  logit_omega = 0,
  logit_p = 0,
  ln_tauO = log(1),
  ln_kappa = log(1),
  omega_s = rep(0, mesh$n),
  S = matrix(K, nrow = M, ncol = T - 1),
  N = matrix(K, nrow = M, ncol = T)
)

f <- function(par) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(dat, par, warn = FALSE)
  gamma <- exp(log_gamma)
  omega <- plogis(logit_omega)
  p <- plogis(logit_p)

  Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
  jnll <- -dgmrf(omega_s,
    mu = 0, Q = Q, log = TRUE,
    scale = 1 / exp(ln_tauO)
  )
  lambda_i <- exp(mu_lambda + (A_is %*% omega_s)[,1])

  jnll <- jnll - sum(dbinom(S, N[, 1:(T - 1)], omega, log = TRUE), na.rm = TRUE)
  jnll <- jnll - sum(dpois(N[, 1], lambda_i, log = TRUE), na.rm = TRUE)
  for (t in 1:(T - 1)) {
    G <- N[, t + 1] - S[, t]
    jnll <- jnll - sum(dpois(G, gamma, log = TRUE), na.rm = TRUE)
  }
  jnll <- jnll - sum(dbinom(y, size = N, prob = p, log = TRUE), na.rm = TRUE)
  jnll
}
f(par)
K

obj <- MakeADFun(f, par,
  random = c("omega_s", "N", "S"),
  integrate = list(
    S = TMB::SR(0:K, discrete = TRUE),
    N = TMB::SR(0:K, discrete = TRUE)
  ),
  silent = TRUE
)
opt <- nlminb(obj$par, obj$fn, obj$gr,
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
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
mesh_sfc <- fm_as_sfc(mesh)

field_dat <- rbind(
  data.frame(x = coords[, "x"], y = coords[, "y"], omega = omega_i, panel = "true"),
  data.frame(x = coords[, "x"], y = coords[, "y"], omega = omega_est_sites, panel = "estimated")
)

plot_field <- function(df, title) {
#  ggplot(df, aes(x = x, y = y, color = omega)) +
#    geom_point(size = 2.2) +
#    scale_color_gradient2(
#      low = "#2166ac", mid = "white", high = "#b2182b",
#      midpoint = 0, limits = c(-lim, lim), name = expression(omega[i])
#    ) +
#    coord_equal() +
#    labs(title = title, x = "x", y = "y") +
#    theme_minimal(base_size = 11) +
#    theme(panel.grid.minor = element_blank())
  ggplot(df, aes(x = x, y = y)) +
    # mesh triangles drawn first, underneath the points
    geom_sf(
      data = mesh_sfc,
      inherit.aes = FALSE,
      fill = NA,
      color = "grey70",
      linewidth = 0.15
    ) +
    geom_point(aes(color = omega), size = 2.2) +
    scale_color_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, limits = c(-lim, lim), name = expression(omega[i])
    ) +
    coord_sf(
      #xlim = range(df$x), ylim = range(df$y),
      default_crs = sf::st_crs(mesh_sfc), expand = TRUE
    ) +
    labs(title = title, x = "x", y = "y") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

fig_field <- plot_field(field_dat[field_dat$panel == "true", ], expression(omega[i] ~ "true")) +
  plot_field(field_dat[field_dat$panel == "estimated", ], expression(omega[i] ~ "estimated")) +
  plot_layout(guides = "collect")

ggsave("figures/spde_field.png", fig_field, width = 9, height = 4.5, dpi = 150)
