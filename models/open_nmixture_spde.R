# -------------------
# open N-mixture (Dail-Madsen) + spatial GMRF random effect on lambda
# via SPDE, self-test
# -------------------
#
# N_i1 ~ Poisson(lambda_i),  log(lambda_i) = mu_lambda + epsilon_i
# S_it ~ Binomial(N_it, omega)
# G_it ~ Poisson(gamma)
# N_i,t+1 = S_it + G_it
# y_it ~ Binomial(N_it, p)
# epsilon_s ~ GMRF(0, Q / tau^2)
#
# epsilon_s is continuous, integrated out via RTMB's Laplace approximation;
# N_it and S_it are discrete, marginalized jointly via Sequential Reduction
# (SR) over {0, ..., K}
#-------------------

library(RTMB)
library(fmesher)
library(ggplot2)
library(patchwork)

# run parameters
# set by `make NSIM=X`; set nsim manually here instead if running standalone
nsim <- as.integer(Sys.getenv("NSIM"))
# set by `make SEED=X`; set seed manually here instead if running standalone
seed <- as.integer(Sys.getenv("SEED"))
set.seed(seed)
M <- 200
T <- 4

# true parameters
mu_lambda_true <- log(2.5)
omega_true <- 0.7
gamma_true <- 1.5
p_true <- 0.4
epsilon_range <- 0.3
epsilon_SD <- 0.5
ln_kappa_true <- log(sqrt(8) / epsilon_range)
ln_tau_true <- log(1 / (epsilon_SD * exp(ln_kappa_true) * sqrt(4 * pi)))

truth <- c(
  mu_lambda = mu_lambda_true, gamma = gamma_true, omega = omega_true,
  p = p_true, ln_tau = ln_tau_true, ln_kappa = ln_kappa_true
)

res <- matrix(NA_real_,
  nrow = nsim, ncol = length(truth),
  dimnames = list(NULL, names(truth))
)
converged <- rep(FALSE, nsim)
elapsed_sec <- rep(NA, nsim)
# Snapshot of the fields the post-loop spatial-field plot needs, taken from
# whichever replicate last converged - NOT whatever the loop variables
# happen to hold when it ends, since the final replicate can still fail to
# fit (MakeADFun error or non-convergence) and leave A_is/mesh mismatched
# with obj/epsilon_i from an earlier replicate.
last_fit <- NULL

for (s in 1:nsim) {
  cat("replicate", s, "of", nsim, "\n")

  # -------------------------------------------------------------
  # simulate: new site locations, mesh, and spatial field each rep
  # -------------------------------------------------------------
  coords <- matrix(runif(M * 2),
    ncol = 2,
    dimnames = list(NULL, c("x", "y"))
  )
  mesh <- fm_mesh_2d(coords,
    cutoff = 0.05,
    refine = list(max.edge = c(0.1, 0.3))
  )
  spde <- fm_fem(mesh, order = 2)
  A_is <- fm_evaluator(mesh, loc = coords)$proj$A

  Q_sim <- exp(4 * ln_kappa_true) * spde$c0 +
    2 * exp(2 * ln_kappa_true) * spde$g1 +
    spde$g2

  L <- chol(as.matrix(Q_sim))
  epsilon_v <- backsolve(L, rnorm(mesh$n)) / exp(ln_tau_true)
  epsilon_i <- as.numeric(A_is %*% epsilon_v)
  lambda_i <- exp(mu_lambda_true + epsilon_i)
  gamma_i <- gamma_true * exp(epsilon_i)

  N <- matrix(NA, M, T)
  N[, 1] <- rpois(M, lambda_i)
  for (t in 1:(T - 1)) {
    N[, t + 1] <- rbinom(M, N[, t], omega_true) + rpois(M, gamma_i)
  }
  y <- matrix(rbinom(M * T, N, p_true), M, T)
  K <- max(y) * 3

  # leave Jim's control look in here for now...
  # if (K > 30) {
  #   warning("replicate ", s, ": K = ", K, " too large, skipping")
  #   next
  # }

  # -------------------------------------------------------------
  # fit
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
    ln_tau = log(1),
    ln_kappa = log(1),
    epsilon_s = rep(0, mesh$n),
    SN = matrix(K, nrow = M, ncol = 2 * T - 1)
  )

  f <- function(par) {
    getAll(dat, par, warn = FALSE)
    S <- SN[, seq_len(T - 1)]
    N <- SN[, T - 1 + seq_len(T)]
    gamma <- exp(log_gamma)
    omega <- plogis(logit_omega)
    p <- plogis(logit_p)

    Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
    jnll <- -dgmrf(epsilon_s,
      mu = 0, Q = Q, log = TRUE,
      scale = 1 / exp(ln_tau)
    )
    lambda_i <- exp(mu_lambda + (A_is %*% epsilon_s)[, 1])
    gamma_i <- gamma * exp(A_is %*% epsilon_s)[, 1]
    jnll <- jnll - sum(dbinom(S, N[, 1:(T - 1)], omega, log = TRUE), na.rm = TRUE)
    jnll <- jnll - sum(dpois(N[, 1], lambda_i, log = TRUE), na.rm = TRUE)
    for (t in 1:(T - 1)) {
      G <- N[, t + 1] - S[, t]
      jnll <- jnll - sum(dpois(G, gamma_i, log = TRUE), na.rm = TRUE)
    }
    jnll <- jnll - sum(dbinom(y, size = N, prob = p, log = TRUE), na.rm = TRUE)
    jnll
  }

  t0 <- Sys.time()

  obj <- tryCatch(
    MakeADFun(f, par,
      random = c("epsilon_s", "SN"),
      integrate = list(
        SN = TMB::SR(0:K, discrete = TRUE)
      ),
      silent = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(obj)) next

  opt <- tryCatch(
    nlminb(obj$par, obj$fn, obj$gr,
      control = list(eval.max = 1e4, iter.max = 1e4)
    ),
    error = function(e) NULL
  )

  elapsed_sec[s] <- as.numeric(Sys.time() - t0, units = "secs")
  if (is.null(opt) || opt$convergence != 0) next
  converged[s] <- TRUE

  res[s, "mu_lambda"] <- opt$par["mu_lambda"]
  res[s, "gamma"] <- exp(opt$par["log_gamma"])
  res[s, "omega"] <- plogis(opt$par["logit_omega"])
  res[s, "p"] <- plogis(opt$par["logit_p"])
  res[s, "ln_tau"] <- opt$par["ln_tau"]
  res[s, "ln_kappa"] <- opt$par["ln_kappa"]

  last_fit <- list(
    obj = obj, A_is = A_is, mesh = mesh,
    coords = coords, epsilon_i = epsilon_i
  )
}

# -------------------------------------------------------------
# results
# -------------------------------------------------------------
results <- data.frame(sim = 1:nsim, res)
for (param in names(truth)) {
  results[[paste0(param, "_true")]] <- truth[[param]]
}
results$converged <- converged
results$elapsed_sec <- elapsed_sec

saveRDS(results, "results/open_nmixture_spde_sim.rds")
# results <- readRDS("results/open_nmixture_spde_sim.rds")

print(results, digits = 3, row.names = FALSE)
cat(sum(converged), "of", nsim, "replicates converged\n")

results <- results[results$converged, ]

# -------------------------------------------------------------
# plot: estimated vs true spatial field (last converged replicate)
# -------------------------------------------------------------

epsilon_est_sites <- as.numeric(last_fit$A_is %*% last_fit$obj$env$parList()$epsilon_s)
lim <- max(abs(c(last_fit$epsilon_i, epsilon_est_sites)))
mesh_sfc <- fm_as_sfc(last_fit$mesh)

field_dat <- rbind(
  data.frame(
    x = last_fit$coords[, "x"], y = last_fit$coords[, "y"],
    omega = last_fit$epsilon_i, panel = "true"
  ),
  data.frame(
    x = last_fit$coords[, "x"], y = last_fit$coords[, "y"],
    omega = epsilon_est_sites, panel = "estimated"
  )
)

plot_field <- function(df, title) {
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
      default_crs = sf::st_crs(mesh_sfc), expand = TRUE
    ) +
    labs(title = title, x = "x", y = "y") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

fig_field <- plot_field(
  field_dat[field_dat$panel == "true", ],
  expression(omega[i] ~ "true")
) +
  plot_field(
    field_dat[field_dat$panel == "estimated", ],
    expression(omega[i] ~ "estimated")
  ) +
  plot_layout(guides = "collect")


# -------------------------------------------------------------
# violin plot: relative bias of estimates vs. truth, by parameter
# (style matches R/plotting.R::plot_estimate_recovery)
# -------------------------------------------------------------

relbias_long <- do.call(rbind, lapply(names(truth), function(param) {
  data.frame(
    parameter = param,
    rel_bias = (results[[param]] - results[[paste0(param, "_true")]]) /
      results[[paste0(param, "_true")]]
  )
}))
# Ordering mirrors "Open N-mixture" in R/plotting.R::plot_estimate_recovery
# (lambda, p, gamma, omega), with the spatial-field params new to this
# script (ln_tau, ln_kappa) appended after.
relbias_long$parameter <- factor(relbias_long$parameter,
  levels = c("mu_lambda", "p", "gamma", "omega", "ln_tau", "ln_kappa")
)

fig_bias <- ggplot(relbias_long, aes(x = parameter, y = rel_bias)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_violin(
    fill = "#8E44AD", alpha = 0.85, linewidth = 0.3,
    trim = FALSE, na.rm = TRUE
  ) +
  geom_boxplot(
    width = 0.12, outlier.shape = NA, fill = "white",
    alpha = 0.9, na.rm = TRUE
  ) +
  scale_x_discrete(labels = c(
    mu_lambda = expression(lambda),
    p = expression(italic(p)),
    gamma = expression(gamma),
    omega = expression(omega),
    # these are ln_tau/ln_kappa, not tau/kappa themselves
    ln_tau = expression(log(tau)),
    ln_kappa = expression(log(kappa))
  )) +
  labs(x = NULL, y = "Relative bias  (estimate - truth) / truth") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())


# -------------------------------------------------------------
# violin plot: per-fit run time
# -------------------------------------------------------------

fig_time <- ggplot(results, aes(x = "", y = elapsed_sec, fill = "RTMB")) +
  geom_violin(alpha = 0.85, linewidth = 0.3, trim = FALSE, na.rm = TRUE) +
  geom_boxplot(
    width = 0.12, outlier.shape = NA, fill = "white",
    alpha = 0.9, na.rm = TRUE
  ) +
  scale_fill_manual(values = c(RTMB = "#8E44AD"), name = NULL) +
  scale_y_log10() +
  labs(x = NULL, y = "Time per fit (s, log scale)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# -------------------------------------------------------------
# combined: truth vs. estimates | timing / field
# -------------------------------------------------------------

fig_all <- fig_field / (fig_bias | fig_time) &
  theme(aspect.ratio = 1)


ggsave("figures/open_nmixture_spde.png",
  fig_all,
  width = 10,
  height = 8,
  dpi = 300
)
