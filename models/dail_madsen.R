# Dail-Madsen Model (Dail & Madsen 2011) - Open Population N-mixture
#
# Ecological process (initial abundance):
#   N_i1 ~ Poisson(lambda),  i = 1, ..., M
#
# Population dynamics (survival and recruitment):
#   S_it | N_it ~ Binomial(N_it, omega),  t = 1, ..., T-1
#   G_it ~ Poisson(gamma)
#   N_it+1 = S_it + G_it
#
# Observation process (detection):
#   y_it | N_it ~ Binomial(N_it, p),  t = 1, ..., T
#
# where:
#   lambda = expected initial abundance at each site
#   omega  = apparent survival probability
#   gamma  = recruitment rate (births and immigration, constant across sites)
#   p      = detection probability
#   N_it   = true (latent) abundance at site i, time t
#   S_it   = survivors at site i from time t to t+1
#   G_it   = recruits at site i from time t to t+1
#   y_it   = count at site i, time t
#
# N_it and S_it are latent random effects, marginalized jointly by RTMB via
# Sequential Reduction (SR) over {0, ..., K}:
#   L_i = sum_{N_i1} ... sum_{N_iT} sum_{S_i1} ... sum_{S_i,T-1}
#           Poisson(N_i1|lambda) * prod_t Binomial(S_it|N_it,omega) *
#           Poisson(N_it+1-S_it|gamma) * Binomial(y_it|N_it,p)
#
# JAGS fits the same model via MCMC, sampling N_it and S_it explicitly.
# Note: JAGS cannot use SR - it samples the full latent state space
# directly, which is substantially slower. This difference in wall-time is
# part of the story this repository tells.

library(RTMB)
library(unmarked)
library(R2jags)
source("R/utils.R")

# JAGS model as a string - written to a temp file at runtime
jags_model_dm <- "
model {
  # Priors
  lambda ~ dunif(0, 30)
  gamma  ~ dunif(0, 30)
  omega  ~ dunif(0, 1)
  p      ~ dunif(0, 1)

  # Likelihood
  for (i in 1:M) {
    # Initial abundance
    N[i, 1] ~ dpois(lambda)
    y[i, 1] ~ dbin(p, N[i, 1])

    for (t in 2:T) {
      # Survivors from t-1 to t
      S[i, t] ~ dbin(omega, N[i, t-1])
      # Recruits at t
      G[i, t] ~ dpois(gamma)
      # Total abundance at t
      N[i, t] <- S[i, t] + G[i, t]
      # Observation
      y[i, t] ~ dbin(p, N[i, t])
    }
  }
}
"

sim_dm <- function(M, T, lambda, gamma, omega, p) {
  N <- matrix(NA, M, T)
  N[, 1] <- rpois(M, lambda)
  for (t in 1:(T - 1)) {
    S <- rbinom(M, N[, t], omega)
    G <- rpois(M, gamma)
    N[, t + 1] <- S + G
  }
  matrix(rbinom(M * T, N, p), M, T)
}

fit_rtmb <- function(y) {
  M <- nrow(y)
  T <- ncol(y)
  K <- max(y) * 2
  dat <- list(y = y, M = M, T = T)

  f <- function(par) {
    "[<-" <- ADoverload("[<-")
    "c" <- ADoverload("c")
    getAll(par, dat)
    lambda <- exp(log_lambda)
    gamma <- exp(log_gamma)
    omega <- plogis(logit_omega)
    p <- plogis(logit_p)

    jnll <- -sum(dbinom(S, N[, 1:(T - 1)], omega, log = TRUE), na.rm = TRUE)
    jnll <- jnll - sum(dpois(N[, 1], lambda, log = TRUE), na.rm = TRUE)
    for (t in 1:(T - 1)) {
      G <- N[, t + 1] - S[, t]
      jnll <- jnll - sum(dpois(G, gamma, log = TRUE), na.rm = TRUE)
    }
    jnll <- jnll - sum(dbinom(y, size = N, prob = p, log = TRUE), na.rm = TRUE)
    jnll
  }

  par <- list(
    log_lambda  = log(mean(y[, 1]) + 0.1),
    log_gamma   = log(1.5),
    logit_omega = 0,
    logit_p     = 0,
    S = matrix(K, nrow = M, ncol = T - 1),
    N = matrix(K, nrow = M, ncol = T)
  )

  obj <- tryCatch(
    MakeADFun(f, par,
      random = c("N", "S"),
      integrate = list(
        S = TMB::SR(0:K, discrete = TRUE),
        N = TMB::SR(0:K, discrete = TRUE)
      ),
      silent = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(obj)) {
    return(c(lambda = NA, gamma = NA, omega = NA, p = NA))
  }

  opt <- tryCatch(
    nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 1e5, eval.max = 1e5)),
    error = function(e) NULL
  )
  if (is.null(opt) || opt$convergence != 0) {
    return(c(lambda = NA, gamma = NA, omega = NA, p = NA))
  }

  result <- c(
    exp(opt$par["log_lambda"]),
    exp(opt$par["log_gamma"]),
    plogis(opt$par["logit_omega"]),
    plogis(opt$par["logit_p"])
  )
  names(result) <- c("lambda", "gamma", "omega", "p")
  result
}

fit_unm <- function(y) {
  K <- max(y) * 2
  umf <- unmarkedFramePCO(y = y, numPrimary = ncol(y))
  fit <- tryCatch(
    suppressWarnings(pcountOpen(~1, ~1, ~1, ~1, data = umf, K = K, se = FALSE)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(c(lambda = NA, gamma = NA, omega = NA, p = NA))
  }

  result <- c(
    exp(coef(fit, "lambda")),
    exp(coef(fit, "gamma")),
    plogis(coef(fit, "omega")),
    plogis(coef(fit, "det"))
  )
  names(result) <- c("lambda", "gamma", "omega", "p")
  result
}

fit_jags_dm <- function(y, n.chains, n.iter, n.burnin, n.thin,
                        lambda_true) {
  M <- nrow(y)
  T <- ncol(y)

  model_file <- file.path(tempdir(), "jags_model_dm.txt")
  writeLines(jags_model_dm, model_file)

  # N[i,t] for t>=2 is a deterministic node (N <- S+G), cannot be initialised.
  # Only stochastic nodes need inits: N[i,1], S[i,t>=2], G[i,t>=2].
  N1_init <- pmax(apply(y, 1, max) + 2L, y[, 1] + 1L)
  S_init <- matrix(NA_integer_, M, T)
  G_init <- matrix(NA_integer_, M, T)
  N_prev <- N1_init
  for (t in 2:T) {
    S_init[, t] <- pmin(N_prev, round(N_prev * 0.8))
    G_init[, t] <- pmax(1L, y[, t] + 1L)
    N_prev <- S_init[, t] + G_init[, t]
  }

  jags_data <- list(y = y, M = M, T = T)
  jags_inits <- function() {
    list(
      lambda = runif(1, 1, lambda_true * 2),
      gamma  = runif(1, 0.5, 3),
      omega  = runif(1, 0.5, 0.99),
      p      = runif(1, 0.1, 0.9),
      N      = cbind(N1_init, matrix(NA_integer_, M, T - 1)),
      S      = S_init,
      G      = G_init
    )
  }

  fit <- tryCatch(
    suppressWarnings(
      jags(
        data = jags_data,
        inits = jags_inits,
        parameters.to.save = c("lambda", "gamma", "omega", "p"),
        model.file = model_file,
        n.chains = n.chains,
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.thin = n.thin,
        progress.bar = "none"
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(
      estimates      = c(lambda = NA, gamma = NA, omega = NA, p = NA),
      jags_converged = NA
    ))
  }

  sums <- fit$BUGSoutput$summary
  rhats <- sums[rownames(sums) != "deviance", "Rhat"]
  converged <- all(rhats < 1.1, na.rm = TRUE)

  list(
    estimates = c(
      lambda = sums["lambda", "mean"],
      gamma  = sums["gamma", "mean"],
      omega  = sums["omega", "mean"],
      p      = sums["p", "mean"]
    ),
    jags_converged = converged
  )
}

run_dailmadsen <- function(nsim = 1, M = 100, T = 5,
                           lambda_true = 4, gamma_true = 1.5,
                           omega_true = 0.8, p_true = 0.5,
                           seed = 333, mc.cores = 1,
                           n.chains = 3, n.iter = 5000,
                           n.burnin = 2500, n.thin = 1) {
  datasets <- lapply(1:nsim, function(s) {
    set.seed(seed + s)
    sim_dm(M, T, lambda_true, gamma_true, omega_true, p_true)
  })

  # Each dataset timed individually for every framework, so per-sim times
  # are always available and mean times are true per-fit times rather than
  # wall-clock-divided-by-nsim.
  rtmb_raw <- lapply_maybe(datasets, function(y) {
    t <- system.time({
      out <- tryCatch(fit_rtmb(y),
        error = function(e) c(lambda = NA, gamma = NA, omega = NA, p = NA)
      )
    })
    list(estimate = out, time = t["elapsed"])
  }, mc.cores = mc.cores)

  unm_raw <- lapply_maybe(datasets, function(y) {
    t <- system.time({
      out <- tryCatch(fit_unm(y),
        error = function(e) c(lambda = NA, gamma = NA, omega = NA, p = NA)
      )
    })
    list(estimate = out, time = t["elapsed"])
  }, mc.cores = mc.cores)

  jags_raw <- lapply_maybe(datasets, function(y) {
    t <- system.time({
      out <- fit_jags_dm(y,
        n.chains = n.chains,
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.thin = n.thin,
        lambda_true = lambda_true
      )
    })
    list(
      estimates = out$estimates,
      jags_converged = out$jags_converged,
      time_jags = t["elapsed"]
    )
  }, mc.cores = mc.cores)

  res_rtmb <- as.data.frame(do.call(rbind, lapply(rtmb_raw, function(x) x$estimate)))
  res_unm <- as.data.frame(do.call(rbind, lapply(unm_raw, function(x) x$estimate)))
  res_jags <- as.data.frame(
    do.call(rbind, lapply(jags_raw, function(x) x$estimates))
  )
  names(res_rtmb) <- c("lambda", "gamma", "omega", "p")
  names(res_unm) <- c("lambda", "gamma", "omega", "p")
  names(res_jags) <- c("lambda", "gamma", "omega", "p")

  rtmb_times <- sapply(rtmb_raw, function(x) x$time)
  unm_times <- sapply(unm_raw, function(x) x$time)
  jags_times <- sapply(jags_raw, function(x) x$time_jags)

  # Drop the same simulation indices from all three frameworks so estimates
  # always correspond to the same dataset - critical for fair comparison
  ok <- complete.cases(res_rtmb) & complete.cases(res_unm) & complete.cases(res_jags)
  res_rtmb <- res_rtmb[ok, ]
  res_unm <- res_unm[ok, ]
  res_jags <- res_jags[ok, ]

  conv_rate <- mean(sapply(jags_raw, function(x) x$jags_converged), na.rm = TRUE)

  list(
    model = "dailmadsen",
    estimates_rtmb = res_rtmb,
    estimates_unm = res_unm,
    estimates_jags = res_jags,
    times = data.frame(
      rtmb = sum(rtmb_times, na.rm = TRUE),
      unm  = sum(unm_times, na.rm = TRUE),
      jags = sum(jags_times, na.rm = TRUE)
    ),
    rtmb_times_per_sim = rtmb_times,
    unm_times_per_sim = unm_times,
    jags_times_per_sim = jags_times,
    jags_conv_rate = conv_rate,
    truth = c(
      lambda = lambda_true, gamma = gamma_true,
      omega = omega_true, p = p_true
    ),
    nsim = nsim, M = M, T = T
  )
}

# Run standalone if called directly
if (sys.nframe() == 0) {
  res <- run_dailmadsen()
  dir.create("results", showWarnings = FALSE)
  saveRDS(res, "results/dail_madsen.rds")
  cat("Total RTMB time:", res$times$rtmb, "s\n")
  cat("Total unmarked time:", res$times$unm, "s\n")
  cat("Total JAGS time:", res$times$jags, "s\n")
  cat("JAGS convergence rate:", round(res$jags_conv_rate, 3), "\n")
}
