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
# Marginal likelihood integrates out N_it via forward algorithm:
#   alpha_t(m) = P(y_i1,...,y_it, N_it=m)
#             = P(y_it|N_it=m) * sum_n P(N_it=m|N_it-1=n) * alpha_t-1(n)
#   L_i = sum_m alpha_T(m)
#
# JAGS fits the same model via MCMC, sampling N_it and S_it explicitly.
# Note: JAGS cannot use the forward algorithm - it samples the full latent
# state space directly, which is substantially slower. This difference in
# wall-time is part of the story this repository tells.

library(RTMB)
library(unmarked)
library(R2jags)

# JAGS model as a string - written to a temp file at runtime
jags_model_dm <- "
model {
  # Priors
  lambda ~ dunif(0, 30)
  gamma  ~ dunif(0, 100)
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
  ks <- 0:K
  dat <- list(y = y, M = M, T = T, K = K, ks = ks)

  f <- function(par) {
    getAll(par, dat)
    lambda <- exp(log_lambda)
    gamma <- exp(log_gamma)
    omega <- plogis(logit_omega)
    p <- plogis(logit_p)

    # build transition matrix P(N_t+1=m | N_t=n)
    # rows = N_t=n, cols = N_t+1=m
    # P(N_t+1=m|N_t=n) = sum_s P(S=s|n,omega) * P(G=m-s|gamma)
    trans <- matrix(0, K + 1, K + 1)
    for (n in 0:K) {
      surv_probs <- dbinom(0:n, n, omega)
      recr_probs <- dpois(0:K, gamma)
      for (m in 0:K) {
        tp <- 0
        for (s in 0:min(n, m)) {
          tp <- tp + surv_probs[s + 1] * recr_probs[m - s + 1]
        }
        trans[n + 1, m + 1] <- tp
      }
    }

    nll <- 0
    for (i in 1:M) {
      alpha <- dpois(ks, lambda) * dbinom(y[i, 1], ks, p)
      for (t in 2:T) {
        alpha <- as.vector(t(trans) %*% alpha) * dbinom(y[i, t], ks, p)
      }
      nll <- nll - log(sum(alpha))
    }
    nll
  }

  par <- list(
    log_lambda  = log(mean(y[, 1]) + 0.1),
    log_gamma   = log(1.5),
    logit_omega = 0,
    logit_p     = 0
  )

  obj <- tryCatch(MakeADFun(f, par), error = function(e) NULL)
  if (is.null(obj)) {
    return(c(lambda = NA, gamma = NA, omega = NA, p = NA))
  }

  opt <- tryCatch(nlminb(obj$par, obj$fn, obj$gr), error = function(e) NULL)
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
                           seed = 333,
                           n.chains = 3, n.iter = 5000,
                           n.burnin = 2500, n.thin = 1) {
  set.seed(seed)
  datasets <- lapply(1:nsim, function(s) {
    sim_dm(M, T, lambda_true, gamma_true, omega_true, p_true)
  })

  # RTMB - all datasets timed together (matches original behaviour)
  time_rtmb <- system.time({
    res_rtmb <- lapply(datasets, function(y) {
      tryCatch(fit_rtmb(y),
        error = function(e) c(lambda = NA, gamma = NA, omega = NA, p = NA)
      )
    })
  })

  # unmarked - all datasets timed together
  time_unm <- system.time({
    res_unm <- lapply(datasets, function(y) {
      tryCatch(fit_unm(y),
        error = function(e) c(lambda = NA, gamma = NA, omega = NA, p = NA)
      )
    })
  })

  # JAGS - each dataset timed individually so per-sim times are available
  jags_raw <- lapply(datasets, function(y) {
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
  })

  res_rtmb <- as.data.frame(do.call(rbind, res_rtmb))
  res_unm <- as.data.frame(do.call(rbind, res_unm))
  res_jags <- as.data.frame(
    do.call(rbind, lapply(jags_raw, function(x) x$estimates))
  )
  names(res_rtmb) <- c("lambda", "gamma", "omega", "p")
  names(res_unm) <- c("lambda", "gamma", "omega", "p")
  names(res_jags) <- c("lambda", "gamma", "omega", "p")

  # Drop the same simulation indices from all three frameworks so estimates
  # always correspond to the same dataset - critical for fair comparison
  ok <- complete.cases(res_rtmb) & complete.cases(res_unm) & complete.cases(res_jags)
  res_rtmb <- res_rtmb[ok, ]
  res_unm <- res_unm[ok, ]
  res_jags <- res_jags[ok, ]

  jags_times <- sapply(jags_raw, function(x) x$time_jags)
  conv_rate <- mean(sapply(jags_raw, function(x) x$jags_converged), na.rm = TRUE)

  list(
    model = "dailmadsen",
    estimates_rtmb = res_rtmb,
    estimates_unm = res_unm,
    estimates_jags = res_jags,
    times = data.frame(
      rtmb = time_rtmb["elapsed"],
      unm  = time_unm["elapsed"],
      jags = sum(jags_times, na.rm = TRUE) # total wall time, matches rtmb/unm convention
    ),
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
