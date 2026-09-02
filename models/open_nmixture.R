# Dail-Madsen Model (Dail & Madsen 2011) - Open Population N-mixture
#
# Ecological process (initial abundance):
#   N_j1 ~ Poisson(lambda),  j = 1, ..., J
#
# Population dynamics (survival and recruitment):
#   S_jt | N_jt ~ Binomial(N_jt, omega),  t = 1, ..., T-1
#   G_jt ~ Poisson(gamma)
#   N_j,t+1 = S_jt + G_jt
#
# Observation process (detection):
#   y_jt | N_jt ~ Binomial(N_jt, p),  t = 1, ..., T
#
# where:
#   lambda = expected initial abundance at each site
#   omega  = apparent survival probability
#   gamma  = recruitment rate (births and immigration, constant across sites)
#   p      = detection probability
#   N_jt   = true (latent) abundance at site j, time t
#   S_jt   = survivors at site j from time t to t+1
#   G_jt   = recruits at site j from time t to t+1
#   y_jt   = count at site j, time t
#

library(RTMB)
library(unmarked)
library(R2jags)

# run parameters
# set by `make NSIM=X`; set nsim manually here instead if running standalone
nsim <- as.integer(Sys.getenv("NSIM"))
# set by `make SEED=X`; set seed manually here instead if running standalone
seed <- as.integer(Sys.getenv("SEED"))
J <- 100
T <- 5
lambda_true <- 4
gamma_true <- 1.5
omega_true <- 0.8
p_true <- 0.5
n.chains <- 4
n.iter <- 20000
n.burnin <- 10000
n.thin <- 1

# JAGS model
jags_model_dm <- "
model {
  # Priors
  lambda ~ dunif(0, 30)
  gamma  ~ dunif(0, 30)
  omega  ~ dunif(0, 1)
  p      ~ dunif(0, 1)

  # Likelihood
  for (j in 1:J) {
    # Initial abundance
    N[j, 1] ~ dpois(lambda)
    y[j, 1] ~ dbin(p, N[j, 1])

    for (t in 2:T) {
      # Survivors from t-1 to t
      S[j, t] ~ dbin(omega, N[j, t-1])
      # Recruits at t
      G[j, t] ~ dpois(gamma)
      # Total abundance at t
      N[j, t] <- S[j, t] + G[j, t]
      # Observation
      y[j, t] ~ dbin(p, N[j, t])
    }
  }
}
"

sim_dm <- function(J, T, lambda, gamma, omega, p) {
  N <- matrix(NA, J, T)
  N[, 1] <- rpois(J, lambda)
  for (t in 1:(T - 1)) {
    S <- rbinom(J, N[, t], omega)
    G <- rpois(J, gamma)
    N[, t + 1] <- S + G
  }
  matrix(rbinom(J * T, N, p), J, T)
}

fit_unm <- function(y) {
  K <- max(y) * 2
  umf <- unmarkedFramePCO(y = y, numPrimary = ncol(y))
  fit <- tryCatch(
    suppressWarnings(pcountOpen(~1, ~1, ~1, ~1, data = umf, K = K, se = TRUE)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(estimate = c(
      lambda = NA, gamma = NA,
      omega = NA, p = NA
    ), converged = NA))
  }

  se_unm <- tryCatch(SE(fit), error = function(e) NULL)
  converged <- fit@opt$convergence == 0 && !is.null(se_unm) && all(is.finite(se_unm))

  result <- c(
    exp(coef(fit, "lambda")),
    exp(coef(fit, "gamma")),
    plogis(coef(fit, "omega")),
    plogis(coef(fit, "det"))
  )
  names(result) <- c("lambda", "gamma", "omega", "p")
  list(estimate = result, converged = converged)
}

fit_jags_dm <- function(y, n.chains, n.iter, n.burnin, n.thin,
                        lambda_true) {
  J <- nrow(y)
  T <- ncol(y)

  model_file <- file.path(tempdir(), "jags_model_dm.txt")
  writeLines(jags_model_dm, model_file)

  # N[j,t] for t>=2 is a deterministic node (N <- S+G), cannot be initialised.
  # Only stochastic nodes need inits: N[j,1], S[j,t>=2], G[j,t>=2].
  N1_init <- pmax(apply(y, 1, max) + 2L, y[, 1] + 1L)
  S_init <- matrix(NA_integer_, J, T)
  G_init <- matrix(NA_integer_, J, T)
  N_prev <- N1_init
  for (t in 2:T) {
    S_init[, t] <- pmin(N_prev, round(N_prev * 0.8))
    G_init[, t] <- pmax(1L, y[, t] + 1L)
    N_prev <- S_init[, t] + G_init[, t]
  }

  jags_data <- list(y = y, J = J, T = T)
  jags_inits <- function() {
    list(
      lambda = runif(1, 1, lambda_true * 2),
      gamma  = runif(1, 0.5, 3),
      omega  = runif(1, 0.5, 0.99),
      p      = runif(1, 0.1, 0.9),
      N      = cbind(N1_init, matrix(NA_integer_, J, T - 1)),
      S      = S_init,
      G      = G_init
    )
  }

  # jags.parallel() defaults jags.seed to a hardcoded 123, which - unlike
  # sequential jags() inheriting R's naturally-advancing RNG - makes every
  # chain's starting values IDENTICAL across every replicate and every run,
  # since seeds <- jags.seed + seq_len(n.chains) never changes. Draw a fresh
  # one from R's already-seeded stream so inits actually vary..
  jags.seed <- sample.int(1e6, 1)

  fit <- tryCatch(
    suppressWarnings(
      jags.parallel(
        data = jags_data,
        inits = jags_inits,
        parameters.to.save = c("lambda", "gamma", "omega", "p"),
        model.file = model_file,
        n.chains = n.chains,
        n.cluster = n.chains,
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.thin = n.thin,
        jags.seed = jags.seed,
        envir = environment(),
        export_obj_names = c(
          "N1_init", "S_init", "G_init", "lambda_true",
          "n.chains", "n.iter", "n.burnin", "n.thin", "jags.seed"
        )
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

fit_rtmb <- function(y) {
  J <- nrow(y)
  T <- ncol(y)
  K <- max(y) * 2
  dat <- list(y = y, J = J, T = T)

  f <- function(par) {
    getAll(par, dat)
    S <- SN[, seq_len(T - 1)]
    N <- SN[, T - 1 + seq_len(T)]
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
    log_lambda = log(mean(y[, 1]) + 0.1),
    log_gamma = log(1.5),
    logit_omega = 0,
    logit_p = 0,
    SN = matrix(K, nrow = J, ncol = 2 * T - 1)
  )

  obj <- tryCatch(
    MakeADFun(f, par,
      random = c("SN"),
      integrate = list(
        SN = TMB::SR(0:K, discrete = TRUE)
      ),
      silent = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(obj)) {
    return(list(estimate = c(
      lambda = NA, gamma = NA,
      omega = NA, p = NA
    ), converged = NA))
  }

  opt <- tryCatch(
    nlminb(obj$par, obj$fn, obj$gr,
      control = list(iter.max = 1e5, eval.max = 1e5)
    ),
    error = function(e) NULL
  )
  if (is.null(opt) || opt$convergence != 0) {
    return(list(estimate = c(
      lambda = NA,
      gamma = NA, omega = NA, p = NA
    ), converged = NA))
  }

  sdr <- tryCatch(sdreport(obj), error = function(e) NULL)
  rtmb_converged <- !is.null(sdr) && sdr$pdHess

  result <- c(
    exp(opt$par["log_lambda"]),
    exp(opt$par["log_gamma"]),
    plogis(opt$par["logit_omega"]),
    plogis(opt$par["logit_p"])
  )
  names(result) <- c("lambda", "gamma", "omega", "p")
  list(estimate = result, converged = rtmb_converged)
}

run_dailmadsen <- function(nsim, J, T,
                           lambda_true, gamma_true,
                           omega_true, p_true,
                           seed,
                           n.chains, n.iter,
                           n.burnin, n.thin) {
  datasets <- lapply(1:nsim, function(s) {
    set.seed(seed + s)
    sim_dm(J, T, lambda_true, gamma_true, omega_true, p_true)
  })

  unm_raw <- vector("list", nsim)
  for (i in seq_along(datasets)) {
    y <- datasets[[i]]
    t <- system.time({
      out <- tryCatch(fit_unm(y),
        error = function(e) {
          list(estimate = c(
            lambda = NA,
            gamma = NA, omega = NA, p = NA
          ), converged = NA)
        }
      )
    })
    unm_raw[[i]] <- list(
      estimate = out$estimate,
      unm_converged = out$converged, time = t["elapsed"]
    )
  }

  jags_raw <- vector("list", nsim)
  for (i in seq_along(datasets)) {
    y <- datasets[[i]]
    t <- system.time({
      out <- fit_jags_dm(y,
        n.chains = n.chains,
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.thin = n.thin,
        lambda_true = lambda_true
      )
    })
    jags_raw[[i]] <- list(
      estimates = out$estimates,
      jags_converged = out$jags_converged,
      time_jags = t["elapsed"]
    )
  }

  rtmb_raw <- vector("list", nsim)
  for (i in seq_along(datasets)) {
    y <- datasets[[i]]
    t <- system.time({
      out <- tryCatch(fit_rtmb(y),
        error = function(e) {
          list(estimate = c(
            lambda = NA,
            gamma = NA, omega = NA, p = NA
          ), converged = NA)
        }
      )
    })
    rtmb_raw[[i]] <- list(
      estimate = out$estimate,
      rtmb_converged = out$converged, time = t["elapsed"]
    )
  }

  res_unm <- as.data.frame(do.call(rbind, lapply(unm_raw, function(x) x$estimate)))
  res_jags <- as.data.frame(
    do.call(rbind, lapply(jags_raw, function(x) x$estimates))
  )
  res_rtmb <- as.data.frame(do.call(rbind, lapply(rtmb_raw, function(x) x$estimate)))
  names(res_unm) <- c("lambda", "gamma", "omega", "p")
  names(res_jags) <- c("lambda", "gamma", "omega", "p")
  names(res_rtmb) <- c("lambda", "gamma", "omega", "p")

  unm_times <- sapply(unm_raw, function(x) x$time)
  jags_times <- sapply(jags_raw, function(x) x$time_jags)
  rtmb_times <- sapply(rtmb_raw, function(x) x$time)

  jags_conv <- sapply(jags_raw, function(x) x$jags_converged)
  rtmb_conv <- sapply(rtmb_raw, function(x) x$rtmb_converged)
  unm_conv <- sapply(unm_raw, function(x) x$unm_converged)

  # Drop the same simulation indices from all three frameworks so estimates
  # always correspond to the same dataset for fair comparison. JAGS
  # convergence is not required: this model is hard enough for JAGS that
  # requiring Rhat < 1.1 leaves too few replicates
  ok <- complete.cases(res_rtmb) & complete.cases(res_unm) & complete.cases(res_jags) &
    rtmb_conv %in% TRUE & unm_conv %in% TRUE
  res_rtmb <- res_rtmb[ok, ]
  res_unm <- res_unm[ok, ]
  res_jags <- res_jags[ok, ]

  conv_rate <- mean(jags_conv, na.rm = TRUE)
  rtmb_conv_rate <- mean(rtmb_conv, na.rm = TRUE)
  unm_conv_rate <- mean(unm_conv, na.rm = TRUE)

  list(
    model = "dailmadsen",
    estimates_unm = res_unm,
    estimates_jags = res_jags,
    estimates_rtmb = res_rtmb,
    convergence = data.frame(
      sim = seq_len(nsim),
      jags = jags_conv, rtmb = rtmb_conv, unm = unm_conv
    ),
    times = data.frame(
      rtmb = sum(rtmb_times, na.rm = TRUE),
      unm  = sum(unm_times, na.rm = TRUE),
      jags = sum(jags_times, na.rm = TRUE)
    ),
    unm_times_per_sim = unm_times[ok],
    jags_times_per_sim = jags_times[ok],
    rtmb_times_per_sim = rtmb_times[ok],
    jags_conv_rate = conv_rate,
    rtmb_conv_rate = rtmb_conv_rate,
    unm_conv_rate = unm_conv_rate,
    truth = c(
      lambda = lambda_true, gamma = gamma_true,
      omega = omega_true, p = p_true
    ),
    nsim = nsim, J = J, T = T
  )
}

res <- run_dailmadsen(
  nsim = nsim, J = J, T = T,
  lambda_true = lambda_true, gamma_true = gamma_true,
  omega_true = omega_true, p_true = p_true,
  seed = seed,
  n.chains = n.chains, n.iter = n.iter,
  n.burnin = n.burnin, n.thin = n.thin
)
saveRDS(res, "results/open_nmixture.rds")
cat("Total RTMB time:", res$times$rtmb, "s\n")
cat("Total unmarked time:", res$times$unm, "s\n")
cat("Total JAGS time:", res$times$jags, "s\n")
cat("JAGS convergence rate:", round(res$jags_conv_rate, 3), "\n")
cat("RTMB Hessian convergence rate:", round(res$rtmb_conv_rate, 3), "\n")
cat("unmarked Hessian convergence rate:", round(res$unm_conv_rate, 3), "\n")
