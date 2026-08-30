# Dynamic (multi-season) Occupancy Model (MacKenzie et al. 2003)
#
# Ecological process (latent occupancy):
#   Z_i1      ~ Bernoulli(psi)
#   Z_i,t+1   ~ Bernoulli(Z_it * omega + (1 - Z_it) * lambda),  t = 1, ..., T - 1
#
# Observation process:
#   y_itj | Z_it ~ Bernoulli(Z_it * p),  j = 1, ..., nrep
#
# where:
#   psi    = initial occupancy probability
#   omega  = persistence (survival) probability
#   lambda = colonization probability
#   p      = detection probability
#   Z_it   = true (latent) occupancy at site i, season t
#   y_itj  = detection/non-detection at site i, season t, occasion j
#

library(RTMB)
library(unmarked)
library(R2jags)

# run parameters
# set by `make NSIM=X`; set nsim manually here instead if running standalone
nsim <- as.integer(Sys.getenv("NSIM"))
# set by `make SEED=X`; set seed manually here instead if running standalone
seed <- as.integer(Sys.getenv("SEED"))
M <- 200
T <- 4
nrep <- 4
psi_true <- 0.6
omega_true <- 0.8
lambda_true <- 0.2
p_true <- 0.5
n.chains <- 4
n.iter <- 5000
n.burnin <- 2500
n.thin <- 1

# JAGS model as a string - written to a temp file at runtime
jags_model_dynocc <- "
model {
  # Priors
  psi    ~ dunif(0, 1)
  omega  ~ dunif(0, 1)
  lambda ~ dunif(0, 1)
  p      ~ dunif(0, 1)

  # Likelihood
  for (i in 1:M) {
    z[i, 1] ~ dbern(psi)
    for (t in 1:(T - 1)) {
      psi_t[i, t] <- z[i, t] * omega + (1 - z[i, t]) * lambda
      z[i, t + 1] ~ dbern(psi_t[i, t])
    }
    for (t in 1:T) {
      for (j in 1:nrep) {
        y[i, t, j] ~ dbern(z[i, t] * p)
      }
    }
  }
}
"

sim_data_dynocc <- function(M, T, nrep, psi, omega, lambda, p) {
  Z <- matrix(NA, M, T)
  Z[, 1] <- rbinom(M, 1, psi)
  for (t in 1:(T - 1)) {
    survived <- Z[, t] * omega
    psi_t <- survived + (1 - survived) * lambda
    Z[, t + 1] <- rbinom(M, 1, psi_t)
  }

  y <- array(NA, dim = c(M, T, nrep))
  for (t in 1:T) {
    p_detect <- p * Z[, t]
    for (i in 1:nrep) {
      y[, t, i] <- rbinom(M, 1, p_detect)
    }
  }
  list(y = y, M = M, T = T, nrep = nrep)
}

na_estimates_dynocc <- function() {
  c(
    psi_rtmb = NA, omega_rtmb = NA, lambda_rtmb = NA, p_rtmb = NA,
    psi_unm = NA, omega_unm = NA, lambda_unm = NA, p_unm = NA,
    psi_jags = NA, omega_jags = NA, lambda_jags = NA, p_jags = NA
  )
}

fit_all_dynocc <- function(seed, M, T, nrep, psi_true, omega_true,
                           lambda_true, p_true,
                           n.chains = 3, n.iter = 5000,
                           n.burnin = 2500, n.thin = 1) {
  gc() # get rid of residual JAGS misery
  set.seed(seed)
  dat_sim <- sim_data_dynocc(
    M, T, nrep, psi_true,
    omega_true, lambda_true, p_true
  )
  y <- dat_sim$y

  # ------------------------------------------------------------------
  # unmarked - colext() is the dynamic/multi-season occupancy model:
  # col = colonization (lambda), ext = extinction (1 - omega)
  # ------------------------------------------------------------------
  y_wide <- matrix(aperm(y, c(1, 3, 2)), nrow = M) # season-major: M x (nrep * T)
  umf <- unmarkedMultFrame(y = y_wide, numPrimary = T)
  time_unm <- system.time({
    fit_unm <- tryCatch(
      suppressWarnings(colext(~1, ~1, ~1, ~1, data = umf)),
      error = function(e) NULL
    )
  })

  if (is.null(fit_unm)) {
    return(list(
      estimates = na_estimates_dynocc(),
      time_rtmb = NA, time_unm = NA, time_jags = NA,
      jags_converged = NA, rtmb_converged = NA, unm_converged = NA
    ))
  }

  unm_coef <- coef(fit_unm)
  se_unm <- tryCatch(SE(fit_unm), error = function(e) NULL)
  unm_converged <- fit_unm@opt$convergence == 0 && !is.null(se_unm) && all(is.finite(se_unm))

  # ------------------------------------------------------------------
  # JAGS - same dataset, explicit latent Z_it sampled by MCMC
  # ------------------------------------------------------------------
  model_file <- file.path(tempdir(), "jags_model_dynocc.txt")
  writeLines(jags_model_dynocc, model_file)

  jags_data <- list(y = y, M = M, T = T, nrep = nrep)
  jags_inits <- function() {
    list(
      psi    = runif(1, 0.1, 0.9),
      omega  = runif(1, 0.1, 0.9),
      lambda = runif(1, 0.1, 0.9),
      p      = runif(1, 0.1, 0.9),
      z      = apply(y, MARGIN = 1:2, FUN = max) # initialise z at observed max
    )
  }

  jags.seed <- sample.int(1e6, 1)
  time_jags <- system.time({
    fit_jags <- tryCatch(
      suppressWarnings(
        jags.parallel(
          data = jags_data,
          inits = jags_inits,
          parameters.to.save = c("psi", "omega", "lambda", "p"),
          model.file = model_file,
          n.chains = n.chains,
          n.cluster = n.chains,
          n.iter = n.iter,
          n.burnin = n.burnin,
          n.thin = n.thin,
          jags.seed = jags.seed,
          envir = environment(),
          export_obj_names = c("y", "n.chains", "n.iter", "n.burnin", "n.thin", "jags.seed")
        )
      ),
      error = function(e) NULL
    )
  })
  if (is.null(fit_jags)) {
    jags_psi <- NA
    jags_omega <- NA
    jags_lambda <- NA
    jags_p <- NA
    jags_converged <- NA
  } else {
    sums <- fit_jags$BUGSoutput$summary
    jags_psi <- sums["psi", "mean"]
    jags_omega <- sums["omega", "mean"]
    jags_lambda <- sums["lambda", "mean"]
    jags_p <- sums["p", "mean"]
    # Rhat < 1.1 for all monitored params (excluding deviance)
    rhats <- sums[rownames(sums) != "deviance", "Rhat"]
    jags_converged <- all(rhats < 1.1, na.rm = TRUE)
  }

  # ------------------------------------------------------------------
  # RTMB
  # ------------------------------------------------------------------
  dat <- list(y = y, M = M, T = T)

  par <- list(
    logit_psi = 0,
    logit_omega = 0,
    logit_lambda = 0,
    logit_p = 0,
    Z = matrix(1, nrow = M, ncol = T)
  )

  f <- function(par) {
    getAll(dat, par, warn = FALSE)
    psi <- plogis(logit_psi)
    omega <- plogis(logit_omega)
    lambda <- plogis(logit_lambda)
    p <- plogis(logit_p)
    jnll <- -sum(dbinom(Z[, 1], 1, psi, log = TRUE))
    for (t in 1:(T - 1)) {
      survived <- Z[, t] * omega
      psi_t <- survived + (1 - survived) * lambda
      jnll <- jnll - sum(dbinom(Z[, t + 1], 1, psi_t, log = TRUE))
    }
    p_detect <- p * Z
    jnll <- jnll - sum(dbinom(y, 1, p_detect, log = TRUE))
    jnll
  }

  time_rtmb <- system.time({
    obj <- tryCatch(
      MakeADFun(
        f, par,
        random = "Z",
        map = list(Z = factor(ifelse(
          apply(y, MARGIN = 1:2, FUN = \(x) any(x == 1)),
          NA,
          seq_len(prod(dim(y)))
        ))),
        integrate = list(Z = TMB::SR(0:1, discrete = TRUE))
      ),
      error = function(e) NULL
    )
    if (!is.null(obj)) {
      opt <- tryCatch(
        nlminb(obj$par, obj$fn, obj$gr,
          control = list(eval.max = 1e4, iter.max = 1e4)
        ),
        error = function(e) NULL
      )
    }
  })

  if (is.null(obj) || is.null(opt) || opt$convergence != 0) {
    return(list(
      estimates = na_estimates_dynocc(),
      time_rtmb = NA, time_unm = NA, time_jags = NA,
      jags_converged = NA, rtmb_converged = NA, unm_converged = NA
    ))
  }

  sdr <- tryCatch(sdreport(obj), error = function(e) NULL)
  rtmb_converged <- !is.null(sdr) && sdr$pdHess

  list(
    estimates = c(
      psi_rtmb    = unname(plogis(opt$par["logit_psi"])),
      omega_rtmb  = unname(plogis(opt$par["logit_omega"])),
      lambda_rtmb = unname(plogis(opt$par["logit_lambda"])),
      p_rtmb      = unname(plogis(opt$par["logit_p"])),
      psi_unm     = unname(plogis(unm_coef["psi(Int)"])),
      omega_unm   = unname(1 - plogis(unm_coef["ext(Int)"])),
      lambda_unm  = unname(plogis(unm_coef["col(Int)"])),
      p_unm       = unname(plogis(unm_coef["p(Int)"])),
      psi_jags    = jags_psi,
      omega_jags  = jags_omega,
      lambda_jags = jags_lambda,
      p_jags      = jags_p
    ),
    time_rtmb = time_rtmb["elapsed"],
    time_unm = time_unm["elapsed"],
    time_jags = time_jags["elapsed"],
    jags_converged = jags_converged,
    rtmb_converged = rtmb_converged,
    unm_converged = unm_converged
  )
}

run_dynamic_occupancy <- function(nsim, M, T, nrep,
                                  psi_true, omega_true,
                                  lambda_true, p_true,
                                  seed,
                                  n.chains, n.iter,
                                  n.burnin, n.thin) {
  raw <- vector("list", nsim)
  for (s in 1:nsim) {
    raw[[s]] <- fit_all_dynocc(seed + s, M, T, nrep, psi_true,
      omega_true, lambda_true, p_true,
      n.chains = n.chains, n.iter = n.iter,
      n.burnin = n.burnin, n.thin = n.thin
    )
  }

  estimates_all <- do.call(rbind, lapply(raw, function(x) x$estimates))
  times_all <- data.frame(
    rtmb = sapply(raw, function(x) x$time_rtmb),
    unm  = sapply(raw, function(x) x$time_unm),
    jags = sapply(raw, function(x) x$time_jags)
  )

  jags_conv <- sapply(raw, function(x) x$jags_converged)
  rtmb_conv <- sapply(raw, function(x) x$rtmb_converged)
  unm_conv <- sapply(raw, function(x) x$unm_converged)

  # Keep only sims where all three frameworks succeeded and converged
  ok <- complete.cases(estimates_all) & complete.cases(times_all) &
    jags_conv %in% TRUE & rtmb_conv %in% TRUE & unm_conv %in% TRUE
  estimates <- as.data.frame(estimates_all[ok, , drop = FALSE])
  times <- times_all[ok, ]

  conv_rate <- mean(jags_conv, na.rm = TRUE)
  rtmb_conv_rate <- mean(rtmb_conv, na.rm = TRUE)
  unm_conv_rate <- mean(unm_conv, na.rm = TRUE)

  list(
    model = "dynamic_occupancy",
    estimates = estimates,
    times = times,
    convergence = data.frame(
      sim = seq_len(nsim),
      jags = jags_conv, rtmb = rtmb_conv, unm = unm_conv
    ),
    truth = c(
      psi = psi_true, omega = omega_true,
      lambda = lambda_true, p = p_true
    ),
    jags_conv_rate = conv_rate,
    rtmb_conv_rate = rtmb_conv_rate,
    unm_conv_rate = unm_conv_rate,
    nsim = nsim, M = M, T = T, nrep = nrep
  )
}

res <- run_dynamic_occupancy(
  nsim = nsim, M = M, T = T, nrep = nrep,
  psi_true = psi_true, omega_true = omega_true,
  lambda_true = lambda_true, p_true = p_true,
  seed = seed,
  n.chains = n.chains, n.iter = n.iter,
  n.burnin = n.burnin, n.thin = n.thin
)
saveRDS(res, "results/dynamic_occupancy.rds")
cat("Mean RTMB time:", round(mean(res$times$rtmb, na.rm = TRUE), 3), "s\n")
cat("Mean unmarked time:", round(mean(res$times$unm, na.rm = TRUE), 3), "s\n")
cat("Mean JAGS time:", round(mean(res$times$jags, na.rm = TRUE), 3), "s\n")
cat("JAGS convergence rate:", round(res$jags_conv_rate, 3), "\n")
cat("RTMB Hessian convergence rate:", round(res$rtmb_conv_rate, 3), "\n")
cat("unmarked Hessian convergence rate:", round(res$unm_conv_rate, 3), "\n")
