# Occupancy Model (MacKenzie et al. 2002) - Monte Carlo Simulation
#
# Ecological process (latent occupancy):
#   z_i ~ Bernoulli(psi),  i = 1, ..., R
#
# Observation process:
#   y_ij | z_i ~ Bernoulli(z_i * p),  j = 1, ..., T
#
# where:
#   psi  = occupancy probability
#   p    = detection probability
#   z_i  = true (latent) occupancy at site i
#   y_ij = detection/non-detection at site i, occasion j
#
# The marginal likelihood automatically marginalizes out z_i via SR:
#   L(psi, p) = prod_i sum_{k in {0,1}} Bernoulli(k|psi) * prod_j Bernoulli(y_ij|k*p)
#
# JAGS fits the same hierarchical model via MCMC, sampling z_i explicitly.

library(RTMB)
library(unmarked)
library(R2jags)

# JAGS model as a string - written to a temp file at runtime
jags_model_occ <- "
model {
  # Priors
  psi ~ dunif(0, 1)
  p   ~ dunif(0, 1)

  # Likelihood
  for (i in 1:R) {
    z[i] ~ dbern(psi)
    for (j in 1:T) {
      y[i, j] ~ dbern(z[i] * p)
    }
  }

  # Derived
  # (none needed - psi and p are monitored directly)
}
"

sim_data_occ <- function(R, T, psi, p) {
  z_true <- rbinom(R, 1, psi)
  y <- matrix(rbinom(R * T, 1, z_true * p), R, T)
  list(y = y, R = R, T = T)
}

fit_all_occ <- function(s, R, T, psi_true, p_true,
                     n.chains = 3, n.iter = 5000, n.burnin = 2500, n.thin = 1) {
  dat <- sim_data_occ(R, T, psi_true, p_true)
  y   <- dat$y

  # ------------------------------------------------------------------
  # RTMB
  # ------------------------------------------------------------------
  f <- function(par) {
    getAll(par, dat)
    psi <- plogis(logit_psi)
    p   <- plogis(logit_p)
    nll <- 0
    for (i in 1:R) {
      nll <- nll - dbinom(z[i], size = 1, prob = psi, log = TRUE)
      for (j in 1:T) {
        nll <- nll - dbinom(y[i, j], size = 1, prob = z[i] * p, log = TRUE)
      }
    }
    nll
  }
  par <- list(logit_psi = 0, logit_p = 0, z = rep(1, R))

  time_rtmb <- system.time({
    obj <- tryCatch(
      MakeADFun(f, par, random = "z",
                integrate = list(z = TMB::SR(c(0, 1), discrete = TRUE))),
      error = function(e) NULL
    )
    if (!is.null(obj))
      opt <- tryCatch(nlminb(obj$par, obj$fn, obj$gr), error = function(e) NULL)
  })

  if (is.null(obj) || is.null(opt) || opt$convergence != 0) {
    return(list(
      estimates = c(psi_rtmb = NA, p_rtmb = NA,
                    psi_unm  = NA, p_unm  = NA,
                    psi_jags = NA, p_jags = NA),
      time_rtmb = NA, time_unm = NA, time_jags = NA,
      jags_converged = NA
    ))
  }

  # ------------------------------------------------------------------
  # unmarked
  # ------------------------------------------------------------------
  umf      <- unmarkedFrameOccu(y = y)
  time_unm <- system.time({
    fit_unm <- tryCatch(
      suppressWarnings(occu(~1 ~ 1, data = umf)),
      error = function(e) NULL
    )
  })

  if (is.null(fit_unm)) {
    return(list(
      estimates = c(psi_rtmb = NA, p_rtmb = NA,
                    psi_unm  = NA, p_unm  = NA,
                    psi_jags = NA, p_jags = NA),
      time_rtmb = NA, time_unm = NA, time_jags = NA,
      jags_converged = NA
    ))
  }

  # ------------------------------------------------------------------
  # JAGS - same dataset, explicit latent z_i sampled by MCMC
  # ------------------------------------------------------------------
  model_file <- file.path(tempdir(), "jags_model_occ.txt")
  writeLines(jags_model_occ, model_file)

  jags_data  <- list(y = y, R = R, T = T)
  jags_inits <- function() list(
    psi = runif(1, 0.1, 0.9),
    p   = runif(1, 0.1, 0.9),
    z   = apply(y, 1, max)          # initialise z at observed max
  )

  time_jags <- system.time({
    fit_jags <- tryCatch(
      suppressWarnings(
        jags(data               = jags_data,
                      inits              = jags_inits,
                      parameters.to.save = c("psi", "p"),
                      model.file         = model_file,
                      n.chains           = n.chains,
                      n.iter             = n.iter,
                      n.burnin           = n.burnin,
                      n.thin             = n.thin,
                      progress.bar       = "none")
      ),
      error = function(e) NULL
    )
  })

  if (is.null(fit_jags)) {
    jags_psi       <- NA
    jags_p         <- NA
    jags_converged <- NA
  } else {
    sums           <- fit_jags$BUGSoutput$summary
    jags_psi       <- sums["psi", "mean"]
    jags_p         <- sums["p",   "mean"]
    # Rhat < 1.1 for all monitored params (excluding deviance)
    rhats          <- sums[rownames(sums) != "deviance", "Rhat"]
    jags_converged <- all(rhats < 1.1, na.rm = TRUE)
  }

  list(
    estimates = c(
      psi_rtmb = unname(plogis(opt$par["logit_psi"])),
      p_rtmb   = unname(plogis(opt$par["logit_p"])),
      psi_unm  = unname(plogis(coef(fit_unm)[1])),
      p_unm    = unname(plogis(coef(fit_unm)[2])),
      psi_jags = jags_psi,
      p_jags   = jags_p
    ),
    time_rtmb      = time_rtmb["elapsed"],
    time_unm       = time_unm["elapsed"],
    time_jags      = time_jags["elapsed"],
    jags_converged = jags_converged
  )
}

run_occupancy <- function(nsim = 1, R = 200, T = 5,
                          psi_true = 0.2, p_true = 0.5,
                          seed = 1123,
                          n.chains = 3, n.iter = 5000,
                          n.burnin = 2500, n.thin = 1) {
  set.seed(seed)
  raw <- lapply(1:nsim, function(s)
    fit_all_occ(s, R, T, psi_true, p_true,
             n.chains = n.chains, n.iter = n.iter,
             n.burnin = n.burnin, n.thin = n.thin))

  estimates <- as.data.frame(na.omit(
    do.call(rbind, lapply(raw, function(x) x$estimates))
  ))

  times <- na.omit(data.frame(
    rtmb = sapply(raw, function(x) x$time_rtmb),
    unm  = sapply(raw, function(x) x$time_unm),
    jags = sapply(raw, function(x) x$time_jags)
  ))

  conv_rate <- mean(sapply(raw, function(x) x$jags_converged), na.rm = TRUE)

  list(
    model         = "occupancy",
    estimates     = estimates,
    times         = times,
    truth         = c(psi = psi_true, p = p_true),
    jags_conv_rate = conv_rate,
    nsim = nsim, R = R, T = T
  )
}

# Run standalone if called directly
if (sys.nframe() == 0) {
  res <- run_occupancy()
  dir.create("results", showWarnings = FALSE)
  saveRDS(res, "results/occupancy.rds")
  cat("Mean RTMB time:",  round(mean(res$times$rtmb,  na.rm = TRUE), 3), "s\n")
  cat("Mean unmarked time:", round(mean(res$times$unm, na.rm = TRUE), 3), "s\n")
  cat("Mean JAGS time:",  round(mean(res$times$jags,  na.rm = TRUE), 3), "s\n")
  cat("JAGS convergence rate:", round(res$jags_conv_rate, 3), "\n")
}
