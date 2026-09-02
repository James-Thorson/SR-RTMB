# Occupancy Model (MacKenzie et al. 2002)
#
# Ecological process (latent occupancy):
#   z_j ~ Bernoulli(psi),  j = 1, ..., J
#
# Observation process:
#   y_ji | z_j ~ Bernoulli(z_j * p),  i = 1, ..., I
#
# where:
#   psi  = occupancy probability
#   p    = detection probability
#   z_j  = true (latent) occupancy at site j
#   y_ji = detection/non-detection at site j, sample i
#

library(RTMB)
library(unmarked)
library(R2jags)

# run parameters
# set by `make NSIM=X`; set nsim manually here instead if running standalone
nsim <- as.integer(Sys.getenv("NSIM"))
# set by `make SEED=X`; set seed manually here instead if running standalone
seed <- as.integer(Sys.getenv("SEED"))
J <- 200
I <- 5
psi_true <- 0.2
p_true <- 0.5
n.chains <- 4
n.iter <- 5000
n.burnin <- 2500
n.thin <- 1

# JAGS model as a string - written to a temp file at runtime
jags_model_occ <- "
model {
  # Priors
  psi ~ dunif(0, 1)
  p   ~ dunif(0, 1)

  # Likelihood
  for (j in 1:J) {
    z[j] ~ dbern(psi)
    for (i in 1:I) {
      y[j, i] ~ dbern(z[j] * p)
    }
  }

  # Derived
  # (none needed - psi and p are monitored directly)
}
"

sim_data_occ <- function(J, I, psi, p) {
  z_true <- rbinom(J, 1, psi)
  y <- matrix(rbinom(J * I, 1, z_true * p), J, I)
  list(y = y, J = J, I = I)
}

fit_all_occ <- function(seed, J, I, psi_true, p_true,
                        n.chains = 3, n.iter = 5000, n.burnin = 2500, n.thin = 1) {
  gc() # get rid of residual JAGS misery
  set.seed(seed)
  dat <- sim_data_occ(J, I, psi_true, p_true)
  y <- dat$y

  # ------------------------------------------------------------------
  # unmarked
  # ------------------------------------------------------------------
  umf <- unmarkedFrameOccu(y = y)
  time_unm <- system.time({
    fit_unm <- tryCatch(
      suppressWarnings(occu(~1 ~ 1, data = umf)),
      error = function(e) NULL
    )
  })

  if (is.null(fit_unm)) {
    return(list(
      estimates = c(
        psi_rtmb = NA, p_rtmb = NA,
        psi_unm = NA, p_unm = NA,
        psi_jags = NA, p_jags = NA
      ),
      time_rtmb = NA, time_unm = NA, time_jags = NA,
      jags_converged = NA, rtmb_converged = NA, unm_converged = NA
    ))
  }

  se_unm <- tryCatch(SE(fit_unm), error = function(e) NULL)
  unm_converged <- fit_unm@opt$convergence == 0 && !is.null(se_unm) && all(is.finite(se_unm))

  # ------------------------------------------------------------------
  # JAGS - same dataset, explicit latent z_j sampled by MCMC
  # ------------------------------------------------------------------
  model_file <- file.path(tempdir(), "jags_model_occ.txt")
  writeLines(jags_model_occ, model_file)

  jags_data <- list(y = y, J = J, I = I)
  jags_inits <- function() {
    list(
      psi = runif(1, 0.1, 0.9),
      p   = runif(1, 0.1, 0.9),
      z   = apply(y, 1, max) # initialise z at observed max
    )
  }

  jags.seed <- sample.int(1e6, 1)
  time_jags <- system.time({
    fit_jags <- tryCatch(
      suppressWarnings(
        jags.parallel(
          data = jags_data,
          inits = jags_inits,
          parameters.to.save = c("psi", "p"),
          model.file = model_file,
          n.chains = n.chains,
          n.cluster = n.chains,
          n.iter = n.iter,
          n.burnin = n.burnin,
          n.thin = n.thin,
          jags.seed = jags.seed,
          envir = environment(),
          export_obj_names = c(
            "y", "n.chains", "n.iter",
            "n.burnin", "n.thin", "jags.seed"
          )
        )
      ),
      error = function(e) NULL
    )
  })
  if (is.null(fit_jags)) {
    jags_psi <- NA
    jags_p <- NA
    jags_converged <- NA
  } else {
    sums <- fit_jags$BUGSoutput$summary
    jags_psi <- sums["psi", "mean"]
    jags_p <- sums["p", "mean"]
    # Rhat < 1.1 for all monitored params (excluding deviance)
    rhats <- sums[rownames(sums) != "deviance", "Rhat"]
    jags_converged <- all(rhats < 1.1, na.rm = TRUE)
  }

  # ------------------------------------------------------------------
  # RTMB
  # ------------------------------------------------------------------
  f <- function(par) {
    getAll(par, dat)
    psi <- plogis(logit_psi)
    p <- plogis(logit_p)
    nll <- 0
    for (j in 1:J) {
      nll <- nll - dbinom(z[j], size = 1, prob = psi, log = TRUE)
      for (i in 1:I) {
        nll <- nll - dbinom(y[j, i], size = 1, prob = z[j] * p, log = TRUE)
      }
    }
    nll
  }
  par <- list(logit_psi = 0, logit_p = 0, z = rep(1, J))

  time_rtmb <- system.time({
    obj <- tryCatch(
      MakeADFun(f, par,
        random = "z",
        integrate = list(z = TMB::SR(c(0, 1), discrete = TRUE))
      ),
      error = function(e) NULL
    )
    if (!is.null(obj)) {
      opt <- tryCatch(nlminb(obj$par, obj$fn, obj$gr), error = function(e) NULL)
    }
  })

  if (is.null(obj) || is.null(opt) || opt$convergence != 0) {
    return(list(
      estimates = c(
        psi_rtmb = NA, p_rtmb = NA,
        psi_unm = NA, p_unm = NA,
        psi_jags = NA, p_jags = NA
      ),
      time_rtmb = NA, time_unm = NA, time_jags = NA,
      jags_converged = NA, rtmb_converged = NA, unm_converged = NA
    ))
  }

  sdr <- tryCatch(sdreport(obj), error = function(e) NULL)
  rtmb_converged <- !is.null(sdr) && sdr$pdHess

  list(
    estimates = c(
      psi_rtmb = unname(plogis(opt$par["logit_psi"])),
      p_rtmb   = unname(plogis(opt$par["logit_p"])),
      psi_unm  = unname(plogis(coef(fit_unm)[1])),
      p_unm    = unname(plogis(coef(fit_unm)[2])),
      psi_jags = jags_psi,
      p_jags   = jags_p
    ),
    time_rtmb = time_rtmb["elapsed"],
    time_unm = time_unm["elapsed"],
    time_jags = time_jags["elapsed"],
    jags_converged = jags_converged,
    rtmb_converged = rtmb_converged,
    unm_converged = unm_converged
  )
}

run_occupancy <- function(nsim, J, I,
                          psi_true, p_true,
                          seed,
                          n.chains, n.iter,
                          n.burnin, n.thin) {
  raw <- vector("list", nsim)
  for (s in 1:nsim) {
    raw[[s]] <- fit_all_occ(seed + s, J, I, psi_true, p_true,
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
    model = "occupancy",
    estimates = estimates,
    times = times,
    convergence = data.frame(
      sim = seq_len(nsim), jags = jags_conv, rtmb = rtmb_conv, unm = unm_conv
    ),
    truth = c(psi = psi_true, p = p_true),
    jags_conv_rate = conv_rate,
    rtmb_conv_rate = rtmb_conv_rate,
    unm_conv_rate = unm_conv_rate,
    nsim = nsim, J = J, I = I
  )
}

res <- run_occupancy(
  nsim = nsim, J = J, I = I,
  psi_true = psi_true, p_true = p_true,
  seed = seed,
  n.chains = n.chains, n.iter = n.iter,
  n.burnin = n.burnin, n.thin = n.thin
)
saveRDS(res, "results/occupancy.rds")
cat("Mean RTMB time:", round(mean(res$times$rtmb, na.rm = TRUE), 3), "s\n")
cat("Mean unmarked time:", round(mean(res$times$unm, na.rm = TRUE), 3), "s\n")
cat("Mean JAGS time:", round(mean(res$times$jags, na.rm = TRUE), 3), "s\n")
cat("JAGS convergence rate:", round(res$jags_conv_rate, 3), "\n")
cat("RTMB Hessian convergence rate:", round(res$rtmb_conv_rate, 3), "\n")
cat("unmarked Hessian convergence rate:", round(res$unm_conv_rate, 3), "\n")
