# N-mixture Model (Royle 2004) - Monte Carlo Simulation
#
# Ecological process (latent abundance):
#   N_i ~ Poisson(lambda),  i = 1, ..., R
#
# Observation process (detection):
#   y_ij | N_i ~ Binomial(N_i, p),  j = 1, ..., T
#
# where:
#   lambda = expected abundance at each site
#   p      = probability of detecting an individual
#   N_i    = true (latent) abundance at site i
#   y_ij   = count at site i, occasion j
#
# The marginal likelihood integrates out N_i via sequential reduction (SR):
#   L(lambda, p) = prod_i sum_{k=0}^{K} Poisson(k|lambda) * prod_j Binomial(y_ij|k,p)
#
# JAGS fits the same hierarchical model via MCMC, sampling N_i explicitly.

library(RTMB)
library(unmarked)
library(R2jags)

# run parameters
# set by `make NSIM=X`; set nsim manually here instead if running standalone
nsim <- as.integer(Sys.getenv("NSIM"))
# set by `make SEED=X`; set seed manually here instead if running standalone
seed <- as.integer(Sys.getenv("SEED"))
R <- 100
T <- 5
lambda_true <- 32
p_true <- 0.25
n.chains <- 4
n.iter <- 20000
n.burnin <- 10000
n.thin <- 1

# JAGS model as a string - written to a temp file at runtime
jags_model_nmix <- "
model {
  # Priors
  lambda ~ dunif(0, 80)
  p      ~ dunif(0, 1)

  # Likelihood
  for (i in 1:R) {
    N[i] ~ dpois(lambda)
    for (j in 1:T) {
      y[i, j] ~ dbin(p, N[i])
    }
  }
}
"

sim_data_nmix <- function(R, T, lambda, p) {
  N_true <- rpois(R, lambda)
  y <- matrix(rbinom(R * T, rep(N_true, T), p), R, T)
  list(y = y, R = R, T = T)
}

fit_all_nmix <- function(seed, R, T, lambda_true, p_true,
                         n.chains = 3, n.iter = 5000, n.burnin = 2500, n.thin = 1) {
  set.seed(seed)
  dat <- sim_data_nmix(R, T, lambda_true, p_true)
  y <- dat$y
  K <- max(y) * 3 # buffer around K

  # ------------------------------------------------------------------
  # unmarked
  # ------------------------------------------------------------------
  umf <- unmarkedFramePCount(y = y)
  time_unm <- system.time({
    fit_unm <- tryCatch(
      pcount(~1 ~ 1, data = umf, K = K),
      error = function(e) NULL
    )
  })

  if (is.null(fit_unm)) {
    return(list(
      estimates = c(
        lambda_rtmb = NA, p_rtmb = NA,
        lambda_unm = NA, p_unm = NA,
        lambda_jags = NA, p_jags = NA
      ),
      time_rtmb = NA, time_unm = NA, time_jags = NA,
      jags_converged = NA
    ))
  }

  # ------------------------------------------------------------------
  # JAGS - same dataset, N_i sampled explicitly by MCMC
  # ------------------------------------------------------------------
  model_file <- file.path(tempdir(), "jags_model_nmix.txt")
  writeLines(jags_model_nmix, model_file)

  jags_data <- list(y = y, R = R, T = T)
  jags_inits <- function() {
    list(
      lambda = runif(1, 1, lambda_true * 2),
      p      = runif(1, 0.1, 0.9),
      N      = apply(y, 1, max) + 1L # initialise above observed max
    )
  }

  # This is a timing comparison against RTMB/unmarked, so JAGS chains run
  # in parallel (one per core) rather than sequentially.
  jags.seed <- sample.int(1e6, 1)
  time_jags <- system.time({
    fit_jags <- tryCatch(
      suppressWarnings(
        jags.parallel(
          data = jags_data,
          inits = jags_inits,
          parameters.to.save = c("lambda", "p"),
          model.file = model_file,
          n.chains = n.chains,
          n.cluster = n.chains,
          n.iter = n.iter,
          n.burnin = n.burnin,
          n.thin = n.thin,
          jags.seed = jags.seed,
          envir = environment(),
          export_obj_names = c("y", "lambda_true", "n.chains", "n.iter", "n.burnin", "n.thin", "jags.seed")
        )
      ),
      error = function(e) NULL
    )
  })

  if (is.null(fit_jags)) {
    jags_lambda <- NA
    jags_p <- NA
    jags_converged <- NA
  } else {
    sums <- fit_jags$BUGSoutput$summary
    jags_lambda <- sums["lambda", "mean"]
    jags_p <- sums["p", "mean"]
    rhats <- sums[rownames(sums) != "deviance", "Rhat"]
    jags_converged <- all(rhats < 1.1, na.rm = TRUE)
  }

  # ------------------------------------------------------------------
  # RTMB
  # ------------------------------------------------------------------
  f <- function(par) {
    getAll(par, dat)
    nll <- 0
    for (i in 1:R) {
      nll <- nll - dpois(N[i], exp(log_lambda), log = TRUE)
      for (j in 1:T) {
        nll <- nll - dbinom(y[i, j], N[i], plogis(logit_p), log = TRUE)
      }
    }
    nll
  }
  par <- list(log_lambda = log(mean(y) + 0.1), logit_p = 0, N = apply(y, 1, max))

  time_rtmb <- system.time({
    obj <- tryCatch(
      MakeADFun(f, par,
        random = "N",
        integrate = list(N = TMB::SR(0:K, discrete = TRUE))
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
        lambda_rtmb = NA, p_rtmb = NA,
        lambda_unm = NA, p_unm = NA,
        lambda_jags = NA, p_jags = NA
      ),
      time_rtmb = NA, time_unm = NA, time_jags = NA,
      jags_converged = NA
    ))
  }

  list(
    estimates = c(
      lambda_rtmb = unname(exp(opt$par["log_lambda"])),
      p_rtmb      = unname(plogis(opt$par["logit_p"])),
      lambda_unm  = unname(exp(coef(fit_unm, "state"))),
      p_unm       = unname(plogis(coef(fit_unm, "det"))),
      lambda_jags = jags_lambda,
      p_jags      = jags_p
    ),
    time_rtmb = time_rtmb["elapsed"],
    time_unm = time_unm["elapsed"],
    time_jags = time_jags["elapsed"],
    jags_converged = jags_converged
  )
}

run_nmixture <- function(nsim, R, T,
                         lambda_true, p_true,
                         seed,
                         n.chains, n.iter,
                         n.burnin, n.thin) {
  raw <- vector("list", nsim)
  for (s in 1:nsim) {
    raw[[s]] <- fit_all_nmix(seed + s, R, T, lambda_true, p_true,
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

  # Keep only sims where all three frameworks succeeded
  ok <- complete.cases(estimates_all) & complete.cases(times_all)
  estimates <- as.data.frame(estimates_all[ok, , drop = FALSE])
  times <- times_all[ok, ]

  conv_rate <- mean(sapply(raw, function(x) x$jags_converged), na.rm = TRUE)

  list(
    model = "nmixture",
    estimates = estimates,
    times = times,
    truth = c(lambda = lambda_true, p = p_true),
    jags_conv_rate = conv_rate,
    nsim = nsim, R = R, T = T
  )
}

res <- run_nmixture(
  nsim = nsim, R = R, T = T,
  lambda_true = lambda_true, p_true = p_true,
  seed = seed,
  n.chains = n.chains, n.iter = n.iter,
  n.burnin = n.burnin, n.thin = n.thin
)
saveRDS(res, "results/nmixture.rds")
cat("Mean RTMB time:", round(mean(res$times$rtmb, na.rm = TRUE), 3), "s\n")
cat("Mean unmarked time:", round(mean(res$times$unm, na.rm = TRUE), 3), "s\n")
cat("Mean JAGS time:", round(mean(res$times$jags, na.rm = TRUE), 3), "s\n")
cat("JAGS convergence rate:", round(res$jags_conv_rate, 3), "\n")
