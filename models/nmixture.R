# N-mixture Model (Royle 2004) - Monte Carlo Simulation
#
# Ecological process (latent abundance):
#   N_j ~ Poisson(lambda),  j = 1, ..., J
#
# Observation process (detection):
#   y_ji | N_j ~ Binomial(N_j, p),  i = 1, ..., I
#
# where:
#   lambda = expected abundance at each site
#   p      = probability of detecting an individual
#   N_j    = true (latent) abundance at site j
#   y_ji   = count at site j, sample i

library(RTMB)
library(unmarked)
library(R2jags)

# run parameters
# set by `make NSIM=X`; set nsim manually here instead if running standalone
nsim <- as.integer(Sys.getenv("NSIM"))
# set by `make SEED=X`; set seed manually here instead if running standalone
seed <- as.integer(Sys.getenv("SEED"))
J <- 100
I <- 5
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
  for (j in 1:J) {
    N[j] ~ dpois(lambda)
    for (i in 1:I) {
      y[j, i] ~ dbin(p, N[j])
    }
  }
}
"

sim_data_nmix <- function(J, I, lambda, p) {
  N_true <- rpois(J, lambda)
  y <- matrix(rbinom(J * I, rep(N_true, I), p), J, I)
  list(y = y, J = J, I = I)
}

fit_all_nmix <- function(seed, J, I, lambda_true, p_true,
                         n.chains = 3, n.iter = 5000, n.burnin = 2500, n.thin = 1) {
  set.seed(seed)
  dat <- sim_data_nmix(J, I, lambda_true, p_true)
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
      jags_converged = NA, rtmb_converged = NA, unm_converged = NA
    ))
  }

  se_unm <- tryCatch(SE(fit_unm), error = function(e) NULL)
  unm_converged <- fit_unm@opt$convergence == 0 && !is.null(se_unm) && all(is.finite(se_unm))

  # ------------------------------------------------------------------
  # JAGS - same dataset, N_j sampled explicitly by MCMC
  # ------------------------------------------------------------------
  model_file <- file.path(tempdir(), "jags_model_nmix.txt")
  writeLines(jags_model_nmix, model_file)

  jags_data <- list(y = y, J = J, I = I)
  jags_inits <- function() {
    list(
      lambda = runif(1, 1, lambda_true * 2),
      p      = runif(1, 0.1, 0.9),
      N      = apply(y, 1, max) + 1L # initialise above observed max
    )
  }

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
          export_obj_names = c(
            "y", "lambda_true", "n.chains",
            "n.iter", "n.burnin", "n.thin", "jags.seed"
          )
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
    for (j in 1:J) {
      nll <- nll - dpois(N[j], exp(log_lambda), log = TRUE)
      for (i in 1:I) {
        nll <- nll - dbinom(y[j, i], N[j], plogis(logit_p), log = TRUE)
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
      jags_converged = NA, rtmb_converged = NA, unm_converged = NA
    ))
  }

  sdr <- tryCatch(sdreport(obj), error = function(e) NULL)
  rtmb_converged <- !is.null(sdr) && sdr$pdHess

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
    jags_converged = jags_converged,
    rtmb_converged = rtmb_converged,
    unm_converged = unm_converged
  )
}

run_nmixture <- function(nsim, J, I,
                         lambda_true, p_true,
                         seed,
                         n.chains, n.iter,
                         n.burnin, n.thin) {
  raw <- vector("list", nsim)
  for (s in 1:nsim) {
    raw[[s]] <- fit_all_nmix(seed + s, J, I, lambda_true, p_true,
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

  # JAGS convergence is not required here: this model is hard enough for
  # JAGS that requiring Rhat < 1.1 leaves too few replicates
  ok <- complete.cases(estimates_all) & complete.cases(times_all) &
    rtmb_conv %in% TRUE & unm_conv %in% TRUE
  estimates <- as.data.frame(estimates_all[ok, , drop = FALSE])
  times <- times_all[ok, ]

  conv_rate <- mean(jags_conv, na.rm = TRUE)
  rtmb_conv_rate <- mean(rtmb_conv, na.rm = TRUE)
  unm_conv_rate <- mean(unm_conv, na.rm = TRUE)

  list(
    model = "nmixture",
    estimates = estimates,
    times = times,
    convergence = data.frame(
      sim = seq_len(nsim), jags = jags_conv,
      rtmb = rtmb_conv, unm = unm_conv
    ),
    truth = c(lambda = lambda_true, p = p_true),
    jags_conv_rate = conv_rate,
    rtmb_conv_rate = rtmb_conv_rate,
    unm_conv_rate = unm_conv_rate,
    nsim = nsim, J = J, I = I
  )
}

res <- run_nmixture(
  nsim = nsim, J = J, I = I,
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
cat("RTMB Hessian convergence rate:", round(res$rtmb_conv_rate, 3), "\n")
cat("unmarked Hessian convergence rate:", round(res$unm_conv_rate, 3), "\n")
