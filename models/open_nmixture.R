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


library(RTMB)
library(unmarked)
library(R2jags)

# run parameters
# set by `make NSIM=X`; set nsim manually here instead if running standalone
nsim <- as.integer(Sys.getenv("NSIM"))
# set by `make SEED=X`; set seed manually here instead if running standalone
seed <- as.integer(Sys.getenv("SEED"))
M <- 100
T <- 5
lambda_true <- 4
gamma_true <- 1.5
omega_true <- 0.8
p_true <- 0.5
n.chains <- 4
n.iter <- 20000
n.burnin <- 10000
n.thin <- 1

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

  # jags.parallel() defaults jags.seed to a hardcoded 123, which - unlike
  # sequential jags() inheriting R's naturally-advancing RNG - makes every
  # chain's starting values IDENTICAL across every replicate and every run,
  # since seeds <- jags.seed + seq_len(n.chains) never changes. Draw a fresh
  # one from R's already-seeded stream so inits actually vary per replicate.
  jags.seed <- sample.int(1e6, 1)

  # JAGS samples the full latent state explicitly here (no SR), making it
  # the slowest fit in the pipeline - run its chains in parallel, one per
  # core, instead of sequentially in a single process.
  # jags.parallel() re-resolves these names on each worker (its internal
  # .runjags() does eval(expression(n.iter)) etc., and inits() is called
  # there too) rather than just inheriting fit_jags_dm's closure, so they
  # all need to be listed explicitly here or the workers fail to find them.
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
  M <- nrow(y)
  T <- ncol(y)
  K <- max(y) * 2
  dat <- list(y = y, M = M, T = T)

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
    SN = matrix(K, nrow = M, ncol = 2 * T - 1)
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

run_dailmadsen <- function(nsim, M, T,
                           lambda_true, gamma_true,
                           omega_true, p_true,
                           seed,
                           n.chains, n.iter,
                           n.burnin, n.thin) {
  datasets <- lapply(1:nsim, function(s) {
    set.seed(seed + s)
    sim_dm(M, T, lambda_true, gamma_true, omega_true, p_true)
  })

  # Each dataset timed individually for every framework, so per-sim times
  # are always available and mean times are true per-fit times rather than
  # wall-clock-divided-by-nsim.
  unm_raw <- vector("list", nsim)
  for (i in seq_along(datasets)) {
    y <- datasets[[i]]
    t <- system.time({
      out <- tryCatch(fit_unm(y),
        error = function(e) c(lambda = NA, gamma = NA, omega = NA, p = NA)
      )
    })
    unm_raw[[i]] <- list(estimate = out, time = t["elapsed"])
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
        error = function(e) c(lambda = NA, gamma = NA, omega = NA, p = NA)
      )
    })
    rtmb_raw[[i]] <- list(estimate = out, time = t["elapsed"])
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

  # Drop the same simulation indices from all three frameworks so estimates
  # always correspond to the same dataset - critical for fair comparison
  ok <- complete.cases(res_rtmb) & complete.cases(res_unm) & complete.cases(res_jags)
  res_rtmb <- res_rtmb[ok, ]
  res_unm <- res_unm[ok, ]
  res_jags <- res_jags[ok, ]

  conv_rate <- mean(sapply(jags_raw, function(x) x$jags_converged), na.rm = TRUE)

  list(
    model = "dailmadsen",
    estimates_unm = res_unm,
    estimates_jags = res_jags,
    estimates_rtmb = res_rtmb,
    times = data.frame(
      rtmb = sum(rtmb_times, na.rm = TRUE),
      unm  = sum(unm_times, na.rm = TRUE),
      jags = sum(jags_times, na.rm = TRUE)
    ),
    unm_times_per_sim = unm_times,
    jags_times_per_sim = jags_times,
    rtmb_times_per_sim = rtmb_times,
    jags_conv_rate = conv_rate,
    truth = c(
      lambda = lambda_true, gamma = gamma_true,
      omega = omega_true, p = p_true
    ),
    nsim = nsim, M = M, T = T
  )
}

res <- run_dailmadsen(
  nsim = nsim, M = M, T = T,
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
