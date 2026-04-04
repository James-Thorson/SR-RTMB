# run_all_parallel.R
# Parallel version of run_all.R - distributes simulations across cores.
# Each simulation runs on its own core; per-sim seeds are set deterministically
# so results are statistically equivalent to the sequential version.
# Usage: Rscript run_all_parallel.R
#   or:  make parallel

source("models/occupancy.R")
source("models/nmixture.R")
source("models/dail_madsen.R")

library(parallel)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# -----------------------------------------------------------
# Simulation settings - edit here to control all models
# -----------------------------------------------------------
nsim     <- 5          # set to 100 for full run
mc.cores <- detectCores() - 1

# Occupancy
occ_R        <- 200
occ_T        <- 5
occ_psi      <- 0.2
occ_p        <- 0.5
occ_seed     <- 1123

# N-mixture
nmix_R       <- 100
nmix_T       <- 5
nmix_lambda  <- 32
nmix_p       <- 0.25
nmix_seed    <- 333

# Dail-Madsen
dm_M         <- 100
dm_T         <- 5
dm_lambda    <- 4
dm_gamma     <- 1.5
dm_omega     <- 0.8
dm_p         <- 0.5
dm_seed      <- 333

# JAGS MCMC settings - applied to occupancy and N-mixture
jags_chains  <- 4
jags_iter    <- 20000
jags_burnin  <- 10000
jags_thin    <- 1

# Dail-Madsen needs longer chains to converge
dm_jags_iter   <- 50000
dm_jags_burnin <- 25000

# -----------------------------------------------------------
# Parallel wrappers
# Each sim gets seed + s so results are fully reproducible
# -----------------------------------------------------------

run_occupancy_par <- function(nsim, R, T, psi_true, p_true, seed,
                               n.chains, n.iter, n.burnin, n.thin, mc.cores) {
  raw <- mclapply(1:nsim, function(s) {
    set.seed(seed + s)
    fit_all_occ(s, R, T, psi_true, p_true,
                n.chains = n.chains, n.iter = n.iter,
                n.burnin = n.burnin, n.thin = n.thin)
  }, mc.cores = mc.cores)

  estimates <- as.data.frame(na.omit(
    do.call(rbind, lapply(raw, function(x) x$estimates))
  ))
  times <- na.omit(data.frame(
    rtmb = sapply(raw, function(x) x$time_rtmb),
    unm  = sapply(raw, function(x) x$time_unm),
    jags = sapply(raw, function(x) x$time_jags)
  ))
  conv_rate <- mean(sapply(raw, function(x) x$jags_converged), na.rm = TRUE)
  list(model = "occupancy", estimates = estimates, times = times,
       truth = c(psi = psi_true, p = p_true),
       jags_conv_rate = conv_rate, nsim = nsim, R = R, T = T)
}

run_nmixture_par <- function(nsim, R, T, lambda_true, p_true, seed,
                              n.chains, n.iter, n.burnin, n.thin, mc.cores) {
  raw <- mclapply(1:nsim, function(s) {
    set.seed(seed + s)
    fit_all_nmix(s, R, T, lambda_true, p_true,
                 n.chains = n.chains, n.iter = n.iter,
                 n.burnin = n.burnin, n.thin = n.thin)
  }, mc.cores = mc.cores)

  estimates <- as.data.frame(na.omit(
    do.call(rbind, lapply(raw, function(x) x$estimates))
  ))
  times <- na.omit(data.frame(
    rtmb = sapply(raw, function(x) x$time_rtmb),
    unm  = sapply(raw, function(x) x$time_unm),
    jags = sapply(raw, function(x) x$time_jags)
  ))
  conv_rate <- mean(sapply(raw, function(x) x$jags_converged), na.rm = TRUE)
  list(model = "nmixture", estimates = estimates, times = times,
       truth = c(lambda = lambda_true, p = p_true),
       jags_conv_rate = conv_rate, nsim = nsim, R = R, T = T)
}

run_dailmadsen_par <- function(nsim, M, T, lambda_true, gamma_true, omega_true,
                                p_true, seed, n.chains, n.iter, n.burnin, n.thin,
                                mc.cores) {
  datasets <- lapply(1:nsim, function(s) {
    set.seed(seed + s)
    sim_dm(M, T, lambda_true, gamma_true, omega_true, p_true)
  })

  time_rtmb <- system.time({
    res_rtmb <- mclapply(datasets, function(y)
      tryCatch(fit_rtmb(y),
               error = function(e) c(lambda = NA, gamma = NA, omega = NA, p = NA)),
      mc.cores = mc.cores)
  })

  time_unm <- system.time({
    res_unm <- mclapply(datasets, function(y)
      tryCatch(fit_unm(y),
               error = function(e) c(lambda = NA, gamma = NA, omega = NA, p = NA)),
      mc.cores = mc.cores)
  })

  jags_raw <- mclapply(datasets, function(y) {
    t <- system.time({
      out <- fit_jags_dm(y, n.chains = n.chains, n.iter = n.iter,
                         n.burnin = n.burnin, n.thin = n.thin,
                         lambda_true = lambda_true)
    })
    list(estimates = out$estimates, jags_converged = out$jags_converged,
         time_jags = t["elapsed"])
  }, mc.cores = mc.cores)

  res_rtmb <- na.omit(as.data.frame(do.call(rbind, res_rtmb)))
  res_unm  <- na.omit(as.data.frame(do.call(rbind, res_unm)))
  res_jags <- na.omit(as.data.frame(
    do.call(rbind, lapply(jags_raw, function(x) x$estimates))
  ))
  names(res_rtmb) <- c("lambda", "gamma", "omega", "p")
  names(res_unm)  <- c("lambda", "gamma", "omega", "p")
  names(res_jags) <- c("lambda", "gamma", "omega", "p")

  jags_times <- sapply(jags_raw, function(x) x$time_jags)
  conv_rate  <- mean(sapply(jags_raw, function(x) x$jags_converged), na.rm = TRUE)

  list(model = "dailmadsen",
       estimates_rtmb = res_rtmb, estimates_unm = res_unm, estimates_jags = res_jags,
       times = data.frame(
         rtmb = time_rtmb["elapsed"],
         unm  = time_unm["elapsed"],
         jags = sum(jags_times, na.rm = TRUE)
       ),
       jags_times_per_sim = jags_times,
       jags_conv_rate = conv_rate,
       truth = c(lambda = lambda_true, gamma = gamma_true,
                 omega  = omega_true,  p     = p_true),
       nsim = nsim, M = M, T = T)
}

# -----------------------------------------------------------
# Run models
# -----------------------------------------------------------
cat("Running occupancy model (", mc.cores, "cores )...\n")
res_occ <- run_occupancy_par(nsim     = nsim,
                              R        = occ_R,
                              T        = occ_T,
                              psi_true = occ_psi,
                              p_true   = occ_p,
                              seed     = occ_seed,
                              n.chains = jags_chains,
                              n.iter   = jags_iter,
                              n.burnin = jags_burnin,
                              n.thin   = jags_thin,
                              mc.cores = mc.cores)
saveRDS(res_occ, "results/occupancy.rds")
cat("Done.",
    "Mean RTMB time:",     round(mean(res_occ$times$rtmb, na.rm = TRUE), 3), "s |",
    "Mean unmarked time:", round(mean(res_occ$times$unm,  na.rm = TRUE), 3), "s |",
    "Mean JAGS time:",     round(mean(res_occ$times$jags, na.rm = TRUE), 3), "s |",
    "JAGS conv rate:",     round(res_occ$jags_conv_rate, 3), "\n\n")

cat("Running N-mixture model (", mc.cores, "cores )...\n")
res_nmix <- run_nmixture_par(nsim        = nsim,
                              R           = nmix_R,
                              T           = nmix_T,
                              lambda_true = nmix_lambda,
                              p_true      = nmix_p,
                              seed        = nmix_seed,
                              n.chains    = jags_chains,
                              n.iter      = jags_iter,
                              n.burnin    = jags_burnin,
                              n.thin      = jags_thin,
                              mc.cores    = mc.cores)
saveRDS(res_nmix, "results/nmixture.rds")
cat("Done.",
    "Mean RTMB time:",     round(mean(res_nmix$times$rtmb, na.rm = TRUE), 3), "s |",
    "Mean unmarked time:", round(mean(res_nmix$times$unm,  na.rm = TRUE), 3), "s |",
    "Mean JAGS time:",     round(mean(res_nmix$times$jags, na.rm = TRUE), 3), "s |",
    "JAGS conv rate:",     round(res_nmix$jags_conv_rate, 3), "\n\n")

cat("Running Dail-Madsen model (", mc.cores, "cores )...\n")
res_dm <- run_dailmadsen_par(nsim        = nsim,
                              M           = dm_M,
                              T           = dm_T,
                              lambda_true = dm_lambda,
                              gamma_true  = dm_gamma,
                              omega_true  = dm_omega,
                              p_true      = dm_p,
                              seed        = dm_seed,
                              n.chains    = jags_chains,
                              n.iter      = dm_jags_iter,
                              n.burnin    = dm_jags_burnin,
                              n.thin      = jags_thin,
                              mc.cores    = mc.cores)
saveRDS(res_dm, "results/dail_madsen.rds")
cat("Done.",
    "Total RTMB time:",    round(res_dm$times$rtmb, 3), "s |",
    "Total unmarked time:", round(res_dm$times$unm,  3), "s |",
    "Total JAGS time:",    round(res_dm$times$jags, 3), "s |",
    "JAGS conv rate:",     round(res_dm$jags_conv_rate, 3), "\n\n")

# -----------------------------------------------------------
# Timing summary (identical to run_all.R)
# -----------------------------------------------------------
timing <- data.frame(
  Model           = c("Occupancy", "N-mixture", "Dail-Madsen"),
  RTMB_mean_s     = round(c(mean(res_occ$times$rtmb,  na.rm = TRUE),
                             mean(res_nmix$times$rtmb, na.rm = TRUE),
                             res_dm$times$rtmb / res_dm$nsim), 3),
  unmarked_mean_s = round(c(mean(res_occ$times$unm,   na.rm = TRUE),
                             mean(res_nmix$times$unm,  na.rm = TRUE),
                             res_dm$times$unm / res_dm$nsim), 3),
  JAGS_mean_s     = round(c(mean(res_occ$times$jags,  na.rm = TRUE),
                             mean(res_nmix$times$jags, na.rm = TRUE),
                             res_dm$times$jags / res_dm$nsim), 3)
)
timing$speedup_vs_unm  <- round(timing$unmarked_mean_s / timing$RTMB_mean_s, 1)
timing$speedup_vs_jags <- round(timing$JAGS_mean_s     / timing$RTMB_mean_s, 1)

cat("=== Timing Summary ===\n")
print(timing, row.names = FALSE)

cat("\n=== JAGS Convergence Rates (Rhat < 1.1) ===\n")
cat("Occupancy:   ", round(res_occ$jags_conv_rate,  3), "\n")
cat("N-mixture:   ", round(res_nmix$jags_conv_rate, 3), "\n")
cat("Dail-Madsen: ", round(res_dm$jags_conv_rate,   3), "\n")

timing_md <- paste0(
  "| Model | RTMB (s) | unmarked (s) | JAGS (s) | RTMB vs unmarked | RTMB vs JAGS |\n",
  "|---|---|---|---|---|---|\n",
  paste0("| ", timing$Model,
         " | ", timing$RTMB_mean_s,
         " | ", timing$unmarked_mean_s,
         " | ", timing$JAGS_mean_s,
         " | ", timing$speedup_vs_unm,  "x",
         " | ", timing$speedup_vs_jags, "x |",
         collapse = "\n")
)
writeLines(timing_md, "results/timing_summary.md")
cat("\nTiming saved to results/timing_summary.md\n")

readme <- readLines("README.md")
start  <- grep("\\| Occupancy \\|",   readme)[1]
end    <- grep("\\| Dail-Madsen \\|", readme)[1]
if (!is.na(start) && !is.na(end)) {
  readme[start]   <- paste0("| Occupancy | ",   timing$RTMB_mean_s[1],
                             " | ", timing$unmarked_mean_s[1],
                             " | ", timing$JAGS_mean_s[1],
                             " | ", timing$speedup_vs_unm[1],  "x",
                             " | ", timing$speedup_vs_jags[1], "x |")
  readme[start+1] <- paste0("| N-mixture | ",  timing$RTMB_mean_s[2],
                             " | ", timing$unmarked_mean_s[2],
                             " | ", timing$JAGS_mean_s[2],
                             " | ", timing$speedup_vs_unm[2],  "x",
                             " | ", timing$speedup_vs_jags[2], "x |")
  readme[end]     <- paste0("| Dail-Madsen | ", timing$RTMB_mean_s[3],
                             " | ", timing$unmarked_mean_s[3],
                             " | ", timing$JAGS_mean_s[3],
                             " | ", timing$speedup_vs_unm[3],  "x",
                             " | ", timing$speedup_vs_jags[3], "x |")
  writeLines(readme, "README.md")
  cat("README timing table updated\n")
}

# -----------------------------------------------------------
# Plots (identical to run_all.R)
# -----------------------------------------------------------

plot_hist <- function(x, truth, col, main, xlab = "Estimate") {
  hist(x, col = col, border = "white", main = main,
       xlab = xlab, ylab = "Frequency", las = 1)
  abline(v = truth, col = "#DD4444", lwd = 2, lty = 2)
}

col_rtmb <- "#4C72B0"
col_unm  <- "#55A868"
col_jags <- "#CC8400"

# Occupancy - columns = parameters (psi, p), rows = method (RTMB, unmarked, JAGS)
png("figures/occupancy.png", width = 1600, height = 1350, res = 150)
op <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 1), bg = "white", family = "sans")
est_occ <- res_occ$estimates
plot_hist(est_occ$psi_rtmb, res_occ$truth["psi"], col_rtmb,
          expression(bold(paste("RTMB  ", psi))))
plot_hist(est_occ$p_rtmb,   res_occ$truth["p"],   col_rtmb,
          expression(bold(paste("RTMB  ", italic(p)))))
plot_hist(est_occ$psi_unm,  res_occ$truth["psi"], col_unm,
          expression(bold(paste("unmarked  ", psi))))
plot_hist(est_occ$p_unm,    res_occ$truth["p"],   col_unm,
          expression(bold(paste("unmarked  ", italic(p)))))
plot_hist(est_occ$psi_jags, res_occ$truth["psi"], col_jags,
          expression(bold(paste("JAGS  ", psi))))
plot_hist(est_occ$p_jags,   res_occ$truth["p"],   col_jags,
          expression(bold(paste("JAGS  ", italic(p)))))
par(op)
dev.off()

# N-mixture - columns = parameters (lambda, p), rows = method (RTMB, unmarked, JAGS)
png("figures/nmixture.png", width = 1600, height = 1350, res = 150)
op <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 1), bg = "white", family = "sans")
est_nmix <- res_nmix$estimates
plot_hist(est_nmix$lambda_rtmb, res_nmix$truth["lambda"], col_rtmb,
          expression(bold(paste("RTMB  ", lambda))))
plot_hist(est_nmix$p_rtmb,      res_nmix$truth["p"],      col_rtmb,
          expression(bold(paste("RTMB  ", italic(p)))))
plot_hist(est_nmix$lambda_unm,  res_nmix$truth["lambda"], col_unm,
          expression(bold(paste("unmarked  ", lambda))))
plot_hist(est_nmix$p_unm,       res_nmix$truth["p"],      col_unm,
          expression(bold(paste("unmarked  ", italic(p)))))
plot_hist(est_nmix$lambda_jags, res_nmix$truth["lambda"], col_jags,
          expression(bold(paste("JAGS  ", lambda))))
plot_hist(est_nmix$p_jags,      res_nmix$truth["p"],      col_jags,
          expression(bold(paste("JAGS  ", italic(p)))))
par(op)
dev.off()

# Dail-Madsen - columns = parameters (lambda, gamma, omega, p), rows = method
png("figures/dail_madsen.png", width = 2400, height = 1350, res = 150)
op <- par(mfrow = c(3, 4), mar = c(4, 4, 3, 1), bg = "white", family = "sans")
est_dm_rtmb <- res_dm$estimates_rtmb
est_dm_unm  <- res_dm$estimates_unm
est_dm_jags <- res_dm$estimates_jags

# Row 1: RTMB
plot_hist(est_dm_rtmb$lambda, res_dm$truth["lambda"], col_rtmb,
          expression(bold(paste("RTMB  ", lambda))))
plot_hist(est_dm_rtmb$gamma,  res_dm$truth["gamma"],  col_rtmb,
          expression(bold(paste("RTMB  ", gamma))))
plot_hist(est_dm_rtmb$omega,  res_dm$truth["omega"],  col_rtmb,
          expression(bold(paste("RTMB  ", omega))))
plot_hist(est_dm_rtmb$p,      res_dm$truth["p"],      col_rtmb,
          expression(bold(paste("RTMB  ", italic(p)))))

# Row 2: unmarked
plot_hist(est_dm_unm$lambda,  res_dm$truth["lambda"], col_unm,
          expression(bold(paste("unmarked  ", lambda))))
plot_hist(est_dm_unm$gamma,   res_dm$truth["gamma"],  col_unm,
          expression(bold(paste("unmarked  ", gamma))))
plot_hist(est_dm_unm$omega,   res_dm$truth["omega"],  col_unm,
          expression(bold(paste("unmarked  ", omega))))
plot_hist(est_dm_unm$p,       res_dm$truth["p"],      col_unm,
          expression(bold(paste("unmarked  ", italic(p)))))

# Row 3: JAGS
plot_hist(est_dm_jags$lambda, res_dm$truth["lambda"], col_jags,
          expression(bold(paste("JAGS  ", lambda))))
plot_hist(est_dm_jags$gamma,  res_dm$truth["gamma"],  col_jags,
          expression(bold(paste("JAGS  ", gamma))))
plot_hist(est_dm_jags$omega,  res_dm$truth["omega"],  col_jags,
          expression(bold(paste("JAGS  ", omega))))
plot_hist(est_dm_jags$p,      res_dm$truth["p"],      col_jags,
          expression(bold(paste("JAGS  ", italic(p)))))

par(op)
dev.off()

cat("\nAll figures saved to figures/\n")
