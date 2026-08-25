# run_all.R
# Runs all three models and saves results/*.rds (estimates + per-sim timings).
# Usage:
#   Rscript run_all.R            # sequential (mc.cores = 1)
#   Rscript run_all.R parallel   # parallel (mc.cores = detectCores() - 1)
#   or: make test / make run
# Follow with `make plots` (Rscript plots.R) to generate figures.
source("models/occupancy.R")
source("models/nmixture.R")
source("models/dail_madsen.R")

args <- commandArgs(trailingOnly = TRUE)
mc.cores <- if (identical(args[1], "parallel")) parallel::detectCores() - 1 else 1

Sys.setenv(OMP_NUM_THREADS = "4")
dir.create("results", showWarnings = FALSE)
# -----------------------------------------------------------
# Simulation settings - edit here to control all models
# -----------------------------------------------------------
nsim <- 100
# Occupancy
occ_R <- 200
occ_T <- 5
occ_psi <- 0.2
occ_p <- 0.5
occ_seed <- 1123
# N-mixture
nmix_R <- 100
nmix_T <- 5
nmix_lambda <- 32
nmix_p <- 0.25
nmix_seed <- 333
# Dail-Madsen
dm_M <- 100
dm_T <- 5
dm_lambda <- 4
dm_gamma <- 1.5
dm_omega <- 0.8
dm_p <- 0.5
dm_seed <- 333
# JAGS MCMC settings - applied to occupancy and N-mixture
jags_chains <- 4
jags_iter <- 40000
jags_burnin <- 20000
jags_thin <- 1
# Dail-Madsen needs longer chains to converge
dm_jags_iter <- 60000
dm_jags_burnin <- 50000
# -----------------------------------------------------------
# Run models
# -----------------------------------------------------------
cat("Running occupancy model (mc.cores =", mc.cores, ")...\n")
res_occ <- run_occupancy(
  nsim = nsim, R = occ_R, T = occ_T,
  psi_true = occ_psi, p_true = occ_p, seed = occ_seed, mc.cores = mc.cores,
  n.chains = jags_chains, n.iter = jags_iter,
  n.burnin = jags_burnin, n.thin = jags_thin
)
saveRDS(res_occ, "results/occupancy.rds")
cat(
  "Done.",
  "Mean RTMB time:", round(mean(res_occ$times$rtmb, na.rm = TRUE), 3), "s |",
  "Mean unmarked time:", round(mean(res_occ$times$unm, na.rm = TRUE), 3), "s |",
  "Mean JAGS time:", round(mean(res_occ$times$jags, na.rm = TRUE), 3), "s |",
  "JAGS conv rate:", round(res_occ$jags_conv_rate, 3), "\n\n"
)
cat("Running N-mixture model (mc.cores =", mc.cores, ")...\n")
res_nmix <- run_nmixture(
  nsim = nsim, R = nmix_R, T = nmix_T,
  lambda_true = nmix_lambda, p_true = nmix_p, seed = nmix_seed, mc.cores = mc.cores,
  n.chains = jags_chains, n.iter = jags_iter,
  n.burnin = jags_burnin, n.thin = jags_thin
)
saveRDS(res_nmix, "results/nmixture.rds")
cat(
  "Done.",
  "Mean RTMB time:", round(mean(res_nmix$times$rtmb, na.rm = TRUE), 3), "s |",
  "Mean unmarked time:", round(mean(res_nmix$times$unm, na.rm = TRUE), 3), "s |",
  "Mean JAGS time:", round(mean(res_nmix$times$jags, na.rm = TRUE), 3), "s |",
  "JAGS conv rate:", round(res_nmix$jags_conv_rate, 3), "\n\n"
)
cat("Running Dail-Madsen model (mc.cores =", mc.cores, ")...\n")
res_dm <- run_dailmadsen(
  nsim = nsim, M = dm_M, T = dm_T,
  lambda_true = dm_lambda, gamma_true = dm_gamma,
  omega_true = dm_omega, p_true = dm_p, seed = dm_seed, mc.cores = mc.cores,
  n.chains = jags_chains, n.iter = dm_jags_iter,
  n.burnin = dm_jags_burnin, n.thin = jags_thin
)
saveRDS(res_dm, "results/dail_madsen.rds")
cat(
  "Done.",
  "Mean RTMB time:", round(mean(res_dm$rtmb_times_per_sim, na.rm = TRUE), 3), "s |",
  "Mean unmarked time:", round(mean(res_dm$unm_times_per_sim, na.rm = TRUE), 3), "s |",
  "Mean JAGS time:", round(mean(res_dm$jags_times_per_sim, na.rm = TRUE), 3), "s |",
  "JAGS conv rate:", round(res_dm$jags_conv_rate, 3), "\n\n"
)
cat("All models complete. Results saved to results/. Run `make plots` to generate figures.\n")
