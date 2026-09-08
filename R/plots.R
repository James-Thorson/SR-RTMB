# plots.R
# Loads results/*.rds (written by models/occupancy.R, dynamic_occupancy.R,
# nmixture.R, open_nmixture.R) and generates all figures.
# Usage: Rscript R/plots.R   or   make plots
source("R/plotting.R")

res_occ <- readRDS("results/occupancy.rds")
res_dynocc <- readRDS("results/dynamic_occupancy.rds")
res_nmix <- readRDS("results/nmixture.rds")
res_om <- readRDS("results/open_nmixture.rds")

models <- list(
  Occupancy = list(
    kind = "wide",
    estimates = res_occ$estimates,
    times = list(rtmb = res_occ$times$rtmb, unm = res_occ$times$unm, jags = res_occ$times$jags),
    truth = res_occ$truth,
    nsim = res_occ$nsim,
    jags_conv_rate = res_occ$jags_conv_rate,
    rtmb_conv_rate = res_occ$rtmb_conv_rate,
    unm_conv_rate = res_occ$unm_conv_rate
  ),
  "Dynamic Occupancy" = list(
    kind = "wide",
    estimates = res_dynocc$estimates,
    times = list(rtmb = res_dynocc$times$rtmb, unm = res_dynocc$times$unm, jags = res_dynocc$times$jags),
    truth = res_dynocc$truth,
    nsim = res_dynocc$nsim,
    jags_conv_rate = res_dynocc$jags_conv_rate,
    rtmb_conv_rate = res_dynocc$rtmb_conv_rate,
    unm_conv_rate = res_dynocc$unm_conv_rate
  ),
  "N-mixture" = list(
    kind = "wide",
    estimates = res_nmix$estimates,
    times = list(rtmb = res_nmix$times$rtmb, unm = res_nmix$times$unm, jags = res_nmix$times$jags),
    truth = res_nmix$truth,
    nsim = res_nmix$nsim,
    jags_conv_rate = res_nmix$jags_conv_rate,
    rtmb_conv_rate = res_nmix$rtmb_conv_rate,
    unm_conv_rate = res_nmix$unm_conv_rate
  ),
  "Open N-mixture" = list(
    kind = "split",
    estimates = list(rtmb = res_om$estimates_rtmb, unm = res_om$estimates_unm, jags = res_om$estimates_jags),
    times = list(rtmb = res_om$rtmb_times_per_sim, unm = res_om$unm_times_per_sim, jags = res_om$jags_times_per_sim),
    truth = res_om$truth,
    nsim = res_om$nsim,
    jags_conv_rate = res_om$jags_conv_rate,
    rtmb_conv_rate = res_om$rtmb_conv_rate,
    unm_conv_rate = res_om$unm_conv_rate
  )
)

timing <- build_timing_table(models)
cat("=== Timing Summary ===\n")
print(timing, row.names = FALSE)
cat("\n=== Convergence per Model and Framework ===\n")
print(build_convergence_table(models), row.names = FALSE)

timings_long <- write_timing_csv(models, "results/timings.csv")
plot_timing_violin(timings_long, "figures/timing_violin.png")

relbias_long <- build_relbias_long(models)
plot_estimate_recovery(relbias_long, "figures/estimate_recovery.png")

cat("\nAll figures saved to figures/\n")
