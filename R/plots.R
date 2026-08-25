# plots.R
# Loads results/*.rds (written by run_all.R) and generates all figures.
# Usage: Rscript R/plots.R   or   make plots
source("R/plotting.R")

dir.create("figures", showWarnings = FALSE)

res_occ <- readRDS("results/occupancy.rds")
res_nmix <- readRDS("results/nmixture.rds")
res_dm <- readRDS("results/dail_madsen.rds")

timing <- build_timing_table(res_occ, res_nmix, res_dm)
cat("=== Timing Summary ===\n")
print(timing, row.names = FALSE)
cat("\n=== JAGS Convergence Rates (Rhat < 1.1) ===\n")
cat("Occupancy:   ", round(res_occ$jags_conv_rate, 3), "\n")
cat("N-mixture:   ", round(res_nmix$jags_conv_rate, 3), "\n")
cat("Dail-Madsen: ", round(res_dm$jags_conv_rate, 3), "\n")

timings_long <- write_timing_csv(res_occ, res_nmix, res_dm, "results/timings.csv")
plot_timing_violin(timings_long, "figures/timing_violin.png")

relbias_long <- build_relbias_long(res_occ, res_nmix, res_dm)
plot_estimate_recovery(relbias_long, "figures/estimate_recovery.png")

cat("\nAll figures saved to figures/\n")
