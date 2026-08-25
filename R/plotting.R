# Shared plotting/timing helpers. Sourced by plots.R after results/*.rds
# have been written by run_all.R.

library(ggplot2)

col_rtmb <- "#1b9e77"
col_unm <- "#d95f02"
col_jags <- "#7570b3"

# -----------------------------------------------------------
# Timing summary table (mean seconds per fit, by model x framework)
# -----------------------------------------------------------

build_timing_table <- function(res_occ, res_nmix, res_dm) {
  timing <- data.frame(
    Model = c("Occupancy", "N-mixture", "Dail-Madsen"),
    RTMB_mean_s = round(c(
      mean(res_occ$times$rtmb, na.rm = TRUE),
      mean(res_nmix$times$rtmb, na.rm = TRUE),
      mean(res_dm$rtmb_times_per_sim, na.rm = TRUE)
    ), 3),
    unmarked_mean_s = round(c(
      mean(res_occ$times$unm, na.rm = TRUE),
      mean(res_nmix$times$unm, na.rm = TRUE),
      mean(res_dm$unm_times_per_sim, na.rm = TRUE)
    ), 3),
    JAGS_mean_s = round(c(
      mean(res_occ$times$jags, na.rm = TRUE),
      mean(res_nmix$times$jags, na.rm = TRUE),
      mean(res_dm$jags_times_per_sim, na.rm = TRUE)
    ), 3)
  )
  timing$speedup_vs_unm <- round(timing$unmarked_mean_s / timing$RTMB_mean_s, 1)
  timing$speedup_vs_jags <- round(timing$JAGS_mean_s / timing$RTMB_mean_s, 1)
  timing
}

# -----------------------------------------------------------
# Per-sim timing CSV + violin plot
# -----------------------------------------------------------

times_to_long <- function(times_wide, model) {
  data.frame(
    model = model,
    framework = rep(c("RTMB", "unmarked", "JAGS"), each = nrow(times_wide)),
    seconds = c(times_wide$rtmb, times_wide$unm, times_wide$jags)
  )
}

write_timing_csv <- function(res_occ, res_nmix, res_dm, path = "results/timings.csv") {
  dm_times <- data.frame(
    rtmb = res_dm$rtmb_times_per_sim,
    unm  = res_dm$unm_times_per_sim,
    jags = res_dm$jags_times_per_sim
  )
  timings_long <- rbind(
    times_to_long(res_occ$times, "Occupancy"),
    times_to_long(res_nmix$times, "N-mixture"),
    times_to_long(dm_times, "Dail-Madsen")
  )
  write.csv(timings_long, path, row.names = FALSE)
  cat("Per-sim timings saved to", path, "\n")
  timings_long
}

plot_timing_violin <- function(timings_long, path = "plots/timing_violin.png") {
  timings_long$framework <- factor(timings_long$framework, levels = c("RTMB", "unmarked", "JAGS"))
  timings_long$model <- factor(timings_long$model, levels = c("Occupancy", "N-mixture", "Dail-Madsen"))

  fig <- ggplot(timings_long, aes(x = framework, y = seconds, fill = framework)) +
    geom_violin(trim = FALSE, alpha = 0.85, linewidth = 0.3) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.6) +
    scale_y_log10() +
    scale_fill_manual(
      values = c(RTMB = col_rtmb, unmarked = col_unm, JAGS = col_jags),
      guide = "none"
    ) +
    facet_wrap(~model, scales = "free_y") +
    labs(
      x = NULL, y = "Time per fit (s, log scale)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  ggsave(path, fig, width = 10, height = 4.5, dpi = 150)
  cat("Timing violin plot saved to", path, "\n")
}

# -----------------------------------------------------------
# Estimate recovery: relative bias by parameter and framework
# -----------------------------------------------------------

# occupancy/nmixture store estimates as one data.frame with columns
# named "<param>_<rtmb|unm|jags>"
relbias_from_wide <- function(estimates, truth, model_name) {
  frameworks <- c(rtmb = "RTMB", unm = "unmarked", jags = "JAGS")
  rows <- list()
  for (param in names(truth)) {
    for (fw_key in names(frameworks)) {
      est <- estimates[[paste0(param, "_", fw_key)]]
      rows[[length(rows) + 1]] <- data.frame(
        model = model_name,
        parameter = param,
        framework = frameworks[[fw_key]],
        rel_bias = (est - truth[[param]]) / truth[[param]]
      )
    }
  }
  do.call(rbind, rows)
}

# dail-madsen stores estimates as three separate data.frames, one per
# framework, each with columns named after the parameters directly
relbias_from_split <- function(estimates_by_framework, truth, model_name) {
  frameworks <- c(rtmb = "RTMB", unm = "unmarked", jags = "JAGS")
  rows <- list()
  for (param in names(truth)) {
    for (fw_key in names(frameworks)) {
      est <- estimates_by_framework[[fw_key]][[param]]
      rows[[length(rows) + 1]] <- data.frame(
        model = model_name,
        parameter = param,
        framework = frameworks[[fw_key]],
        rel_bias = (est - truth[[param]]) / truth[[param]]
      )
    }
  }
  do.call(rbind, rows)
}

build_relbias_long <- function(res_occ, res_nmix, res_dm) {
  occ <- relbias_from_wide(res_occ$estimates, res_occ$truth, "Occupancy")
  nmix <- relbias_from_wide(res_nmix$estimates, res_nmix$truth, "N-mixture")
  dm <- relbias_from_split(
    list(rtmb = res_dm$estimates_rtmb, unm = res_dm$estimates_unm, jags = res_dm$estimates_jags),
    res_dm$truth, "Dail-Madsen"
  )
  rbind(occ, nmix, dm)
}

plot_estimate_recovery <- function(relbias_long, path = "plots/estimate_recovery.png") {
  relbias_long$framework <- factor(relbias_long$framework, levels = c("RTMB", "unmarked", "JAGS"))
  relbias_long$model <- factor(relbias_long$model, levels = c("Occupancy", "N-mixture", "Dail-Madsen"))
  # "p" is shared by all three models but should sort last within each
  # facet, so it goes last in the global level order too.
  relbias_long$parameter <- factor(relbias_long$parameter, levels = c("psi", "lambda", "gamma", "omega", "p"))

  fig <- ggplot(relbias_long, aes(x = parameter, y = rel_bias, fill = framework)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_violin(
      position = position_dodge(width = 0.8), width = 0.75,
      trim = FALSE, alpha = 0.85, linewidth = 0.3
    ) +
    geom_boxplot(
      position = position_dodge(width = 0.8), width = 0.12,
      outlier.shape = NA, fill = "white", alpha = 0.6
    ) +
    scale_fill_manual(
      values = c(RTMB = col_rtmb, unmarked = col_unm, JAGS = col_jags),
      name = NULL
    ) +
    scale_x_discrete(labels = c(
      psi = expression(psi),
      lambda = expression(lambda),
      gamma = expression(gamma),
      omega = expression(omega),
      p = expression(italic(p))
    )) +
    facet_wrap(~model, scales = "free") +
    labs(x = NULL, y = "Relative bias  (estimate - truth) / truth") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  ggsave(path, fig, width = 12, height = 5.5, dpi = 150)
  cat("Estimate recovery plot saved to", path, "\n")
}
