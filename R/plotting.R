# Shared plotting/timing helpers. Sourced by R/plots.R after results/*.rds
# have been written by the model scripts in models/.

library(ggplot2)

framework_colors <- c(RTMB = "#8E44AD", unmarked = "#2C3E50", JAGS = "#E67E22")

# -----------------------------------------------------------
# `models` is a named list, one entry per model, in the order they
# should appear in tables/facets. Each entry has:
#   kind           - "wide"  (one estimates data.frame, columns "<param>_<rtmb|unm|jags>",
#                    like occupancy/n-mixture/dynamic-occupancy)
#                  - "split" (list(rtmb=, unm=, jags=) of per-framework
#                    data.frames, like open n-mixture / dail-madsen)
#   estimates      - shaped per `kind`
#   times          - list(rtmb = <vec>, unm = <vec>, jags = <vec>) of
#                    per-replicate seconds, regardless of `kind`
#   truth          - named truth vector
#   nsim           - replicates requested (before any dropped for NAs)
#   jags_conv_rate - fraction of replicates with JAGS Rhat < 1.1
#   rtmb_conv_rate - fraction of replicates with an invertible RTMB Hessian
#   unm_conv_rate  - fraction of replicates with an invertible unmarked Hessian
# -----------------------------------------------------------

# -----------------------------------------------------------
# Timing summary table (mean seconds per fit, by model x framework)
# -----------------------------------------------------------

build_timing_table <- function(models) {
  timing <- data.frame(
    Model = names(models),
    RTMB_mean_s = round(sapply(models, function(m) mean(m$times$rtmb, na.rm = TRUE)), 3),
    unmarked_mean_s = round(sapply(models, function(m) mean(m$times$unm, na.rm = TRUE)), 3),
    JAGS_mean_s = round(sapply(models, function(m) mean(m$times$jags, na.rm = TRUE)), 3),
    row.names = NULL
  )
  timing$speedup_vs_unm <- round(timing$unmarked_mean_s / timing$RTMB_mean_s, 1)
  timing$speedup_vs_jags <- round(timing$JAGS_mean_s / timing$RTMB_mean_s, 1)
  timing
}

# -----------------------------------------------------------
# Convergence summary table (fraction converged, by model x framework)
# -----------------------------------------------------------

build_convergence_table <- function(models) {
  data.frame(
    Model = names(models),
    RTMB_conv_rate = round(sapply(models, function(m) m$rtmb_conv_rate), 3),
    unmarked_conv_rate = round(sapply(models, function(m) m$unm_conv_rate), 3),
    JAGS_conv_rate = round(sapply(models, function(m) m$jags_conv_rate), 3),
    row.names = NULL
  )
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

build_timing_long <- function(models) {
  do.call(rbind, lapply(names(models), function(name) {
    times_to_long(as.data.frame(models[[name]]$times), name)
  }))
}

write_timing_csv <- function(models, path = "results/timings.csv") {
  timings_long <- build_timing_long(models)
  write.csv(timings_long, path, row.names = FALSE)
  cat("Per-sim timings saved to", path, "\n")
  timings_long
}

plot_timing_violin <- function(timings_long, path = "figures/timing_violin.png", title = NULL) {
  timings_long$framework <- factor(timings_long$framework, levels = c("RTMB", "unmarked", "JAGS"))
  timings_long$model <- factor(timings_long$model, levels = unique(timings_long$model))

  fig <- ggplot(timings_long, aes(x = framework, y = seconds, fill = framework)) +
    geom_violin(trim = FALSE, alpha = 0.85, linewidth = 0.3) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.9) +
    scale_y_log10() +
    scale_fill_manual(
      values = framework_colors,
      guide = "none"
    ) +
    facet_wrap(~model, nrow = 2, dir = "v", scales = "free_y") +
    labs(
      x = NULL, y = "Time per fit (s, log scale)", title = title
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  ggsave(path, fig, width = 11, height = 7.5, dpi = 150)
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

build_relbias_long <- function(models) {
  do.call(rbind, lapply(names(models), function(name) {
    m <- models[[name]]
    if (m$kind == "wide") {
      relbias_from_wide(m$estimates, m$truth, name)
    } else {
      relbias_from_split(m$estimates, m$truth, name)
    }
  }))
}

plot_estimate_recovery <- function(relbias_long, path = "figures/estimate_recovery.png", title = NULL) {
  relbias_long$framework <- factor(relbias_long$framework, levels = c("RTMB", "unmarked", "JAGS"))

  model_order <- unique(relbias_long$model)
  relbias_long$model <- factor(relbias_long$model, levels = model_order)

  # Panels fill column-by-column below (facet_wrap(dir = "v"), nrow =
  # facet_nrow), so pair up models within the same column: the top model's
  # parameter order (its own order of first appearance) is the reference,
  # and the bottom model repeats its shared parameters in that same order
  # before appending whatever parameters are new to it.
  facet_nrow <- 2
  natural_order <- function(model_name) {
    unique(relbias_long$parameter[relbias_long$model == model_name])
  }
  param_order_by_model <- list()
  for (i in seq_along(model_order)) {
    m <- as.character(model_order[i])
    if (i %% facet_nrow == 1) {
      param_order_by_model[[m]] <- natural_order(m)
    } else {
      top_order <- param_order_by_model[[as.character(model_order[i - 1])]]
      own_order <- natural_order(m)
      param_order_by_model[[m]] <- c(intersect(top_order, own_order), setdiff(own_order, top_order))
    }
  }

  sep <- "::"
  panel_key_levels <- unlist(lapply(model_order, function(m) {
    paste(m, param_order_by_model[[as.character(m)]], sep = sep)
  }))
  relbias_long$panel_key <- factor(
    paste(relbias_long$model, relbias_long$parameter, sep = sep),
    levels = panel_key_levels
  )

  param_expr <- c(psi = "psi", omega = "omega", lambda = "lambda", gamma = "gamma", p = "italic(p)")
  key_params <- sub(paste0("^.*", sep), "", panel_key_levels)
  key_label_lookup <- as.expression(lapply(param_expr[key_params], function(e) parse(text = e)[[1]]))
  names(key_label_lookup) <- panel_key_levels
  # facet_wrap(scales = "free") shows only the levels present in each panel,
  # so labels must be looked up by name (not positionally) or they end up
  # misaligned with what's actually displayed in that panel.
  key_labels <- function(breaks) key_label_lookup[breaks]

  # geom_boxplot's own dodge computation drops all but one dodge slot when
  # every x level is shared across all fill groups (as it is here), so the
  # group aesthetic must be set explicitly or most boxes silently vanish.
  fig <- ggplot(relbias_long, aes(
    x = panel_key, y = rel_bias, fill = framework,
    group = interaction(panel_key, framework)
  )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_violin(
      position = position_dodge(width = 0.8), width = 0.75,
      trim = FALSE, alpha = 0.85, linewidth = 0.3
    ) +
    geom_boxplot(
      position = position_dodge(width = 0.8), width = 0.12,
      outlier.shape = NA, fill = "white", alpha = 0.9
    ) +
    scale_fill_manual(
      values = framework_colors,
      name = NULL
    ) +
    scale_x_discrete(labels = key_labels) +
    facet_wrap(~model, nrow = facet_nrow, dir = "v", scales = "free") +
    labs(x = NULL, y = "Relative bias  (estimate - truth) / truth", title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  ggsave(path, fig, width = 11, height = 8.5, dpi = 150)
  cat("Estimate recovery plot saved to", path, "\n")
}
