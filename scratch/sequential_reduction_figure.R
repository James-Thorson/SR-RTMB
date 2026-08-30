# Single-figure explainer for sequential reduction (SR), condensed from
# docs/sequential_reduction_demo.html. Reproduces the site-1 occupancy
# walkthrough (Sections 3-7 of that doc) as three panels:
#   A) Step 1 - tabulate each factor's log-value at z = 0 and z = 1
#   B) Step 2 - cumulative merge (add log-values) as factors fold in
#   C) Step 3 - final merged clique collapsed via log-sum-exp
#
# Usage: Rscript docs/sequential_reduction_figure.R

library(ggplot2)
library(patchwork)

log_sum_exp <- function(x) {
  x_max <- max(x[is.finite(x)])
  x_max + log(sum(exp(x - x_max)))
}

color_absent <- "#9a9890"
color_present <- "#2a78d6"
floor_value <- -10

clip_floor <- function(x, floor = floor_value) {
  pmax(x, floor)
}

# ---- simulate the same occupancy data as the tutorial, trace site 1 ----

set.seed(333)
psi <- 0.6
p <- 0.7
n_sites <- 10
n_visits <- 5

z_true <- rbinom(n_sites, size = 1, prob = psi)
Y <- matrix(
  rbinom(n_sites * n_visits, size = 1, prob = rep(z_true, each = n_visits) * p),
  nrow = n_sites,
  ncol = n_visits,
  byrow = TRUE
)
y_site1 <- Y[1, ]

grid_z <- c(0L, 1L)

# ---- Step 1: one clique per factor ----

clique_A <- dbinom(grid_z, size = 1, prob = psi, log = TRUE)

clique_B <- lapply(seq_len(n_visits), function(j) {
  dbinom(y_site1[j], size = 1, prob = grid_z * p, log = TRUE)
})

factor_names <- c("A", paste0("B", seq_len(n_visits)))
factor_values <- c(list(clique_A), clique_B)

step1_dat <- do.call(rbind, lapply(seq_along(factor_values), function(k) {
  data.frame(
    factor = factor_names[k],
    z = factor(grid_z, levels = c(0, 1), labels = c("z = 0", "z = 1")),
    log_value = factor_values[[k]]
  )
}))
step1_dat$factor <- factor(step1_dat$factor, levels = factor_names)
step1_dat$log_value_clipped <- clip_floor(step1_dat$log_value)
step1_dat$is_impossible <- !is.finite(step1_dat$log_value)

step1_dat$impossible_label <- ifelse(step1_dat$is_impossible, "-Inf", "")

panel_A <- ggplot(step1_dat, aes(x = factor, y = log_value_clipped, fill = z)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = impossible_label, y = log_value_clipped + 0.3),
    position = position_dodge(width = 0.7),
    size = 2.6,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("z = 0" = color_absent, "z = 1" = color_present)) +
  labs(
    title = "A. Tabulate: one clique per factor",
    subtitle = "log-value of each factor at z = 0 and z = 1",
    x = NULL, y = "log-value", fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

# ---- Step 2: cumulative merge as factors fold in ----

cum_z0 <- cumsum(sapply(factor_values, `[`, 1))
cum_z1 <- cumsum(sapply(factor_values, `[`, 2))

step2_dat <- rbind(
  data.frame(factor = factor_names, z = "z = 0", cum_log_value = cum_z0),
  data.frame(factor = factor_names, z = "z = 1", cum_log_value = cum_z1)
)
step2_dat$factor <- factor(step2_dat$factor, levels = factor_names)
step2_dat$cum_log_value_clipped <- clip_floor(step2_dat$cum_log_value)

panel_B <- ggplot(step2_dat, aes(x = factor, y = cum_log_value_clipped, color = z, group = z)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("z = 0" = color_absent, "z = 1" = color_present)) +
  labs(
    title = "B. Merge: fold in factors one at a time",
    subtitle = "cumulative log-value (line clipped at -10)",
    x = NULL, y = "cumulative log-value", color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

# ---- Step 3: integrate out via log-sum-exp ----

merged <- c(z0 = sum(sapply(factor_values, `[`, 1)), z1 = sum(sapply(factor_values, `[`, 2)))
loglik_site1 <- log_sum_exp(merged)

step3_dat <- data.frame(
  z = factor(c("z = 0", "z = 1", "l_i"), levels = c("z = 0", "z = 1", "l_i")),
  log_value = c(clip_floor(merged["z0"]), merged["z1"], loglik_site1),
  label = c("-Inf", round(merged["z1"], 2), round(loglik_site1, 2)),
  kind = c("clique", "clique", "result")
)

panel_C <- ggplot(step3_dat, aes(x = z, y = log_value, fill = kind)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = label), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(clique = color_present, result = "#3a3a38"), guide = "none") +
  labs(
    title = "C. Integrate out: log-sum-exp",
    subtitle = "l_i = log(exp(z=0) + exp(z=1))",
    x = NULL, y = "log-value"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# ---- assemble and save ----

fig <- (panel_A | panel_B | panel_C) +
  plot_annotation(
    title = "Sequential reduction: marginalizing a discrete latent state",
    subtitle = sprintf(
      "Occupancy model, site 1 (visits: %s) - three steps replace a 2^R joint table with R tables of size 2",
      paste(y_site1, collapse = ", ")
    )
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/sequential_reduction.png", fig, width = 12, height = 4.5, dpi = 150)
