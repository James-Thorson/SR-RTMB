# Empirical MakeADFun() construction time vs. truncation bound K, comparing
# combined (SN) and separate (S, N) integrate() declarations for the open
# N-mixture likelihood at T = 12.
# see R/sequential_reduction.R panels (b)/(c) for the elimination-order
# argument this figure illustrates empirically.

library(RTMB)
library(ggplot2)

set.seed(1)
T <- 12
lambda <- 4
gamma <- 1.5
omega <- 0.8
p <- 0.5

N_true <- numeric(T)
N_true[1] <- rpois(1, lambda)
S_true <- numeric(T - 1)
for (t in 1:(T - 1)) {
  S_true[t] <- rbinom(1, N_true[t], omega)
  G <- rpois(1, gamma)
  N_true[t + 1] <- S_true[t] + G
}
y <- rbinom(T, N_true, p)
dat <- list(y = y, T = T, lambda = lambda, gamma = gamma, omega = omega, p = p)

# one shared likelihood; combined=TRUE splits S,N out of a single SN vector
make_f <- function(combined) {
  function(par) {
    getAll(par, dat)
    if (combined) {
      S <- SN[1:(T - 1)]
      N <- SN[T - 1 + seq_len(T)]
    }
    jnll <- -dpois(N[1], lambda, log = TRUE)
    for (t in 1:(T - 1)) {
      jnll <- jnll - dbinom(S[t], N[t], omega, log = TRUE)
      G <- N[t + 1] - S[t]
      jnll <- jnll - dpois(G, gamma, log = TRUE)
    }
    jnll - sum(dbinom(y, N, p, log = TRUE))
  }
}

par_SN <- list(SN = c(S_true, N_true))
par_S_N <- list(S = S_true, N = N_true)

K_max_values <- seq(from = 5, to = 55, by = 10)
time_SN <- numeric(length(K_max_values))
time_S_N <- numeric(length(K_max_values))

for (i in seq_along(K_max_values)) {
  support <- TMB::SR(0:K_max_values[i], discrete = TRUE)

  time_SN[i] <- system.time(
    MakeADFun(make_f(TRUE), par_SN,
      random = "SN",
      integrate = list(SN = support), silent = TRUE
    )
  )["elapsed"]

  time_S_N[i] <- system.time(
    MakeADFun(make_f(FALSE), par_S_N,
      random = c("S", "N"),
      integrate = list(S = support, N = support), silent = TRUE
    )
  )["elapsed"]

  cat(sprintf(
    "K_max=%3d  SN=%7.4fs  S,N=%8.4fs  ratio=%7.1fx\n",
    K_max_values[i], time_SN[i], time_S_N[i], time_S_N[i] / time_SN[i]
  ))
}

timing <- data.frame(
  K_max = rep(K_max_values, 2),
  seconds = c(time_SN, time_S_N),
  parameterization = rep(c("Combined (SN)", "Separate (S, N)"), each = length(K_max_values))
)
timing$parameterization <- factor(
  timing$parameterization,
  levels = c("Combined (SN)", "Separate (S, N)")
)

col_combined <- "#8E44AD" # endpoint-first, no fill-in (panel c)
col_separate <- "#3f3f3f" # interior-first, fill-in (panel b)

fig <- ggplot(timing, aes(x = K_max, y = seconds, color = parameterization)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = seq(5, 55, by = 5), limits = c(5, 55)) +
  scale_y_log10() +
  scale_color_manual(
    values = c("Combined (SN)" = col_combined, "Separate (S, N)" = col_separate),
    name = NULL
  ) +
  labs(
    x = expression(K[max]),
    y = "MakeADFun() construction time (s, log scale)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
fig
ggsave("figures/elimination_order_timing.png", fig, width = 6.5, height = 5, dpi = 150)
