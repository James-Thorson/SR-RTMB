
library(RTMB)
library(fmesher)
set.seed(1234)

M <- 200
T <- 3

# true parameters
mu_lambda_true <- log(3.5)
omega_true <- 0.7
gamma_true <- 1.5
p_true <- 0.4
omega_range <- 0.3
omega_SD <- 1
ln_kappa_true <- log(sqrt(8) / omega_range)
ln_tauO_true <- log(1 / (omega_SD * exp(ln_kappa_true) * sqrt(4*pi)))

# -------------------------------------------------------------
# 1. simulation
# -------------------------------------------------------------

coords <- matrix(runif(M * 2), ncol = 2, dimnames = list(NULL, c("x", "y")))

mesh <- fm_mesh_2d(coords, cutoff = 0.05, refine = list(max.edge = c(0.1, 0.3)))
spde <- fm_fem(mesh, order = 2)
A_is <- fm_evaluator(mesh, loc = coords)$proj$A

Q_sim <- exp(4 * ln_kappa_true) * spde$c0 +
  2 * exp(2 * ln_kappa_true) * spde$g1 +
  spde$g2

# draw spatial field on mesh vertices
L <- chol(as.matrix(Q_sim))
omega_v <- backsolve(L, rnorm(mesh$n)) / exp(ln_tauO_true)
omega_i <- as.numeric(A_is %*% omega_v)
lambda_i <- exp(mu_lambda_true + omega_i)

N <- matrix(NA, M, T)
N[, 1] <- rpois(M, lambda_i)
for (t in 1:(T - 1)) {
  N[, t + 1] <- rbinom(M, N[, t], omega_true) + rpois(M, gamma_true)
}
y <- matrix(rbinom(M * T, N, p_true), M, T)

# -------------------------------------------------------------
# 2. fitting
# -------------------------------------------------------------

K <- max(y) * 2
dat <- list(
  y = y, M = M, T = T,
  A_is = A_is, M0 = spde$c0, M1 = spde$g1, M2 = spde$g2
)

par <- list(
  mu_lambda = log(mean(y[, 1]) + 0.1),
  log_gamma = log(1),
  logit_omega = 0,
  logit_p = 0,
  ln_tauO = log(1),
  ln_kappa = log(1),
  omega_s = rep(0, mesh$n),
  SN = matrix(K, nrow = M, ncol = 2*T-1)
)

f <- function(par) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(dat, par, warn = FALSE)
  S <- SN[, seq_len(T-1)]
  N <- SN[, T-1+seq_len(T)]
  gamma <- exp(log_gamma)
  omega <- plogis(logit_omega)
  p <- plogis(logit_p)

  Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
  jnll <- -dgmrf(omega_s,
    mu = 0, Q = Q, log = TRUE,
    scale = 1 / exp(ln_tauO)
  )
  lambda_i <- exp(mu_lambda + A_is %*% omega_s)

  jnll <- jnll - sum(dbinom(S, N[, 1:(T - 1)], omega, log = TRUE), na.rm = TRUE)
  jnll <- jnll - sum(dpois(N[, 1], lambda_i, log = TRUE), na.rm = TRUE)
  for (t in 1:(T - 1)) {
    G <- N[, t + 1] - S[, t]
    jnll <- jnll - sum(dpois(G, gamma, log = TRUE), na.rm = TRUE)
  }
  jnll <- jnll - sum(dbinom(y, size = N, prob = p, log = TRUE), na.rm = TRUE)
  jnll
}

################
# LA in outer step ... works and performs well
################

obj <- MakeADFun(f, par,
  random = c("omega_s", "SN"),
  integrate = list(
    SN = TMB::SR(0:K, discrete = TRUE)
    #omega_s = TMB:::LA()
  ),
  silent = TRUE
)
opt <- nlminb(obj$par, obj$fn, obj$gr,
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
)

################
# LA in inner step ... gives degenerate GMRF estimates, i.e., performs poorly
################

obj2 <- MakeADFun(f, par,
  random = c("omega_s", "SN"),
  integrate = list(
    SN = TMB::SR(0:K, discrete = TRUE),
    omega_s = TMB:::LA()
  ),
  silent = TRUE
)
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr,
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
)

################
# LA using intern
################

obj3 <- MakeADFun(f, par,
  random = c("omega_s", "SN"),
  integrate = list(
    SN = TMB::SR(0:K, discrete = TRUE)
  ),
  silent = TRUE,
  intern = TRUE
)
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr,
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
)


################
# LA using internwith Kasper fix
################

obj4 <- MakeADFun(f, par,
  random = c("omega_s", "SN"),
  integrate = list(
    SN = TMB::SR(0:K, discrete = TRUE)
  ),
  silent = TRUE,
  intern = TRUE,
  inner.control = list(decompose=FALSE)
)
opt4 <- nlminb(obj4$par, obj4$fn, obj4$gr,
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
)
sdr4 = sdreport(obj4)

################
# integrate = LA() using inner.control ... still degenerate
################

obj5 <- MakeADFun(f, par,
  random = c("omega_s", "SN"),
  integrate = list(
    SN = TMB::SR(0:K, discrete = TRUE),
    omega_s = TMB:::LA()
  ),
  silent = TRUE,
  inner.control = list(decompose=FALSE)
)
opt5 <- nlminb(obj5$par, obj5$fn, obj5$gr,
  control = list(eval.max = 1e4, iter.max = 1e4, trace = 1)
)



