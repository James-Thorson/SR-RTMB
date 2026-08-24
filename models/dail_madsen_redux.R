
library(RTMB)

nsim = 1
M = 50
T = 5
lambda = 4
gamma = 1.5
omega = 0.8
p = 0.5
seed = 333,


# Simulator
N <- matrix(NA, M, T)
N[, 1] <- rpois(M, lambda)
for (t in 1:(T - 1)) {
  S <- rbinom(M, N[, t], omega)
  G <- rpois(M, gamma)
  N[, t + 1] <- S + G
}
y = matrix(rbinom(M * T, N, p), M, T)

# Estimator
K <- max(y) * 2
dat <- list(y = y, M = M, T = T, K = K)

f <- function(par) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(par, dat)
  lambda <- exp(log_lambda)
  gamma <- exp(log_gamma)
  omega <- plogis(logit_omega)
  p <- plogis(logit_p)
  #N = y + Nextra
  #jnll = -1 * sum(dpois(G, gamma, log = TRUE), na.rm= TRUE)
  jnll = -1 * sum(dbinom(S, N[,1:(T-1)], omega, log = TRUE), na.rm= TRUE)
  jnll = jnll - sum(dpois( N[,1], lambda, log = TRUE), na.rm= TRUE)
  for (t in 1:(T - 1)) {
    #S <- N[, t+1] - G[, t]
    #jnll = jnll - sum(dbinom(S, N[,t], omega, log = TRUE), na.rm= TRUE)
    G = N[, t+1] - S[, t]
    jnll = jnll - sum(dpois(G, gamma, log = TRUE), na.rm= TRUE)
  }
  jnll = jnll - sum(dbinom(y, size = N, prob = p, log = TRUE), na.rm= TRUE)
  return(jnll)
}

par <- list(
  log_lambda  = log(mean(y[, 1]) + 0.1),
  log_gamma   = log(1.5),
  logit_omega = 0,
  logit_p     = 0,
  #Nextra = matrix(K, nrow = M, ncol = T ),
  #G = matrix(K, nrow = M, ncol = T-1 ),
  S = matrix(K, nrow = M, ncol = T-1 ),
  N = matrix(K, nrow = M, ncol = T )
)
f(par)

obj <- MakeADFun(
  f,
  par,
  #random = c("N", "G"),
  random = c("N", "S"),
  integrate = list(
    #G = TMB::SR( 0:K, discrete=TRUE),
    S = TMB::SR( 0:K, discrete=TRUE),
    N = TMB::SR( 0:K, discrete=TRUE)
  )
)
obj$fn(obj$par)

opt <- nlminb(
  obj$par, obj$fn, obj$gr,
  control = list(iter.max = 1e5, eval.max = 1e5)
)
sdr = sdreport(obj)
