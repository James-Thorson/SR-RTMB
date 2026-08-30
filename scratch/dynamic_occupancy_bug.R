# !!!NOTE -- this bug has been patched if using latest RTMB!!!
#
# kaskr/RTMB's rolling master branch fixed this issue at aea9b2d —
# "sdreport_xtra: Wrap posfun() around exact variance", made 2026-08-30T13:34:49Z.
#

library(RTMB)
set.seed(1234)
M <- 200
T <- 2
nrep <- 3
psi_true <- 0.6
omega_true <- 0.8
lambda_true <- 0.2
p_true <- 0.5

Z <- matrix(NA, M, T)
Z[, 1] <- rbinom(M, 1, psi_true)
for (t in 1:(T - 1)) {
  survived <- Z[, t] * omega_true
  psi_t <- survived + (1 - survived) * lambda_true
  Z[, t + 1] <- rbinom(M, 1, psi_t)
}

y <- array(NA, dim = c(M, T, nrep))
for (t in 1:T) {
  p_detect <- p_true * Z[, t]
  for (i in 1:nrep) y[, t, i] <- rbinom(M, 1, p_detect)
}

dat <- list(y = y)
par <- list(
  logit_psi = 0, logit_omega = 0, logit_lambda = 0, logit_p = 0,
  Z = matrix(1, nrow = M, ncol = T)
)

f <- function(par) {
  getAll(dat, par, warn = FALSE)
  psi <- plogis(logit_psi)
  omega <- plogis(logit_omega)
  lambda <- plogis(logit_lambda)
  p <- plogis(logit_p)
  jnll <- -sum(dbinom(Z[, 1], 1, psi, log = TRUE))
  survived <- Z[, 1] * omega
  psi_t <- survived + (1 - survived) * lambda
  jnll <- jnll - sum(dbinom(Z[, 2], 1, psi_t, log = TRUE))
  p_detect <- p * Z
  jnll <- jnll - sum(dbinom(y, 1, p_detect, log = TRUE))
  jnll
}

obj <- MakeADFun(f, par, random = "Z", integrate = list(Z = TMB::SR(0:1, discrete = TRUE)))
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
summary(sdr) # returns 4 NaN's...

nan_rows <- which(is.nan(summary(sdr, "random")[, 2]))
nan_sites <- ((nan_rows - 1) %% M) + 1
y[nan_sites, , ] # raw capture history for just those sites
apply(y[nan_sites, , ], 1, any) #  so all animals were detected, which means Z == 1
summary(sdr)[nan_rows + length(opt$par), 1:2] # Z = 1, but SE = NaN?!

# ---- 1. why this shouldn't happen ----
# Z_jt is Bernoulli, so Z_jt in {0, 1} and Z^2 == Z pointwise. That means
#   Var(Z) = E[Z^2] - E[Z]^2 = E[Z] - E[Z]^2 = E[Z] * (1 - E[Z])
# (note: this E[Z] is P(Z=1|data)

# For the sites that fail, Z itself is forced to exactly 1 by a detection, i.e.
# E[Z] = 1 exactly -- not just "close to 1". Plug that in:
#   Var(Z) = 1 * (1 - 1) = 0, exactly, for every Z stuck at Z = 1.
# A NaN here means sqrt() was handed something invalid for a quantity
# that isn't merely non-negative in general -- for THIS Z it should be PROVABLY zero.

# traceback deep dive
options(warn = 2) # change warnings to errors
sdr <- sdreport(obj) # now stops right where "NaNs produced"
traceback() # inspect the call stack

# ---- take a gander at sdreport_xtra()
src <- deparse(RTMB:::sdreport_xtra) # get sdreport_xtra's code
hit <- grep("ans\\$sd <- sqrt", src, fixed = FALSE) # find the line that takes the square root
print(hit)
src[hit] # approximation of total variance

# My understanding of ans$sd <- sqrt(diag.term1 + diag.term2) (line 83)
# diag.term1 = variance contributed by fixed-effect uncertainty
# diag.term2 = the random effect's own variance, given the fitted parameters
# total SE = sqrt(diag.term1 + diag.term2)

# For a normal random effect it uses a Hessian-based formula:
# curvature of the likelihood around the estimate).
# But once you use
# `integrate=` (exact/discrete marginalization, like our Z here), RTMB
# is routing to a DIFFERENT internal function instead: sdreport_xtra().
#   RTMB/R/advector.R:997  xtra1 <- sdreport_xtra(obj, ..., what = "raneffvector")
#   RTMB/R/sdreport.R:77   a_mean <- head(as.vector(obj1$gr(p)), n)  # E[Z],   computed one way
#   RTMB/R/sdreport.R:105  b      <- head(as.vector(obj2$gr(p)), n)  # E[Z^2], computed a SEPARATE way
#   RTMB/R/sdreport.R:115  diag.term2 <- b - a_mean^2                # variance = E[Z^2] - E[Z]^2
#   RTMB/R/sdreport.R:118  ans$sd <- sqrt(diag.term1 + diag.term2)   # then take the square root
# Because Z only takes values 0 or 1 here, Z^2 is always equal to Z --
# so E[Z^2] and E[Z] itself are mathematically the SAME number.
# sdreport_xtra() doesn't know that -- it estimates E[Z] and E[Z^2]
# independently anyway. They can come back bit identical (as they do
# here) and it STILL breaks, because squaring a_mean to compute
# diag.term2 introduces its own extra numerical rounding error.

# ---- can we produce a similar failure, with no RTMB?
a_mean <- 1 + .Machine$double.eps # pick a number ONE teeny bit above 1
b <- a_mean # pretend E[Z^2] came back identical to a_mean like for Z=1
cat("\na_mean =", format(a_mean, digits = 20), "\n") # shows the hidden digits
cat("a_mean^2 =", format(a_mean^2, digits = 20), "\n") # squaring rounds UP one more tiny bit
cat("b - a_mean^2 =", b - a_mean^2, "\n") # this is sdreport.R:115 -- comes out negative
cat("sqrt(...) =", sqrt(b - a_mean^2), "\n") # this is sdreport.R:118 -- sqrt() of a negative number = NaN


###########################################
###########################################
###########################################
# ADDITION FROM JIM 29 AUGUST: WE COULD JUST MAP OFF THESE SPECIFIC Z's VIA map():
map <- list(Z = factor(ifelse(
  apply(y, MARGIN = 1:2, FUN = \(x) any(x == 1)),
  NA,
  seq_len(prod(dim(y)))
)))

# and the model estimates fine...so this is certainly an option if the above isn't worth fixing.
