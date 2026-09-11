# Sequential Reduction for Discrete-Valued Latent Variables in `RTMB`

Occupancy and abundance (e.g., N-mixture) models require integrating over discrete latent states. Discrete-valued latent variables also arise in other evolutionary models, e.g., state-switching for species traits in phylogenetic trait imputation. This repository implements them in `RTMB` using the Sequential Reduction (SR) algorithm to automatically marginalize those states, and compares parameter recovery and runtime against `unmarked` and `JAGS` and also includes new spatial and mixed-type extensions.

## Quick start

``` bash
git clone https://github.com/James-Thorson/SR-RTMB.git
cd SR-RTMB
make          # install pkgs, run everything and save outputs/figures
make clean    # remove results/ and figures/ outputs
make NSIM=1   # debugging mode; else runs NSIM = 100 (set in Makefile)
### NOTE ###  # the random seed is set in the Makefile and passed to other scripts
```

## Repository structure

```         
SR-RTMB/
├── data/                    # example datasets (e.g. phylogenetic tree, trait data)
├── figures/                 # figures (make plots / make phylo / make spde)
├── models/                  # one script per model: (RTMB, unmarked, JAGS)  results
├── R/                       # shared helpers (install, plotting, SR explainer)
├── results/                 # saved *.rds / timings.csv from each models/*.R run
├── Makefile
├── main.tex                 # manuscript
└── thorson_citations.bib
```

## Demo: discrete latent occupancy model

Here we demonstrate a minimal example of the occupancy model from MacKenzie et al. (2002). We begin with the ecological process of interest (latent occupancy):

$$z_j \sim \text{Bernoulli}(\psi), \quad j = 1, \dots, J$$

and specify how this influences the observed data $y_{ji}$:

$$y_{ji} \mid z_j \sim \text{Bernoulli}(z_j \, p), \quad i = 1, \dots, I$$

where $\psi$ is occupancy probability, $p$ is detection probability, $z_j$ is the true (latent) occupancy state at site $j$, and $y_{ji}$ is detection/non-detection at site $j$, sample $i$.

``` r
# some fake data
# y[j,i] = detection/non-detection
y <- matrix(
  c(
    1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0
  ), 4, 3,
  byrow = TRUE
)
```

**`unmarked`**

First, we implement this model using purpose built software via the `unmarked` package:

``` r
library(unmarked)
# z[j] ~ Bernoulli(psi)
# y[j,i] | z[j] ~ Bernoulli(z[j] * p)
fit <- occu(~1 ~1, unmarkedFrameOccu(y))
print(fit)
```

**`RTMB`**

We now implement the same model using `RTMB`. For an introduction to `RTMB`, see the [introduction vignette](https://cran.r-project.org/web/packages/RTMB/vignettes/RTMB-introduction.html).

``` r
library(RTMB)

# declare data and parameters as tagged lists:
J <- nrow(y)
I <- ncol(y)
dat <- list(y = y, J = J, I = I)
par <- list(logit_psi = 0, logit_p = 0, z = rep(1, J))

# set up the objective function
f <- function(par) {
  getAll(par, dat)                      # works like attach
  psi <- plogis(logit_psi)
  p <- plogis(logit_p)
  jnll <- 0
  for (j in 1:J) {
    # z[j] ~ Bernoulli(psi)
    jnll <- jnll - dbinom(z[j], 1, psi, log = TRUE)
    for (i in 1:I) {
      # y[j,i] | z[j] ~ Bernoulli(z[j] * p)
      jnll <- jnll - dbinom(y[j, i], 1, z[j] * p, log = TRUE)
    }
  }
  jnll
}

# construct objective function with automatic differentiation
obj <- MakeADFun(f, par,
  random = "z",
  integrate = list(
    z = TMB::SR(c(0, 1), discrete = TRUE)
  )
)

opt <- nlminb(obj$par, obj$fn, obj$gr) # run the optimizer
sdr <- sdreport(obj)                   # get uncertainty
print(sdr)
```

While `RTMB` requires more code than `unmarked`, it is able to generalize to more complicated model structures. For example, see Thorson and Kristensen (2024) and the corresponding `RTMB` code for spatial models specified using sparse matrices at [spacetime-ecologist/spacetime-ecologists-RTMB](https://github.com/spacetime-ecologist/spacetime-ecologists-RTMB) or the open N-mixture model with a spatial Gaussian Markov random field (GMRF) specified in the main text.

## Models

`RTMB` code for the following models can be found in the repository:

- `models/occupancy.R` — occupancy (MacKenzie et al. 2002)
- `models/dynamic_occupancy.R` — dynamic occupancy (MacKenzie et al. 2003)
- `models/nmixture.R` — N-mixture (Royle 2004)
- `models/open_nmixture.R` — open N-mixture / Dail-Madsen (2011)
- `models/open_nmixture_spde.R` — open N-mixture with a spatial GMRF (SPDE)
- `models/phylogenetic_mixed_traits.R` — phylogenetic trait imputation

## References

- Dail, D. and Madsen, L. (2011) Models for Estimating Abundance from Repeated Counts of an Open Metapopulation. *Biometrics* 67:577-587.
- MacKenzie, D.I. et al. (2002) Estimating Site Occupancy Rates When Detection Probabilities Are Less Than One. *Ecology* 83:2248-2255.
- MacKenzie, D.I. et al. (2003) Estimating Site Occupancy, Colonization, and Local Extinction When a Species Is Detected Imperfectly. *Ecology* 84:2200-2207.
- Royle, J.A. (2004) N-Mixture Models for Estimating Population Size from Spatially Replicated Counts. *Biometrics* 60:108-115.
- Thorson, J.T. and K. Kristensen 2024. Spatio-temporal Models for Ecologists. 1st edition. Chapman and Hall/CRC, Boca Raton, FL.
