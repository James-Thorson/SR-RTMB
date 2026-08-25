# Discrete Latent Variable Models in RTMB

Occupancy, N-mixture, and open population models such as the Dail-Madsen are among the most widely used statistical tools in wildlife ecology, yet all share a common challenge: inference requires integrating over discrete latent states that are never directly observed. Automatic differentiation frameworks like TMB and its derivatives have historically struggled with this class of models. This has largely excluded such models from the ADMB/TMB/RTMB ecosystems. Here we implement all three models in RTMB using sequential reduction to marginalize over discrete latent states, and compare parameter recovery and computational performance against `unmarked`, a widely used package purpose-built for this class of models, and `JAGS`, a general-purpose Bayesian sampler.

## Quick Start

``` bash
git clone https://github.com/James-Thorson/RTMBdre.git
cd RTMBdre
make install        # install required R packages (system JAGS must be installed first)
time make           # run
```

This runs the simulation and saves results to `results/`. Run `make plots` afterward to generate figures into `figures/`.

## Models

All discrete latent states below are marginalized via Sequential Reduction (SR) in RTMB; JAGS samples them explicitly via MCMC for comparison.

### Occupancy (MacKenzie et al. 2002)

```         
z_i ~ Bernoulli(psi)
y_ij | z_i ~ Bernoulli(z_i * p)
```

### N-mixture (Royle 2004)

```         
N_i ~ Poisson(lambda)
y_ij | N_i ~ Binomial(N_i, p)
```

### Dail-Madsen (Dail & Madsen 2011)

```         
N_i1 ~ Poisson(lambda)
S_it | N_it ~ Binomial(N_it, omega)
G_it ~ Poisson(gamma)
N_it+1 = S_it + G_it
y_it | N_it ~ Binomial(N_it, p)
```

## Results

### Estimate recovery

Relative bias, `(estimate - truth) / truth`, across simulated datasets, one panel per model and one violin/boxplot per parameter, colored by framework. For RTMB and unmarked these are maximum likelihood estimates; for JAGS these are posterior means. The dashed line at 0 marks unbiased recovery.

![Estimate recovery](figures/estimate_recovery.png)

### Runtime

![Timing](figures/timing_violin.png)

### Convergence

Convergence is assessed via the Gelman-Rubin R-hat statistic (Gelman & Rubin 1992) across 4 chains; a fit is flagged converged when R-hat \< 1.1 for all monitored parameters. Rates are printed to console by `make plots`. Dail-Madsen's lower convergence rate reflects the difficulty of sampling its high-dimensional latent state space directly.

## JAGS MCMC Settings

JAGS settings are configurable at the top of `run_all.R`. Occupancy and N-mixture use shorter chains; Dail-Madsen requires longer runs due to the complexity of its latent state space.

## Repository Structure

```         
models/
  occupancy.R                    # occupancy model; JAGS comparison embedded
  nmixture.R                     # N-mixture model; JAGS comparison embedded
  dail_madsen.R                  # Dail-Madsen model; JAGS comparison embedded
  dail_madsen_spde.R             # spatial Dail-Madsen with GMRF on lambda via SPDE (in progress)
  phylogenetic_mixed_traits.R    # phylogenetic trait imputation (viviparity/body size)
R/
  utils.R                 # shared helpers (e.g. lapply_maybe) sourced by model scripts
  plotting.R               # shared timing/plotting helpers sourced by plots.R
  run_all.R                # runs all models, saves results/*.rds
  plots.R                  # loads results/*.rds, saves figures to figures/
  install.R                # installs required R packages
figures/                # figures (created by plots.R and phylogenetic_mixed_traits.R)
results/                # .rds results and per-sim timing CSV (created by run_all.R)
Makefile
```

Each model file contains RTMB, unmarked, and JAGS fits in a single internal worker function (`fit_all_occ`, `fit_all_nmix`, `fit_jags_dm`) that operates on the same simulated dataset. This ensures all timing comparisons are made on identical data.

## Usage

``` bash
make install    # install required R packages (once)
make            # run all three models in parallel, save results/*.rds
make plots      # build figures from results/*.rds into figures/
```

Other targets:

| Target | Effect |
|------------------------------------|------------------------------------|
| `make test` | run all models sequentially instead of in parallel |
| `make occupancy`, `make nmixture`, `make dailmadsen` | run a single model |
| `make spde` | run the spatial Dail-Madsen extension (in progress) |
| `make phylo` | run the phylogenetic trait imputation model |
| `make clean` | remove everything in `results/` and `figures/` |

## Dependencies

``` bash
make install
```

This runs `R/install.R` which installs all required R packages (`RTMB`, `unmarked`, `R2jags`, `fmesher`, `ggplot2`, `patchwork`, `ape`, `ggforce`, `ggnewscale`, `viridis`, and the Bioconductor package `ggtree`). JAGS also requires a system-level installation before `R2jags` will work:

- **Debian/Ubuntu**: `sudo apt-get install jags`
- **macOS**: `brew install jags`
- **Windows/other**: <https://mcmc-jags.sourceforge.io>

## Planned Extensions

- **In progress**: spatial GMRF on initial abundance via SPDE (`models/dail_madsen_spde.R`)
  - Possibly running into computational challenges with this model???
- Unified timing comparison across RTMB, unmarked, JAGS.

## References

- Dail, D. and Madsen, L. (2011) Models for Estimating Abundance from Repeated Counts of an Open Metapopulation. *Biometrics* 67:577-587.
- MacKenzie, D.I. et al. (2002) Estimating Site Occupancy Rates When Detection Probabilities Are Less Than One. *Ecology* 83:2248-2255.
- Royle, J.A. (2004) N-Mixture Models for Estimating Population Size from Spatially Replicated Counts. *Biometrics* 60:108-115.
- Gelman, A. and Rubin, D.B. (1992) Inference from Iterative Simulation Using Multiple Sequences. *Statistical Science* 7:457-472.
