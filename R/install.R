# install.R
# Installs all R packages required to run the RTMBdre repository.
# Usage: Rscript R/install.R
#   or:  make install
#
# Note: JAGS also requires a system-level installation.
#   Debian/Ubuntu: sudo apt-get install jags
#   macOS:         brew install jags
#   Windows/other: https://mcmc-jags.sourceforge.io

pkgs <- c(
  "RTMB", # automatic differentiation + sequential reduction (SR)
  "unmarked", # MLE benchmark for occupancy, N-mixture, open N-mixture
  "R2jags", # JAGS interface (requires system JAGS)
  "fmesher", # SPDE mesh construction for the spatial open N-mixture model
  "sf", # spatial field plotting for the SPDE model
  "scales", # log-scale axis labels on the timing plot
  "ggplot2", # timing violin plot and SPDE field plot
  "patchwork", # combining ggplot2 panels
  "ape", # phylogenetic tree I/O for trait imputation
  "ggforce", # trait imputation figure
  "ggnewscale", # trait imputation figure
  "viridis" # trait imputation figure
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(to_install) == 0) {
  cat("All required CRAN packages are already installed.\n")
} else {
  cat("Installing:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, repos = "https://cloud.r-project.org")
  cat("Done.\n")
}

# ggtree is Bioconductor-only, used by the trait imputation figure
if (!requireNamespace("ggtree", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("ggtree", update = FALSE, ask = FALSE)
}
