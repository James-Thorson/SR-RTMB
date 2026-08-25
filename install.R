# install.R
# Installs all R packages required to run the RTMBdre repository.
# Usage: Rscript install.R
#   or:  make install
#
# Note: JAGS also requires a system-level installation.
#   Debian/Ubuntu: sudo apt-get install jags
#   macOS:         brew install jags
#   Windows/other: https://mcmc-jags.sourceforge.io

pkgs <- c(
  "RTMB", # automatic differentiation + sequential reduction (SR)
  "unmarked", # MLE benchmark for occupancy, N-mixture, Dail-Madsen
  "R2jags", # JAGS interface (requires system JAGS)
  "fmesher", # SPDE mesh construction for spatial Dail-Madsen
  "parallel",
  "ggplot2", # timing violin plot and SPDE field plot
  "patchwork" # combining ggplot2 panels
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(to_install) == 0) {
  cat("All required packages are already installed.\n")
} else {
  cat("Installing:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, repos = "https://cloud.r-project.org")
  cat("Done.\n")
}
