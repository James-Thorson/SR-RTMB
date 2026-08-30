.PHONY: all clean install occupancy dynamic_occupancy nmixture open_nmixture spde phylo plots

# set global variables that can be passed into all simulation scripts
NSIM ?= 100
SEED ?= 11223344
export NSIM SEED

all: plots spde phylo

# make sure packages are installed
install:
	Rscript R/install.R

# run the 2x2 occupancy style model simulations
occupancy: install
	Rscript models/occupancy.R

dynamic_occupancy: install
	Rscript models/dynamic_occupancy.R

nmixture: install
	Rscript models/nmixture.R

open_nmixture: install
	Rscript models/open_nmixture.R

# run the open N-mixture SPDE simulation
spde: install
	Rscript models/open_nmixture_spde.R

# run the phylogenetic mixed trait analysis
phylo: install
	Rscript models/phylogenetic_mixed_traits.R

# make simulation plots
plots: occupancy dynamic_occupancy nmixture open_nmixture
	Rscript R/plots.R

clean:
	# wipe results and plots
	rm -f results/*.rds results/timings.csv figures/*.png
