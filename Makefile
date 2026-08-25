.PHONY: all run parallel plots occupancy nmixture dailmadsen spde phylo install clean

all: run

# run all models sequentially
test:
	Rscript R/run_all.R

# run all models in parallel (one sim per core)
run:
	Rscript R/run_all.R parallel

# generate figures from results/*.rds (run after run/test)
plots:
	Rscript R/plots.R

# run individual models
occupancy:
	Rscript models/occupancy.R

nmixture:
	Rscript models/nmixture.R

dailmadsen:
	Rscript models/dail_madsen.R

spde:
	Rscript models/dail_madsen_spde.R

phylo:
	Rscript models/phylogenetic_mixed_traits.R

# install required R packages
install:
	Rscript R/install.R

# clean results
clean:
	rm -f results/*.rds results/*.csv figures/*.png
