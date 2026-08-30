# Shared helpers sourced by run_all.R and the model files.

# Run FUN over X sequentially (mc.cores <= 1) or in parallel via
# parallel::mclapply (mc.cores > 1). Centralizing this means every model's
# run_*() function times each replicate individually regardless of mode,
# so per-sim timings are always available and mean timings are always a
# true per-fit time rather than wall-clock-divided-by-nsim.
lapply_maybe <- function(X, FUN, mc.cores = 1, ...) {
  if (mc.cores > 1) {
    parallel::mclapply(X, FUN, mc.cores = mc.cores, ...)
  } else {
    lapply(X, FUN, ...)
  }
}
