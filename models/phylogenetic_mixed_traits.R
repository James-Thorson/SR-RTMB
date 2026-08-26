
#################
# Copied from: C:\Users\jtuth\Desktop\Work\Collab-2026\2026 -- sequential reduction\state_switching_20226-08-23.R
#################

# Make sure to use github version of TMB and RTMB
# pak::pak("kaskr/adcomp/TMB")
# pak::pak("kaskr/RTMB/RTMB")

data_dir = "data"
results_dir = "figures"

library(RTMB)
library(ape)
library(ggplot2)
library(ggtree)
library(ggforce)
library(ggnewscale)
library(viridis)

tree <- read.tree( file.path(data_dir,"Liolaemus_tree_final_parity.tre") )
vroot = tree$Nnode + 1
edge_ez = tree$edge  # parent, child
length_e = tree$edge.length
p_error = 0.001

viriparity = readRDS( file = file.path(data_dir, "viviparity.RDS") )
v_i = match( viriparity$species, c(tree$tip.label,tree$node.label) )

#
y1_i = viriparity$viviparity
x1_v = rep( round(mean(y1_i,na.rm=TRUE)), max(tree$edge) )
x1_v[v_i] = y1_i

y2_i = x2_i = viriparity$size
x2_v = rep( NA, max(tree$edge) )
x2_v[v_i] = y2_i
x2_v = ifelse( is.na(x2_v), mean(y2_i,na.rm=TRUE), x2_v )

# Define map
map = list(
  logr = factor(c(1,1))
)

###################
# Fit in RTMB
# NEXT ... try adding measurement error to check mean vs. mode
###################

#
p = list(
  logr = log( c(0.1,0.1) ),
  x1_v = x1_v,                 # Use real one, so mapped off when data are present (no assignment error)
  x2_v = x2_v,                 # Use real one, so mapped off when data are present (no assignment error)
  logsigma = log(1),
  logsd = log(1),
  finvrho = 0.2,
  beta = 0
)
map$finvrho = factor(NA)

make_onehot <- function(n) {
  old <- TapeConfig()
  on.exit(TapeConfig(old))
  TapeConfig(comparison="tape")
  MakeTape(function(x) x==(1:n), 1)
}
onehot_comparison = make_onehot(2)
onehot <-
function( level,  # Count from zero
          nlevels = 2,
          type = c("abs", "index", "comparison") ){

  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  type = match.arg(type)
  if(type == "abs"){
    vec = 1 - abs(seq_len(nlevels) - level - 1)
    if(nlevels>2) vec = (vec + abs(vec))/2
  }
  if(type == "index"){
    vec = rep(0, nlevels)
    vec[level + 1] = 1
  }
  if(type == "comparison"){
    vec = onehot_comparison(level + 1)
  }
  return(vec)
}

#
onehot_type = c("abs", "index", "comparison")[1]
get_jnll = function(p, what = "jnll"){
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  # Process model
  Q = cbind( c(-exp(p$logr[1]),exp(p$logr[1])), c(exp(p$logr[2]),-exp(p$logr[2])))
  loglik1_v = matrix( 0, nrow = length(p$x1_v), ncol = 2 )
  for(ei in 1:nrow(edge_ez)){
    vchild = edge_ez[ei,2]
    vparent = edge_ez[ei,1]
    M = Matrix::expm( length_e[ei] * Q )
    xhat_child = onehot(p$x1_v[vparent], type = onehot_type) %*% M
    x_child = onehot(p$x1_v[vchild], type = onehot_type)
    loglik1_v[vchild,1] = dmultinom( x_child, prob = xhat_child, size = 1, log = TRUE )
    loglik1_v[vchild,2] = dnorm( p$x2_v[vchild], mean = p$x2_v[vparent], sd = sqrt(length_e[ei]) * exp(p$logsigma), log = TRUE )    # mean = (2*plogis(p$finvrho)-1) * p$x2_t[ti-1],
  }
  # Measurement model
  loglik2_i = matrix( 0, nrow = length(v_i), ncol = 2 )
  for(ii in seq_along(v_i)){
    if( p_error > 0 ){
      loglik2_i[ii,1] = dbinom( y1_i[ii], prob = p_error + (1-2*p_error) * p$x1_v[v_i[ii]], size = 1, log = TRUE )
    }
    loglik2_i[ii,2] = dnorm( y2_i[ii], mean = p$x2_v[v_i[ii]] + p$beta * p$x1_v[v_i[ii]], sd = exp(p$logsd), log = TRUE )
  }
  jll = sum(loglik1_v, na.rm=TRUE) + sum(loglik2_i, na.rm = TRUE)
  REPORT(jll)
  REPORT(Q)
  if( what == "jnll" ){
    out = -1 * jll
  }else{
    out = mget(c("p","loglik1_v","loglik2_i","Q","jll"))
  }
  return(out)
}
get_jnll(p, what = "all")

obj = MakeADFun(
  get_jnll,
  parameters = p,
  integrate = list( x1_v = TMB::SR(c(0,1), discrete=TRUE), x2_v = TMB:::LA() ),
  random = c("x1_v","x2_v"),
  map = map
)
obj$fn(obj$par)

opt = nlminb(
  obj$par, obj$fn, obj$gr,
  control = list(trace = 1)
)
obj$gr(opt$par)

#sdr <- sdreport(obj, type="mode") ## gives exact posterior mode
sdr <- sdreport(obj) ## gives exact posterior mean
x1hat_v = as.list(sdr, what="Estimate")$x1_v
x2hat_v = as.list(sdr, what="Estimate")$x2_v
#prob1_t = sdr$value[names(sdr$value) == "prob1_t"]
cbind( x1_v, y1_i[match(seq_along(x1_v),v_i)], x1hat_v )

# Transition matrix
phat = p
phat$logr[] = opt$par['logr']
get_jnll(phat, what = "all")$Q


# Make tree
p <- ggtree(tree, layout = "circular") +
  geom_tiplab2(size=2.5,offset=2)

# Get coordinates of ALL nodes
pd <- p$data

# Add your trait values
pd$value1 <- x1hat_v[pd$node]
pd$value2 <- x2hat_v[pd$node]

# Plot heatmap cells for tips + internal nodes
#p1 <- p +
#  geom_tile(
#    data = pd,
#    aes(
#      x = x,
#      y = y,
#      fill = value1
#    ),
#    width = 1,
#    height = 0.9
#  ) +
#  scale_fill_viridis_c(
#    option = "D",
#    name = "value",
#    na.value = "white"
#  )
#
#p1


p1 <- p +

  # -------------------------
  # Trait 1: left semicircle
  # -------------------------
  geom_arc_bar(
    data = pd,
    mapping = aes(
      x0 = x,
      y0 = y,
      r0 = 0,
      r = 0.5,
      start = pi / 2,
      end = 3 * pi / 2,
      fill = value1
    ),
    inherit.aes = FALSE,
    color = NA
  ) +

  scale_fill_viridis_c(
    option = "D",
    name = "Pr(viviparous)",
    na.value = "white"
  ) +

  ggnewscale::new_scale_fill() +

  # -------------------------
  # Trait 2: right semicircle
  # -------------------------
  geom_arc_bar(
    data = pd,
    mapping = aes(
      x0 = x,
      y0 = y,
      r0 = 0,
      r = 0.5,
      start = -pi / 2,
      end = pi / 2,
      fill = value2
    ),
    inherit.aes = FALSE,
    color = NA
  ) +

  scale_fill_viridis_c(
    option = "C",
    name = "Body size (log10)",
    na.value = "white"
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

ggsave(
  plot = p1,
  filename= file.path(results_dir, "trait_imputation.png"),
  height = 8, width=6
)

