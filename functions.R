MS.estimates.full <- function(Tree.clade, epsilon){ # e.g., MS.estimates.full(Tree.clade = "A.a", epsilon = 0.9)
  require(geiger)
  clade.n.tips <- tree_list_full[[Tree.clade]][[3]]
  Tree.r <- bd.ms(phy = tree_list_full[[Tree.clade]][[6]], crown = TRUE, epsilon = epsilon)
  clade.crown.r <- bd.ms(time = tree_list_full[[Tree.clade]][[5]], n = clade.n.tips, crown = TRUE, epsilon = epsilon)
  clade.stem.r <- bd.ms(time = tree_list_full[[Tree.clade]][[4]], n = clade.n.tips, crown = FALSE, epsilon = epsilon)
  clade.crown.bounds <- crown.limits(time = tree_list_full[[Tree.clade]][[5]], r = Tree.r, epsilon = epsilon)
  clade.stem.bounds <- stem.limits(time = tree_list_full[[Tree.clade]][[4]], r = Tree.r, epsilon = epsilon)
  clade.crown.prob <- crown.p(time = tree_list_full[[Tree.clade]][[5]], n = clade.n.tips, r = Tree.r, epsilon = epsilon)
  clade.stem.prob <- stem.p(time = tree_list_full[[Tree.clade]][[4]], n = clade.n.tips, r = Tree.r, epsilon = epsilon)
  exceptionally.diverse.crown <- clade.n.tips >= clade.crown.bounds[2]
  exceptionally.depauperate.crown <- clade.n.tips <= clade.crown.bounds[1]
  exceptionally.diverse.stem <- clade.n.tips >= clade.stem.bounds[2]
  exceptionally.depauperate.stem <- clade.n.tips <= clade.stem.bounds[1]
  exceptional.crown.pval <- min(clade.crown.prob, 1 - clade.crown.prob) * 2
  exceptional.stem.pval <- min(clade.stem.prob, 1 - clade.stem.prob) * 2
  estimates <- list(Tree.clade, epsilon, Tree.r, clade.crown.r, clade.stem.r, clade.n.tips, clade.crown.bounds, exceptionally.diverse.crown, exceptionally.depauperate.crown, exceptional.crown.pval, clade.stem.bounds, exceptionally.diverse.stem, exceptionally.depauperate.stem, exceptional.stem.pval)
  names(estimates) <- c("Tree.clade", "epsilon", "Tree.r", "clade.crown.r", "clade.stem.r", "clade.n.tips", "clade.crown.bounds", "exceptionally.diverse.crown", "exceptionally.depauperate.crown", "exceptional.crown.pval", "clade.stem.bounds", "exceptionally.diverse.stem", "exceptionally.depauperate.stem", "exceptional.stem.pval")
  return(estimates)
}

MS.estimates.half <- function(Tree.clade, epsilon){ # e.g., MS.estimates.half(Tree.clade = "G.d", epsilon = 0)
  require(geiger)
  full.n.tips <- tree_list_full[[Tree.clade]][[3]]
  clade.n.tips <- tree_list_half[[Tree.clade]][[3]]
  Tree.r <- bd.ms(phy = tree_list_half[[Tree.clade]][[6]], missing = length(tree_list_full[[Tree.clade]][[6]][[3]]) - length(tree_list_half[[Tree.clade]][[6]][[3]]), crown = TRUE, epsilon = epsilon)
  clade.crown.r <- bd.ms(time = tree_list_half[[Tree.clade]][[5]], n = clade.n.tips, missing = length(tree_list_full[[Tree.clade]][[7]][[3]]) - length(tree_list_half[[Tree.clade]][[7]][[3]]), crown = TRUE, epsilon = epsilon)
  clade.stem.r <- bd.ms(time = tree_list_half[[Tree.clade]][[4]], n = clade.n.tips, missing = length(tree_list_full[[Tree.clade]][[7]][[3]]) - length(tree_list_half[[Tree.clade]][[7]][[3]]), crown = FALSE, epsilon = epsilon)
  clade.crown.bounds <- crown.limits(time = tree_list_half[[Tree.clade]][[5]], r = Tree.r, epsilon = epsilon)
  clade.stem.bounds <- stem.limits(time = tree_list_half[[Tree.clade]][[4]], r = Tree.r, epsilon = epsilon)
  clade.crown.prob <- crown.p(time = tree_list_half[[Tree.clade]][[5]], n = clade.n.tips, r = Tree.r, epsilon = epsilon)
  clade.stem.prob <- stem.p(time = tree_list_half[[Tree.clade]][[4]], n = clade.n.tips, r = Tree.r, epsilon = epsilon)
  exceptionally.diverse.crown <- clade.n.tips >= clade.crown.bounds[2]
  exceptionally.depauperate.crown <- clade.n.tips <= clade.crown.bounds[1]
  exceptionally.diverse.stem <- clade.n.tips >= clade.stem.bounds[2]
  exceptionally.depauperate.stem <- clade.n.tips <= clade.stem.bounds[1]
  exceptional.crown.pval <- min(clade.crown.prob, 1 - clade.crown.prob) * 2
  exceptional.stem.pval <- min(clade.stem.prob, 1 - clade.stem.prob) * 2
  estimates <- list(Tree.clade, epsilon, Tree.r, clade.crown.r, clade.stem.r, clade.n.tips, clade.crown.bounds, exceptionally.diverse.crown, exceptionally.depauperate.crown, exceptional.crown.pval, clade.stem.bounds, exceptionally.diverse.stem, exceptionally.depauperate.stem, exceptional.stem.pval)
  names(estimates) <- c("Tree.clade", "epsilon", "Tree.r", "clade.crown.r", "clade.stem.r", "clade.n.tips", "clade.crown.bounds", "exceptionally.diverse.crown", "exceptionally.depauperate.crown", "exceptional.crown.pval", "clade.stem.bounds", "exceptionally.diverse.stem", "exceptionally.depauperate.stem", "exceptional.stem.pval")
  return(estimates)
}

MS.estimates.quarter <- function(Tree.clade, epsilon){ # e.g., MS.estimates.quarter(Tree.clade = "M.j", epsilon = 0.9)
  require(geiger)
  full.n.tips <- tree_list_full[[Tree.clade]][[3]]
  clade.n.tips <- tree_list_quarter[[Tree.clade]][[3]]
  Tree.r <- bd.ms(phy = tree_list_quarter[[Tree.clade]][[6]], missing = length(tree_list_full[[Tree.clade]][[6]][[3]]) - length(tree_list_quarter[[Tree.clade]][[6]][[3]]), crown = TRUE, epsilon = epsilon)
  clade.crown.r <- bd.ms(time = tree_list_quarter[[Tree.clade]][[5]], n = clade.n.tips, missing = length(tree_list_full[[Tree.clade]][[7]][[3]]) - length(tree_list_quarter[[Tree.clade]][[7]][[3]]), crown = TRUE, epsilon = epsilon)
  clade.stem.r <- bd.ms(time = tree_list_quarter[[Tree.clade]][[4]], n = clade.n.tips, missing = length(tree_list_full[[Tree.clade]][[7]][[3]]) - length(tree_list_quarter[[Tree.clade]][[7]][[3]]), crown = FALSE, epsilon = epsilon)
  clade.crown.bounds <- crown.limits(time = tree_list_quarter[[Tree.clade]][[5]], r = Tree.r, epsilon = epsilon)
  clade.stem.bounds <- stem.limits(time = tree_list_quarter[[Tree.clade]][[4]], r = Tree.r, epsilon = epsilon)
  clade.crown.prob <- crown.p(time = tree_list_quarter[[Tree.clade]][[5]], n = clade.n.tips, r = Tree.r, epsilon = epsilon)
  clade.stem.prob <- stem.p(time = tree_list_quarter[[Tree.clade]][[4]], n = clade.n.tips, r = Tree.r, epsilon = epsilon)
  exceptionally.diverse.crown <- clade.n.tips >= clade.crown.bounds[2]
  exceptionally.depauperate.crown <- clade.n.tips <= clade.crown.bounds[1]
  exceptionally.diverse.stem <- clade.n.tips >= clade.stem.bounds[2]
  exceptionally.depauperate.stem <- clade.n.tips <= clade.stem.bounds[1]
  exceptional.crown.pval <- min(clade.crown.prob, 1 - clade.crown.prob) * 2
  exceptional.stem.pval <- min(clade.stem.prob, 1 - clade.stem.prob) * 2
  estimates <- list(Tree.clade, epsilon, Tree.r, clade.crown.r, clade.stem.r, clade.n.tips, clade.crown.bounds, exceptionally.diverse.crown, exceptionally.depauperate.crown, exceptional.crown.pval, clade.stem.bounds, exceptionally.diverse.stem, exceptionally.depauperate.stem, exceptional.stem.pval)
  names(estimates) <- c("Tree.clade", "epsilon", "Tree.r", "clade.crown.r", "clade.stem.r", "clade.n.tips", "clade.crown.bounds", "exceptionally.diverse.crown", "exceptionally.depauperate.crown", "exceptional.crown.pval", "clade.stem.bounds", "exceptionally.diverse.stem", "exceptionally.depauperate.stem", "exceptional.stem.pval")
  return(estimates)
}
