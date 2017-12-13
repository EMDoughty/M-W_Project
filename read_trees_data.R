library(geiger)
library(phytools)

setwd("~/Dropbox/M&W/")


# Import data
# Import table of crown and stem ages
trees.dat <- read.csv(file = "evo13378-sup-0006-TablesS1_edited.csv")
trees.dat <- trees.dat[, 1:5]

# Import trees
trees.full <- read.tree(file = "evo13378-sup-0012-SuppMat.tre")
trees.half <- read.tree(file = "evo13378-sup-0013-SuppMat.tre")
trees.quarter <- read.tree(file = "evo13378-sup-0014-SuppMat.tre")


# Plot trees
png(file = "trees.full.png")
par(mar = c(1, 1, 1, 1) + 0.1, mfrow = c(4, 5))
for(i in 1:length(trees.full)){
  plot(ladderize(trees.full[[i]]), show.tip.label = F, main = names(trees.full[i]), font.main = 1)
  rm(i)
}
dev.off()

png(file = "trees.half.png")
par(mar = c(1, 1, 1, 1) + 0.1, mfrow = c(4, 5))
for(i in 1:length(trees.half)){
  plot(ladderize(trees.half[[i]]), show.tip.label = F, main = names(trees.half[i]), font.main = 1)
  rm(i)
}
dev.off()

png(file = "trees.quarter.png")
par(mar = c(1, 1, 1, 1) + 0.1, mfrow = c(4, 5))
for(i in 1:length(trees.quarter)){
  plot(ladderize(trees.quarter[[i]]), show.tip.label = F, main = names(trees.quarter[i]), font.main = 1)
  rm(i)
}
dev.off()

# Munge data
# Make list of trees.dat
tree_list_full <- split(trees.dat, list(trees.dat$Backbone.tree, trees.dat$Clade)) # 200 elements: A.a, B.a, ...

# Make lists for incompletely sampled trees
trees.labs <- trees.dat
trees.labs[, 3:5] <- rep(NA, nrow(trees.labs))
tree_list_half <- split(trees.labs, list(trees.labs$Backbone.tree, trees.labs$Clade))
tree_list_quarter <- split(trees.labs, list(trees.labs$Backbone.tree, trees.labs$Clade))

# Add completely sampled trees and clades to list
for(i in 1:length(tree_list_full)){
  tree_list_full[[i]] <- append(tree_list_full[[i]], trees.full[paste("tree_", casefold(tree_list_full[[i]][, 1]), sep = "")])
  clade <- extract.clade(tree_list_full[[i]][[6]], getMRCA(tree_list_full[[i]][[6]], tree_list_full[[i]][[6]][[3]][grep(strsplit(x = names(tree_list_full[i]), split = "\\.")[[1]][2], tree_list_full[[i]][[6]][[3]])]))
  tree_list_full[[i]] <- append(tree_list_full[[i]], list(clade))
  names(tree_list_full[[i]])[7] <- paste("clade_", tree_list_full[[i]][[2]], sep = "")
  rm(i, clade)
}

# Add incompletely sampled trees, clades, tip numbers, crown and stem ages to lists
for(i in 1:length(tree_list_half)){
  tree_list_half[[i]] <- append(tree_list_half[[i]], trees.half[paste("tree_", casefold(tree_list_half[[i]][, 1]), sep = "")])
  clade <- extract.clade(tree_list_half[[i]][[6]], getMRCA(tree_list_half[[i]][[6]], tree_list_half[[i]][[6]][[3]][grep(strsplit(x = names(tree_list_half[i]), split = "\\.")[[1]][2], tree_list_half[[i]][[6]][[3]])]))
  tree_list_half[[i]] <- append(tree_list_half[[i]], list(clade))
  names(tree_list_half[[i]])[7] <- paste("clade_", tree_list_half[[i]][[2]], sep = "")
  ntips <- length(tree_list_half[[i]][[7]][[3]])
  tree_list_half[[i]][[3]] <- ntips
  H <- nodeHeights(tree_list_half[[i]][[6]])
  nn <- findMRCA(tree_list_half[[i]][[6]], tree_list_half[[i]][[7]][[3]])
  stem.age <- 100 - H[tree_list_half[[i]][[6]][[1]] == phytools:::getAncestors(tree_list_half[[i]][[6]], node = nn, type = "parent")][1]
  tree_list_half[[i]][[4]] <- round(stem.age, digits = 3)
  crown.age <- 100 - nodeheight(tree = tree_list_half[[i]][[6]], node = nn)
  tree_list_half[[i]][[5]] <- round(crown.age, digits = 3)
  rm(i, clade, ntips, H, nn, stem.age, crown.age)
}

for(i in 1:length(tree_list_quarter)){
  tree_list_quarter[[i]] <- append(tree_list_quarter[[i]], trees.quarter[paste("tree_", casefold(tree_list_quarter[[i]][, 1]), sep = "")])
  clade <- extract.clade(tree_list_quarter[[i]][[6]], getMRCA(tree_list_quarter[[i]][[6]], tree_list_quarter[[i]][[6]][[3]][grep(strsplit(x = names(tree_list_quarter[i]), split = "\\.")[[1]][2], tree_list_quarter[[i]][[6]][[3]])]))
  tree_list_quarter[[i]] <- append(tree_list_quarter[[i]], list(clade))
  names(tree_list_quarter[[i]])[7] <- paste("clade_", tree_list_quarter[[i]][[2]], sep = "")
  ntips <- length(tree_list_quarter[[i]][[7]][[3]])
  tree_list_quarter[[i]][[3]] <- ntips
  H <- nodeHeights(tree_list_quarter[[i]][[6]])
  nn <- findMRCA(tree_list_quarter[[i]][[6]], tree_list_quarter[[i]][[7]][[3]])
  stem.age <- 100 - H[tree_list_quarter[[i]][[6]][[1]] == phytools:::getAncestors(tree_list_quarter[[i]][[6]], node = nn, type = "parent")][1]
  tree_list_quarter[[i]][[4]] <- round(stem.age, digits = 3)
  crown.age <- 100 - nodeheight(tree = tree_list_quarter[[i]][[6]], node = nn)
  tree_list_quarter[[i]][[5]] <- round(crown.age, digits = 3)
  rm(i, clade, ntips, H, nn, stem.age, crown.age)
}


save.image("~/Dropbox/M&W/trees_data.RData")
