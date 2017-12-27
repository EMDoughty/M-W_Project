library(plyr)

setwd("~/Dropbox/M&W")

load("~/Dropbox/M&W/trees_data.RData")

source("~/Dropbox/M&W/functions.R")

# Low epsilon
MS.tab.full.0.0 <- lapply(names(tree_list_full), MS.estimates.full, epsilon = 0)
MS.tab.full.0.0 <- ldply(MS.tab.full.0.0, data.frame)
row.names(MS.tab.full.0.0) <- MS.tab.full.0.0$Tree.clade

MS.tab.half.0.0 <- lapply(names(tree_list_half), MS.estimates.half, epsilon = 0)
MS.tab.half.0.0 <- ldply(MS.tab.half.0.0, data.frame)
row.names(MS.tab.half.0.0) <- MS.tab.half.0.0$Tree.clade

MS.tab.quarter.0.0 <- lapply(names(tree_list_quarter), MS.estimates.quarter, epsilon = 0)
MS.tab.quarter.0.0 <- ldply(MS.tab.quarter.0.0, data.frame)
row.names(MS.tab.quarter.0.0) <- MS.tab.quarter.0.0$Tree.clade


# Medium epsilon
MS.tab.full.0.5 <- lapply(names(tree_list_full), MS.estimates.full, epsilon = 0.5)
MS.tab.full.0.5 <- ldply(MS.tab.full.0.5, data.frame)
row.names(MS.tab.full.0.5) <- MS.tab.full.0.5$Tree.clade

MS.tab.half.0.5 <- lapply(names(tree_list_half), MS.estimates.half, epsilon = 0.5)
MS.tab.half.0.5 <- ldply(MS.tab.half.0.5, data.frame)
row.names(MS.tab.half.0.5) <- MS.tab.half.0.5$Tree.clade

MS.tab.quarter.0.5 <- lapply(names(tree_list_quarter), MS.estimates.quarter, epsilon = 0.5)
MS.tab.quarter.0.5 <- ldply(MS.tab.quarter.0.5, data.frame)
row.names(MS.tab.quarter.0.5) <- MS.tab.quarter.0.5$Tree.clade


# High epsilon
MS.tab.full.0.9 <- lapply(names(tree_list_full), MS.estimates.full, epsilon = 0.9)
MS.tab.full.0.9 <- ldply(MS.tab.full.0.9, data.frame)
row.names(MS.tab.full.0.9) <- MS.tab.full.0.9$Tree.clade

MS.tab.half.0.9 <- lapply(names(tree_list_half), MS.estimates.half, epsilon = 0.9)
MS.tab.half.0.9 <- ldply(MS.tab.half.0.9, data.frame)
row.names(MS.tab.half.0.9) <- MS.tab.half.0.9$Tree.clade

MS.tab.quarter.0.9 <- lapply(names(tree_list_quarter), MS.estimates.quarter, epsilon = 0.9)
MS.tab.quarter.0.9 <- ldply(MS.tab.quarter.0.9, data.frame)
row.names(MS.tab.quarter.0.9) <- MS.tab.quarter.0.9$Tree.clade


# Null (true) distributions of r
null.tab <- read.csv(file = "evo13378-sup-0006-TablesS1_edited.csv")
null.r <- null.tab$Net.Div.
names(null.r) <- paste(null.tab$Backbone.tree, null.tab$Clade, sep = ".")

######################################################################################################
# Work in progress
# Begin
######################################################################################################

A.rates <- null.r[1:10]
B.rates <- null.r[11:20]
C.rates <- null.r[21:30]
D.rates <- null.r[31:40]
E.rates <- null.r[41:50]
F.rates <- null.r[51:60]
G.rates <- null.r[61:70]
H.rates <- null.r[71:80]
I.rates <- null.r[81:90]
J.rates <- null.r[91:100]
K.rates <- null.r[101:110]
L.rates <- null.r[111:120]
M.rates <- null.r[121:130]
N.rates <- null.r[131:140]
O.rates <- null.r[141:150]
P.rates <- null.r[151:160]
Q.rates <- null.r[161:170]
R.rates <- null.r[171:180]
S.rates <- null.r[181:190]
T.rates <- null.r[191:200]


T.null.dist <- sample(T.rates, size = 100, replace = T)
A.null.r.lower.tail <- quantile(sample(A.rates, size = 10000, replace = T), probs = 0.025)
A.null.r.upper.tail <- quantile(sample(A.rates, size = 10000, replace = T), probs = 0.975)


png("Null.distribution.r.png")
plot(density(null.r$null.r), main = "Null distribution of r", xlab = "Net diversification rate (r)", ylab = "Probability density", font.main = 1) # Looks pretty normal
abline(v = c(null.r.lower.tail, null.r.upper.tail), lty = 2)
dev.off()


######################################################################################################
# End
# Work in progress
######################################################################################################

null.r.lower.tail <- quantile(null.r, probs = 0.025)
null.r.upper.tail <- quantile(null.r, probs = 0.975)

Truly.exceptional <- null.r <= null.r.lower.tail | null.r >= null.r.upper.tail # Determine if null (true) rates are in tails of null distribution

null.r <- data.frame(null.r, Truly.exceptional)

png("Null.distribution.r.png")
plot(density(null.r$null.r), main = "Null distribution of r", xlab = "Net diversification rate (r)", ylab = "Probability density", font.main = 1) # Looks pretty normal
abline(v = c(null.r.lower.tail, null.r.upper.tail), lty = 2)
dev.off()

# Summarize results by sampling level
MS.summary.full <- null.r
MS.summary.full <- cbind(MS.summary.full, 
                         MS.tab.full.0.0$exceptionally.depauperate.crown, MS.tab.full.0.0$exceptionally.diverse.crown, MS.tab.full.0.0$exceptional.crown.pval, MS.tab.full.0.0$exceptionally.depauperate.stem, MS.tab.full.0.0$exceptionally.diverse.stem, MS.tab.full.0.0$exceptional.stem.pval, 
                         MS.tab.full.0.5$exceptionally.depauperate.crown, MS.tab.full.0.5$exceptionally.diverse.crown, MS.tab.full.0.5$exceptional.crown.pval, MS.tab.full.0.5$exceptionally.depauperate.stem, MS.tab.full.0.5$exceptionally.diverse.stem, MS.tab.full.0.5$exceptional.stem.pval,
                         MS.tab.full.0.9$exceptionally.depauperate.crown, MS.tab.full.0.9$exceptionally.diverse.crown, MS.tab.full.0.9$exceptional.crown.pval, MS.tab.full.0.9$exceptionally.depauperate.stem, MS.tab.full.0.9$exceptionally.diverse.stem, MS.tab.full.0.9$exceptional.stem.pval)
colnames(MS.summary.full) <- c("Null.r", "Truly.exceptional", 
                               "Except.depaup.crown.0.0", "Except.diverse.crown.0.0", "Pval.crown.0.0", "Except.depaup.stem.0.0", "Except.diverse.stem.0.0", "Pval.stem.0.0", 
                               "Except.depaup.crown.0.5", "Except.diverse.crown.0.5", "Pval.crown.0.5", "Except.depaup.stem.0.5", "Except.diverse.stem.0.5", "Pval.stem.0.5", 
                               "Except.depaup.crown.0.9", "Except.diverse.crown.0.9", "Pval.crown.0.9", "Except.depaup.stem.0.9", "Except.diverse.stem.0.9", "Pval.stem.0.9")

MS.summary.half <- null.r
MS.summary.half <- cbind(MS.summary.half, 
                         MS.tab.half.0.0$exceptionally.depauperate.crown, MS.tab.half.0.0$exceptionally.diverse.crown, MS.tab.half.0.0$exceptional.crown.pval, MS.tab.half.0.0$exceptionally.depauperate.stem, MS.tab.half.0.0$exceptionally.diverse.stem, MS.tab.half.0.0$exceptional.stem.pval, 
                         MS.tab.half.0.5$exceptionally.depauperate.crown, MS.tab.half.0.5$exceptionally.diverse.crown, MS.tab.half.0.5$exceptional.crown.pval, MS.tab.half.0.5$exceptionally.depauperate.stem, MS.tab.half.0.5$exceptionally.diverse.stem, MS.tab.half.0.5$exceptional.stem.pval,
                         MS.tab.half.0.9$exceptionally.depauperate.crown, MS.tab.half.0.9$exceptionally.diverse.crown, MS.tab.half.0.9$exceptional.crown.pval, MS.tab.half.0.9$exceptionally.depauperate.stem, MS.tab.half.0.9$exceptionally.diverse.stem, MS.tab.half.0.9$exceptional.stem.pval)
colnames(MS.summary.half) <- c("Null.r", "Truly.exceptional", 
                               "Except.depaup.crown.0.0", "Except.diverse.crown.0.0", "Pval.crown.0.0", "Except.depaup.stem.0.0", "Except.diverse.stem.0.0", "Pval.stem.0.0", 
                               "Except.depaup.crown.0.5", "Except.diverse.crown.0.5", "Pval.crown.0.5", "Except.depaup.stem.0.5", "Except.diverse.stem.0.5", "Pval.stem.0.5", 
                               "Except.depaup.crown.0.9", "Except.diverse.crown.0.9", "Pval.crown.0.9", "Except.depaup.stem.0.9", "Except.diverse.stem.0.9", "Pval.stem.0.9")

MS.summary.quarter <- null.r
MS.summary.quarter <- cbind(MS.summary.quarter, 
                         MS.tab.quarter.0.0$exceptionally.depauperate.crown, MS.tab.quarter.0.0$exceptionally.diverse.crown, MS.tab.quarter.0.0$exceptional.crown.pval, MS.tab.quarter.0.0$exceptionally.depauperate.stem, MS.tab.quarter.0.0$exceptionally.diverse.stem, MS.tab.quarter.0.0$exceptional.stem.pval, 
                         MS.tab.quarter.0.5$exceptionally.depauperate.crown, MS.tab.quarter.0.5$exceptionally.diverse.crown, MS.tab.quarter.0.5$exceptional.crown.pval, MS.tab.quarter.0.5$exceptionally.depauperate.stem, MS.tab.quarter.0.5$exceptionally.diverse.stem, MS.tab.quarter.0.5$exceptional.stem.pval,
                         MS.tab.quarter.0.9$exceptionally.depauperate.crown, MS.tab.quarter.0.9$exceptionally.diverse.crown, MS.tab.quarter.0.9$exceptional.crown.pval, MS.tab.quarter.0.9$exceptionally.depauperate.stem, MS.tab.quarter.0.9$exceptionally.diverse.stem, MS.tab.quarter.0.9$exceptional.stem.pval)
colnames(MS.summary.quarter) <- c("Null.r", "Truly.exceptional", 
                               "Except.depaup.crown.0.0", "Except.diverse.crown.0.0", "Pval.crown.0.0", "Except.depaup.stem.0.0", "Except.diverse.stem.0.0", "Pval.stem.0.0", 
                               "Except.depaup.crown.0.5", "Except.diverse.crown.0.5", "Pval.crown.0.5", "Except.depaup.stem.0.5", "Except.diverse.stem.0.5", "Pval.stem.0.5", 
                               "Except.depaup.crown.0.9", "Except.diverse.crown.0.9", "Pval.crown.0.9", "Except.depaup.stem.0.9", "Except.diverse.stem.0.9", "Pval.stem.0.9")

# Save summary tables as .csv files
write.csv(MS.summary.full, file = "MS.summary.full.csv")
write.csv(MS.summary.half, file = "MS.summary.half.csv")
write.csv(MS.summary.quarter, file = "MS.summary.quarter.csv")


# Summarize results by correct and incorrect (i.e., Type I and Type II error) MS inference



