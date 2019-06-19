# Analysis of microsatellites markers used to genotype P. ramosa seeds

setwd("~/MOrDOr_project")
require(readr)
require(poppr)
require(ape)

### Import dataset ----

monpop <- read.genalex("Genotype/data_Genotyping_mordor.csv")
splitStrata(monpop) <- ~field/Year/host
poppr(monpop)

### Calculate MLG ----
P.tab <- mlg.table(monpop)

### Build and extract phylogenetic tree ----

my_tree <- bruvo.boot(monpop, replen = 2, add = TRUE, loss = TRUE, sample = 999,
                       tree = "nj", showtree = TRUE)

write.tree(phy=my_tree, file="tree_genotyping.newick")
