# Build dendrograms for genotype, phenotype and microbiome profiles

setwd("~/MOrDOr_project")
require(readr)
require(data.table)
require(pvclust)
source("MOrDOr_functions.R")

### Genotype ----

#import dataset
genotyping <- read_csv("Genotype/data_Genotyping_mordor.csv")

genotyping = genotyping[3:28,]
genotyping[,2] <- NULL
genotyping <- transpose(genotyping)
colnames(genotyping) <- genotyping[1,]
genotyping = genotyping[2:41,]

#cluster
fit_geno <- pvclust(genotyping, method.hclust="average",
                    method.dist="euclidean",nboot=100)
plot(fit_geno, hang=-1) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit_geno, alpha=.95)

bootstraps <- (round(fit_geno$edges,2)*100)[,1:2]
my_tree <- as.phylo.hclust.with.nodenames(fit_geno$hclust, nodenames=bootstraps[,2])
write.tree(phy=my_tree, file="tree_genotyping_bt.newick",tree.names=TRUE,digits=2)
plot(my_tree,show.node.label=TRUE)

### Phenotype ----
