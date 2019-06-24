# Build dendrograms for genotype, phenotype and microbiome profiles

setwd("~/MOrDOr_project")
require(readr)
require(data.table)
require(pvclust)
require(phyloseq)
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
ed_all <- read_csv("Phenotype/EC50.csv")

ed_matrix <- as.data.frame(levels(as.factor(ed_all$sample)))
colnames(ed_matrix) <- "sample"
rownames(ed_matrix)=ed_matrix$sample

for (i in levels(as.factor(ed_all$sample))) {
  ed_matrix[i,"2PEITC"] <- log10(ed_all$EC50[ed_all$molecules == "2PEITC" & ed_all$sample == i])
  ed_matrix[i,"racGR24"] <- log10(ed_all$EC50[ed_all$molecules == "racGR24" & ed_all$sample == i])
  ed_matrix[i,"+GR24"] <- log10(ed_all$EC50[ed_all$molecules == "+GR24" & ed_all$sample == i])
  ed_matrix[i,"-GR24"] <- log10(ed_all$EC50[ed_all$molecules == "-GR24" & ed_all$sample == i])
  ed_matrix[i,"+eGR24"] <- log10(ed_all$EC50[ed_all$molecules == "+eGR24" & ed_all$sample == i])
  ed_matrix[i,"-eGR24"] <- log10(ed_all$EC50[ed_all$molecules == "-eGR24" & ed_all$sample == i])
  }
ed_matrix = ed_matrix[1:(nrow(ed_matrix)-2),2:ncol(ed_matrix)]
ed_matrix = transpose(ed_matrix)
colnames(ed_matrix) <- levels(as.factor(ed_all$sample))[1:26]

fit_ec50 <- pvclust(ed_matrix, method.hclust="average",
               method.dist="euclidean",nboot=100)
plot(fit_ec50, hang=-1) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit_ec50, alpha=.95)

bootstraps <- (round(fit_ec50$edges,2)*100)[,1:2]
my_tree <- as.phylo.hclust.with.nodenames(fit_ec50$hclust, nodenames=bootstraps[,2])
write.tree(phy=my_tree, file="tree_ec50_bt.newick",tree.names=TRUE,digits=2)
plot(my_tree,show.node.label=TRUE)

### Microbiome profile ----
ps = ps_16S_seed
ps_norm <- transform_sample_counts(ps,function(x){1E6*x/sum(x)})
ps_merged <- merge_samples(ps_norm,
                           group = ps_norm@sam_data$plot_sample)

metadata_otu <- transpose(as.data.frame(ps_merged@otu_table))
colnames(metadata_otu) <- sample_names(ps_merged)

fit <- pvclust(metadata_otu, method.hclust="average",
               method.dist="euclidean",nboot=100)
plot(fit, hang=-1) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)

bootstraps <- (round(fit$edges,2)*100)[,1:2]
my_tree <- as.phylo.hclust.with.nodenames(fit$hclust, nodenames=bootstraps[,2])
write.tree(phy=my_tree, file="tree_ITS_seed_bp.newick",tree.names=TRUE,digits=2)
plot(my_tree,show.node.label=TRUE)
