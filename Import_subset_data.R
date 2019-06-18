# Import and subset data for microbiata analysis

setwd("~/MOrDOr_project")

### Packages ----
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
install.packages("readr")
install.packages("ggplot2")
install.packages("devtools")
install.packages("decontam-master/", repos = NULL, type = "source",
                 dependencies = c("Depends", "Suggests","Imports"))
library(devtools)
devtools::install_github("benjjneb/decontam")

library(phyloseq)
library(readr)
library(ggplot2)
library(decontam)

### Import data ----
abund_16S <- read.csv2("abundance_table_16S.csv", sep = ',',row.names = 1)
taxo_16S <- read.csv2("taxo_16S_silva.csv", sep = ',',row.names = 1, colClasses = "character")
taxo_16S <- as.matrix(taxo_16S)
design <- read_delim("design_16S.csv",",", escape_double = F, trim_ws = T)
rownames(design) <- design$`#SampleID`
tree = read.newick("Qiime2_output/16S/phylogenetic_tree/tree.nwk")

OTU = otu_table(abund_16S, taxa_are_rows = TRUE)
TAX = tax_table(taxo_16S)
DES = sample_data(design)
TREE = phy_tree(tree)

ps_16S = phyloseq(OTU,TAX,DES,TREE)

### Extract ASV id ----
# extract seq ID
ASVid1 = taxa_names(ps_16S)
#set new names
taxa_names(ps_16S) <- paste0("ASV",1:ntaxa(ps_16S),
                             "_",ps_16S@tax_table[,"Genus"],
                             "_",ps_16S@tax_table[,"Species"])
# extract new seq id
ASVid2 = taxa_names(ps_16S)
ASVid = cbind(ASVid1,ASVid2)
write_csv(as.data.frame(ASVid),"ASVid_16S.csv")

### Clean dataset ----
#remove mock and desinfected Pram 120
ps_16S <- prune_samples(!(sample_data(ps_16S)$plot_sample %in% c("Mock",
                                                                 "Pram_120",
                                                                 "Pram_120_d")),
                        ps_16S)
#remove chloroplasts
ps_16S <- subset_taxa(ps_16S,Order != "Chloroplast" )
ps_16S <- prune_taxa(taxa_sums(ps_16S) > 0, ps_16S)

#remove contaminants
#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html#identifying-contaminants-in-marker-gene-and-metagenomics-data
sample_data(ps_16S)$is.neg <- sample_data(ps_16S)$plot_sample %in% c("Control","H2O","PBS")
contamdf.eith <- isContaminant(ps_16S, method="either", 
                               conc="quantif_DNA",
                               neg="is.neg", threshold=0.1)
table(contamdf.eith$contaminant)
head(which(contamdf.eith$contaminant))

ps.pa <- transform_sample_counts(ps_16S, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$plot_sample %in% c("Control","H2O","PBS"), ps.pa)
ps.pa.pos <- prune_samples(!(sample_data(ps.pa)$plot_sample %in% c("Control","H2O","PBS")), ps.pa)

df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.eith$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

plot_frequency(ps_16S,
               taxa_names(ps_16S)[sample(which(contamdf.eith$contaminant),16)],
               conc="quantif_DNA") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

ps_16S_decontam <- prune_taxa(!contamdf.eith$contaminant, ps_16S)
ps_16S_decontam <- prune_samples(!(sample_data(ps_16S_decontam)$plot_sample %in% c("Control","H2O","PBS")),ps_16S_decontam)
ps_16S_decontam <- prune_taxa(taxa_sums(ps_16S_decontam) > 0, ps_16S_decontam)

### Subset data ----
ps_16S_seed <- subset_samples(ps_16S_decontam,Category == "Seed")
ps_16S_seed <- prune_taxa(taxa_sums(ps_16S_seed) > 0, ps_16S_seed)

ps_16S_germ <- subset_samples(ps_16S_decontam,Category == "Germ")
ps_16S_germ <- prune_taxa(taxa_sums(ps_16S_germ) > 0, ps_16S_germ)

ps_16S_soil <- subset_samples(ps_16S_decontam,Category == "Soil")
ps_16S_soil <- prune_taxa(taxa_sums(ps_16S_soil) > 0, ps_16S_soil)

ps_16S_elev <- subset_samples(ps_16S_decontam,
                              plot_sample %in% c("P2S2","P3S2",
                                                 "P4S1","P6S1"))
ps_16S_elev <- prune_taxa(taxa_sums(ps_16S_elev) > 0, ps_16S_elev)                            
