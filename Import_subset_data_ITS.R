# Import and subset data for microbiota analysis

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
abund_ITS <- read.csv2("abundance_table_ITS.csv", sep = ',',row.names = 1)
taxo_ITS <- read.csv2("taxo_ITS.csv", sep = ',',row.names = 1, colClasses = "character")
taxo_ITS <- as.matrix(taxo_ITS)
design <- read_delim("design_ITS.csv",",", escape_double = F, trim_ws = T)
rownames(design) <- design$`#SampleID`
tree = read.newick("Qiime2_output/ITS/phylogenetic_tree/tree.nwk")

OTU = otu_table(abund_ITS, taxa_are_rows = TRUE)
TAX = tax_table(taxo_ITS)
DES = sample_data(design)
TREE = phy_tree(tree)

ps_ITS = phyloseq(OTU,TAX,DES,TREE)

### Extract ASV id ----
# extract seq ID
ASVid1 = taxa_names(ps_ITS)
#set new names
taxa_names(ps_ITS) <- paste0("ASV",1:ntaxa(ps_ITS),
                             "_",ps_ITS@tax_table[,"Genus"],
                             "_",ps_ITS@tax_table[,"Species"])
# extract new seq id
ASVid2 = taxa_names(ps_ITS)
ASVid = cbind(ASVid1,ASVid2)
write_csv(as.data.frame(ASVid),"ASVid_ITS.csv")

### Clean dataset ----
#remove mock and desinfected Pram 120
ps_ITS <- prune_samples(!(sample_data(ps_ITS)$plot_sample %in% c("Mock",
                                                                 "Pram_120",
                                                                 "Pram_120_d")),
                        ps_ITS)

#remove contaminants
#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html#identifying-contaminants-in-marker-gene-and-metagenomics-data
sample_data(ps_ITS)$is.neg <- sample_data(ps_ITS)$plot_sample %in% c("Control","H2O","PBS")
contamdf.eith <- isContaminant(ps_ITS, method="either", 
                               conc="quantif_DNA",
                               neg="is.neg", threshold=0.1)
table(contamdf.eith$contaminant)
head(which(contamdf.eith$contaminant))

ps.pa <- transform_sample_counts(ps_ITS, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$plot_sample %in% c("Control","H2O","PBS"), ps.pa)
ps.pa.pos <- prune_samples(!(sample_data(ps.pa)$plot_sample %in% c("Control","H2O","PBS")), ps.pa)

df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.eith$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

plot_frequency(ps_ITS,
               taxa_names(ps_ITS)[sample(which(contamdf.eith$contaminant),16)],
               conc="quantif_DNA") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

ps_ITS_decontam <- prune_taxa(!contamdf.eith$contaminant, ps_ITS)
ps_ITS_decontam <- prune_samples(!(sample_data(ps_ITS_decontam)$plot_sample %in% c("Control","H2O","PBS")),ps_ITS_decontam)
ps_ITS_decontam <- prune_taxa(taxa_sums(ps_ITS_decontam) > 0, ps_ITS_decontam)

### Subset data ----
ps_ITS_seed <- subset_samples(ps_ITS_decontam,Category == "Seed")
ps_ITS_seed <- prune_taxa(taxa_sums(ps_ITS_seed) > 0, ps_ITS_seed)

ps_ITS_germ <- subset_samples(ps_ITS_decontam,Category == "Germ")
ps_ITS_germ <- prune_taxa(taxa_sums(ps_ITS_germ) > 0, ps_ITS_germ)

ps_ITS_soil <- subset_samples(ps_ITS_decontam,Category == "Soil")
ps_ITS_soil <- prune_taxa(taxa_sums(ps_ITS_soil) > 0, ps_ITS_soil)

ps_ITS_elev <- subset_samples(ps_ITS_decontam,
                              plot_sample %in% c("P2S2","P3S2",
                                                 "P4S1","P6S1"))
ps_ITS_elev <- prune_taxa(taxa_sums(ps_ITS_elev) > 0, ps_ITS_elev)                            
