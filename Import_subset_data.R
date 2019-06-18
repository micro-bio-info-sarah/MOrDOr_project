# Import and subset data for microbiata analysis

setwd("~/MOrDOr_project")

abund_16S <- read.csv2("Sarah_results/abundance_table_16S.csv", sep = ',',row.names = 1)
taxo_16S <- read.csv2("Sarah_results/taxo_16S_silva.csv", sep = ',',row.names = 1, colClasses = "character")
taxo_16S <- as.matrix(taxo_16S)
design <- read_delim("Sarah_results/design_16S.csv",",", escape_double = F, trim_ws = T)
rownames(design) <- design$`#SampleID`
tree = read.newick("Sarah_results/16S/phylogenetic_tree/tree.nwk")

OTU = otu_table(abund_16S, taxa_are_rows = TRUE)
TAX = tax_table(taxo_16S)
DES = sample_data(design)
TREE = phy_tree(tree)

ps_16S = phyloseq(OTU,TAX,DES,TREE)
