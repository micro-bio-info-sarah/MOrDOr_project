# Represent prevalence and abundance of ASV
setwd("~/MOrDOr_project")
require(phyloseq)
require(plyr)
require(ggplot2)

#set the wanted phyloseq object: ps_16S/ITS_seed/soil/germ/elev
ps=ps_16S_seed

# Create table, number of features for each phyla
table(tax_table(ps)[, "Kingdom"], exclude = NULL)
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prev_seed <- apply(X = otu_table(ps),
                   MARGIN = ifelse(taxa_are_rows(ps),
                                   yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})


# Add taxonomy and total read counts to this data.frame
prev_seed <- data.frame(Prevalence = prev_seed,
                        TotalAbundance = taxa_sums(ps),
                        tax_table(ps))
plyr::ddply(prev_seed,
            "Phylum",
            function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#Graphical representation of phylum abundance/prevalence
ggplot(prev_seed, aes(TotalAbundance,
                      Prevalence / nsamples(ps),color=Genus)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + 
  xlab("Total Abundance") + 
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + 
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = -0.2))
