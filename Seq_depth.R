# Estimate sequencing depth
setwd("~/MOrDOr_project")
require(phyloseq)
require(ggplot2)

#set the wanted phyloseq object: ps_16S/ITS_seed/soil/germ/elev
ps=ps_16S_seed

sdt <- data.table(as(sample_data(ps), "data.frame"),
                  TotalReads = sample_sums(ps), 
                  keep.rownames = TRUE)

pSeqDepth <- ggplot(sdt, aes(TotalReads)) + 
  geom_histogram() + 
  ggtitle("Sequencing Depth")

pSeqDepth + labs(x="Sample size (read total)",y= "Number of samples") +
  scale_x_log10()
