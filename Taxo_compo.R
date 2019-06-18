# Taxonomic composition overview

#set the wanted phyloseq object: ps_16S/ITS_seed/soil/germ/elev
ps=ps_16S_seed

#set host as factor
ps@sam_data$host <- as.factor(ps@sam_data$host)
levels(ps@sam_data$host) <- c("Hemp","Oilseed rape","Tobacco")

### Graphical representation: barplot ----

## Plot composition at taxaRank2 level within taxa taxaSet1 at taxaRank1 level
## Restricts plot to numberOfTaxa taxa
plot_composition <- function(physeq, taxaRank1 = "Phylum", taxaSet1 = "Proteobacteria",
                             taxaRank2 = "Family", numberOfTaxa = 9, fill = NULL,
                             x = "Sample", y = "Abundance", facet_grid = NULL) {
  ## Args:
  ## - physeq: phyloseq class object
  ## - taxaRank1: taxonomic level in which to do the first subsetting
  ## - taxaSet1: subset of level taxaRank1 to use
  ## - taxaRank2: taxonomic level used to agglomerate
  ## - numberOfTaxa: number of (top) taxa to keep at level taxaRank2
  ##
  ## Returns:
  ## - ggplot2 graphics
  ggdata <- ggformat(physeq, taxaRank1, taxaSet1, taxaRank2, numberOfTaxa)
  p <- ggplot(ggdata, aes_string(x = x, y = y, fill = fill, color = fill, group = "Sample"))
  ## Manually change color scale to assign grey to "Unknown" (if any)
  if (!is.null(fill) && any(c("Unknown", "Other") %in% unique(ggdata[, fill]))) {
      ranks <- as.character(unique(ggdata[, fill]))
      ranks <- ranks[ ! ranks %in% c("Unknown", "Other")]
      colvals <- c(gg_color_hue(length(ranks)), "grey45", "black")
      names(colvals) <- c(ranks, "Unknown", "Other")
      ## Now add the manually re-scaled layer with Unassigned as grey
      p <- p + scale_fill_manual(values=colvals) + scale_color_manual(values = colvals)

  }
  p <- p + geom_bar(stat = "identity", position = "stack",colour="grey")
  if ( !is.null(facet_grid)) {
    p <- p + facet_grid(facets = facet_grid, scales = "free_x")
  }
  p <- p + theme(axis.text.x=element_text(angle=90), axis.title.x=element_blank())
  p <- p + ggtitle(paste("Composition within", taxaSet1, "(", numberOfTaxa, "top", taxaRank2, ")"))
  return(p)
}

#set the wanted population (Bacteria / Fungi) and taxa rank (Phylum / Genus)
plot_composition(ps, "Kingdom", "Fungi", "Genus", 
                 numberOfTaxa=7, fill="Genus") + 
  facet_wrap(~host, scales="free_x", nrow=1)+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())

### Graphical representation: boxplot ----

# agglomerate taxa and create dataframe from phyloseq object
#set the wanted taxa rank (Phylum / Genus)
taxo_plot <- psmelt(tax_glom(ps, taxrank = "Phylum"))
# convert "tax_rank" to a character vector from a factor
taxo_plot$Phylum <- as.character(taxo_plot$Phylum)

# group dataframe by "tax_rank", calculate abundance median 
medians <- ddply(taxo_plot, ~Phylum,
                 function(x) cbind(median(x$Abundance),
                                   mean(x$Abundance),
                                   sd(x$Abundance),
                                   sum(x$Abundance),
                                   (sum(x$Abundance)/sum(taxa_sums(ps))*100),
                                   sum(x$Abundance > 0),
                                   (sum(x$Abundance > 0)/nsamples(ps)*100)))
colnames(medians) <- c("Phylum","median","mean","sd",
                       "sum","sum_%","prev","prev_%")

# find "tax_rank" whose median is less than 100
remainder <- medians[medians$median <= 100,]$Phylum
# change their name to "Remainder"
taxo_plot[taxo_plot$Phylum %in% remainder,]$Phylum <- 'Remainder'

# boxplot
ggplot(taxo_plot,aes(x=Phylum,y=Abundance,fill= host)) + 
  geom_boxplot(data = taxo_plot, outlier.size=0.5, stat = "boxplot") +
  coord_flip() + scale_y_log10()+
  labs(y="Abundance",fill="Type",title = "Bacteria (16S)")+
  scale_fill_discrete(labels=c("Hemp","Oilseed rape","Tobacco"))+
  theme(legend.title = element_text(size=16, face="bold"),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(face="bold", size=20), 
        axis.text.x  = element_text(vjust=0.5, size=10),
        axis.title.y = element_text(face="bold", size=20), 
        axis.text.y  = element_text(vjust=0.5, size=10),
        strip.text = element_text(angle=45,face="bold", size=16),
        strip.background = element_rect(fill="white", colour="black",size=1.5)) 
