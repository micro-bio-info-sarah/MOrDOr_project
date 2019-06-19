# Estimate alpha-diversity
setwd("~/MOrDOr_project")
require(phyloseq)
require(ggplot2)
source("MOrDOr_functions.R")

#set the wanted phyloseq object: ps_16S/ITS_seed/soil/germ/elev
ps=ps_16S_seed

### Rarefaction ----

#rarefaction curve
rarecurve2(t(otu_table(ps)), step=100,
          xlab="Sample size (read total)",ylab="Detected ASVs", label=FALSE,
          xlim=c(0,13000))
abline(v = 1500)
### 
#1500 pour 16S seed / 2000 pour ITS seed
#2500 pour 16S soil / 2000 pour ITS soil
#2000 pour 16S germ / 1500 pour ITS germ
#1500 pour 16S elev / 1500 pour ITS elev

# rarefy datasets
ps_rare <- rarefy_even_depth(ps, sample.size = 1500, rngseed = 712)

### Graphical viexw of alpha-diversity ----
p <- plot_richness(ps_rare, 
                   color="host", 
                   measures=c("Observed", "InvSimpson", "Shannon"), 
                   x ="plot_sample")

p + geom_boxplot(data=p$data, aes(x=plot_sample,y=value,color=host),
                 alpha=0.1)+
  labs(x="Sample",color="Host:",title = "Bacteria (16S)")+ theme_bw()+
  theme(legend.title = element_text(size=16, face="bold"),
        legend.text = element_text(size = 16),
        legend.position = "bottom",
        axis.title.x = element_text(face="bold", size=20), 
        axis.text.x  = element_text(angle=90,vjust=0.5, size=10),
        axis.title.y = element_text(face="bold", size=20), 
        axis.text.y  = element_text(vjust=0.5, size=10),
        strip.text = element_text(angle=0,face="bold", size=16),
        strip.background = element_rect(fill="white", colour="black",size=1.5))+
  scale_color_discrete(labels=c("Hemp","Oilseed rape","Tobacco"))
