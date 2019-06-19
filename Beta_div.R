# Estimate beta-diversity
setwd("~/MOrDOr_project")
require(phyloseq)
require(ggplot2)
require(vegan)
require(data.table)

#set the wanted phyloseq object: ps_16S/ITS_seed/soil/germ/elev
ps=ps_16S_seed

#first we'll normalize the abundance
ps_norm <- transform_sample_counts(ps,function(x){1E6*x/sum(x)})

### Ordination ----

ordin = ordinate(ps_norm, 
                 method = "PCoA",
                 distance = "bray")

plot_ordination(ps_norm, ordin, type = "samples", axes = 1:2,
                color = "host", shape = "plot",
                label = NULL, title = NULL,
                justDF = FALSE) + 
  scale_shape_manual(values = (1:nlevels(as.factor(ps_norm@sam_data[["plot"]])))) +
  geom_point(size = 5)+ theme_bw()+
  labs(shape="Plot",col="Host",title = "Fungi (ITS)")+
  scale_color_discrete(labels=c("Hemp","Oilseed rape","Tobacco"))+
  stat_ellipse(type = "t", linetype = 2) +
  theme(legend.title = element_text(size=16, face="bold"),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(face="bold", size=20), 
        axis.text.x  = element_text(vjust=0.5, size=20, angle = 90),
        axis.title.y = element_text(face="bold", size=20), 
        axis.text.y  = element_text(vjust=0.5, size=20),
        strip.text = element_text(angle=90,face="bold", size=20),
        strip.background = element_rect(fill="white", 
                                        colour="black",size=1.5))
                                        
### Statistics ----

#transpose otu table
metadata_otu <- as.data.frame(ps_norm@otu_table)
metadata_otu <- transpose(metadata_otu)
rownames(metadata_otu) <- sample_names(ps_norm)
colnames(metadata_otu) <- taxa_names(ps_norm)
#convert sample_data to data.frame
metadata <- as(sample_data(ps_norm), "data.frame")
metadata[,c(1,2,3,4,5,6,8,9)] <- lapply(metadata[,c(1,2,3,4,5,6,8,9)],factor)

dist_16S <- vegan::vegdist(metadata_otu, method = "bray")

# loop to print the r2 from adonis and
#the relative r2 (r2 / nlevels)
for (i in 1:ncol(metadata)){
  if (nlevels(as.factor(metadata[,i])) >= 2) {
    print(" ")
    lv_name = colnames(metadata)[i]
    print(lv_name)
    cap <- capscale(formula = dist_16S  ~ metadata[,lv_name],data = metadata)
    aov <- anova.cca(cap, stop=999)
    print(aov)
    ad <- adonis(formula = dist_16S ~ metadata[,lv_name],data = metadata, perm = 9999)
    r2 = ad$aov.tab[1,"R2"]
    nb_lv =  nrow(ad$coef.sites) #nlevels(as.factor(ps_norm@sam_data[,i]))
    nb_sp = nlevels(as.factor(ps_norm@sam_data$X.SampleID))
    r2_rel = (r2 / nb_lv)*100
    print(r2)
    print(paste0(round(r2,3),"/",nb_lv,"=",round(r2_rel,2)))
  } else { 
    print(paste0(colnames(metadata)[i],"Only one factor level for ",colnames(metadata)[i]))}
  
}
