Estimate alpha-diversity

#set the wanted phyloseq object: ps_16S/ITS_seed/soil/germ/elev
ps=ps_16S_seed

### Rarefaction ----

require(gridBase)

rarecurve2 <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", label = TRUE, col = "black", ...)
  ## See documentation for vegan rarecurve, col is now used to define
  ## custom colors for lines and panels
{
  tot <- rowSums(x)
  S <- specnumber(x)
  nr <- nrow(x)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_len(length(out))) {
    color <- col[((ln-1) %% length(col)) + 1]
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = color, ...)
  }
  if (label) {
    ordilabel(cbind(tot, S), labels = rownames(x), col = col, ...)
  }
  invisible(out)
}

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
