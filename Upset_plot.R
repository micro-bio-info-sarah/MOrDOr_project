# Create upset plot

#https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
#https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html

ps = ps_16S_seed

###  Upset dataframe ----

#merged by plot_sample / plot / host / stimulant (induce NA in sample data)
ps_merged = merge_samples(ps,group = ps@sam_data$plot_sample)

# extraction of otu table
otudf = as.data.frame(as.matrix(otu_table(ps_merged)))
#create vector with samples name
sampnam <- sample_names(ps_merged)
#i need taxa in rows and sample names in columns
otudf <- transpose(otudf)
colnames(otudf) = sampnam
rownames(otudf) = taxa_names(ps_merged)
#create a presence/absence matrix
upsetdf <- as.data.frame(matrix(data = 0,
                                nrow = NROW(otudf),
                                ncol = NCOL(otudf),
                                byrow = FALSE,dimnames = NULL))
rownames(upsetdf) <- rownames(otudf)
colnames(upsetdf) <- colnames(otudf)
upsetdf[otudf > 0] <- 1 
#add taxa sums column
upsetdf$'ASV total abundance' <- taxa_sums(ps_merged)
upsetdf$'ASV total abundance (log10)' <- log10(taxa_sums(ps_merged))
#add prevalence number
upsetdf$prev = rowSums(upsetdf[,1:nsamples(ps_merged)])/nsamples(ps_merged)*100

### Upset plot ----

upsetplot <- upset(upsetdf,
                   sets = colnames(upsetdf)[1:(ncol(upsetdf)-3)],
                   order.by = "degree", 
                   decreasing ="FALSE",
                   nintersects =NA,
                   keep.order = T,
                   queries = list(list(query = intersects, 
                                       params = colnames(upsetdf)[20:26], 
                                       color = "red", active = F)),
                   mainbar.y.label = "ASV count", 
                   boxplot.summary = "ASV total abundance (log10)",
                   sets.x.label = "Total count of ASV")

