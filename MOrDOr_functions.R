# Functions implemented in MOrDOr project

require(gridBase)

### Rarefaction curve ----
#From Martial Briand, IRHS Angers
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

### Taxonomic composition plot ----                   
#From Martial Briand, IRHS Angers
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

### Extract bootstraps for dendrograms ----
#https://stackoverflow.com/questions/22749634/how-to-append-bootstrapped-values-of-clusters-tree-nodes-in-newick-format-in
as.phylo.hclust.with.nodenames <- function (x, nodenames, ...) #We add a nodenames argument
{
  N <- dim(x$merge)[1]
  edge <- matrix(0L, 2 * N, 2)
  edge.length <- numeric(2 * N)
  node <- integer(N)
  node[N] <- N + 2L
  cur.nod <- N + 3L
  j <- 1L
  for (i in N:1) {
    edge[j:(j + 1), 1] <- node[i]
    for (l in 1:2) {
      k <- j + l - 1L
      y <- x$merge[i, l]
      if (y > 0) {
        edge[k, 2] <- node[y] <- cur.nod
        cur.nod <- cur.nod + 1L
        edge.length[k] <- x$height[i] - x$height[y]
      }
      else {
        edge[k, 2] <- -y
        edge.length[k] <- x$height[i]
      }
    }
    j <- j + 2L
  }
  if (is.null(x$labels)) 
    x$labels <- as.character(1:(N + 1))
  node.lab <- nodenames[order(node)] #Here we define our node labels
  obj <- list(edge = edge, edge.length = edge.length/2, tip.label = x$labels, 
              Nnode = N, node.label = node.lab) #And you put them in the final object
  class(obj) <- "phylo"
  reorder(obj)
}
