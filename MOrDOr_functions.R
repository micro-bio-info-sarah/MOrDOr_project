# Functions implemented in MOrDOr project

require(gridBase)

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
