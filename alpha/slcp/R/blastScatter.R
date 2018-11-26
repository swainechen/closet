#' Scatter plot of slcp-formatted blastn results
#'
#' Can work on multiple query contigs and subject hits.
#'
#'
#'
#' @export
#'
blastScatter <- function (x, showbreaks=F, strandcolor=F, border=NA, ...) {
  # scatter plot for the blasts of the data frame
  # concatenate the subject refs
  # this currently also relies on the db data frame
  #
  subjects <- unique(x$Subject)
  refs <- db[subjects,]
  refs <- refs[order(refs$Type, refs$Length),]
  refs$Start <- 0
  refs$End <- cumsum(refs$Length)
  refs$Start <- c(1, refs$End[1:nrow(refs)-1] + 1)
  total <- max(refs$End)
  qs <- queries[queries$Query == unique(x$Query),]
  qs$Start <- 0
  qs$End <- cumsum(qs$Length)
  qs$Start <- c(1, qs$End[1:nrow(qs)-1] + 1)
  totalq <- max(qs$End)
  x$Start <- refs[x$Subject, "Start"] + pmin(x$Sstart, x$Send) - 1
  x$End <- refs[x$Subject, "Start"] + pmax(x$Sstart, x$Send) - 1
  x$Length <- x$End - x$Start + 1
  x$ID <- x$ID/100
  x$QstartA <- qs[paste(x$Query, x$Contig, sep="__"), "Start"] + x$Qstart - 1
  x$QendA <- qs[paste(x$Query, x$Contig, sep="__"), "Start"] + x$Qend - 1
  if (strandcolor) {
    x$StrandColor <- "green"
    x$StrandColor[which(x$Sstart > x$Send)] <- "red"
    x$StrandColor[which(x$Qstart > x$Qend)] <- "red"
  } else {
    x$StrandColor = "black"
  }
  v <- Rle(rep(0, total))
  for(i in 1:nrow(x)) {
    v <- pmax(v, Rle(c(0, x$ID[i], 0),
                     c(x$Start[i]-1, x$Length[i], total-x$End[i])))
  }
  plot.new()
  xlim <- c(0, total)
  ylim <- c(0, totalq)
  plot.window(xlim, ylim, main="Blast Scatter plot")
  rect(x$Start, x$QstartA, x$End, x$QendA, col=x$StrandColor, border=border, ...)
  axis(1, labels=F, at=c(0, refs$End), tick=T)
  axis(1, labels=refs$Subject, at=(refs$Start + refs$End)/2, tick=F)
  axis(2)
  title(main="Blast Scatter Plot", xlab = "Reference sequences", ylab = x[1,"Query"])
  abline(v=refs$End, lty=3)
  if (showbreaks) {
    abline(h=qs$End, lty=3)
  }
  if (strandcolor) {
    legend("topright", title="Alignment strand", legend=c("Plus/Plus", "Plus/Minus"), fill=c("green", "red"))
  }
}
