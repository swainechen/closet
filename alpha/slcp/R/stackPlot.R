#' Stacked bar graph of contig classifications
#'
#'
#' @export
#'
stackPlot <- function(x, border=NA, col=c("gray", "red", "blue"), ...) {
  # this expects a data frame with output from ddply / contig_summary
  # for ex. all_contigs <- ddply(blastdata, .(Query, Contig), contig_summary)
  # columns expected are Chromosome, Plasmid, None, Length
  #
  tobar <- t(as.matrix(x[order(x$Chromosome, x$Plasmid),c("Chromosome", "Plasmid", "None")]))
  width <- x[order(x$Chromosome, x$Plasmid),"Length"]
  barplot(tobar, width=width, border=border, col=col)
  legend("bottomright", fill=col, legend=rownames(tobar))
}
