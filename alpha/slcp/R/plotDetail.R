#' Plot to visualize `contigDetail` or `subjectDetail`
#'
#'
#'
#'
#' @export
#'
plotDetail <- function(x, classification=NULL, border=NULL, minid = 0.7, ...) {
  # classification is like the output of ddply using contig_summary
  # it really needs Query, Contig that match Subject of x (which should have Query__Contig)
  # also classification needs Chromosome, Plasmid, None classifications (0-1)
  # and a total length for each Contig
  #
  # expect a list of 3 Rle elements
  # [[1]] is Plasmid/Chromosome, [[2]] is Identity, [[3]] is the Subject that had the best hit
  # [[4]] is name
  # first get ranges
  full_IR <- IRanges(start = c(start(x[[1]]), start(x[[2]]), start(x[[3]])), end = c(end(x[[1]]), end(x[[2]]), end(x[[3]])))
  segments <- disjoin(full_IR)
  ycoord <- data.frame(Subject = unique(x[[3]]), Y = 1:length(unique(x[[3]])), stringsAsFactors=F)
  ycoord <- ycoord[ycoord$Subject != "",]
  ycoord$maxLengthLeft <- 0
  for(i in 1:nrow(ycoord)) {
    tempintervals <- which(x[[3]]@values == ycoord$Subject[i])
    tempbest <- tempintervals[which.max(x[[3]]@lengths[tempintervals])]
    ycoord$maxLengthLeft[i] <- start(x[[3]])[tempbest]
    ycoord$Total[i] <- sum(x[[3]]@lengths[tempintervals])
  }
  rownames(ycoord) <- ycoord$Subject
# sort by length of subject
#  ycoord <- ycoord[order(-db[ycoord$Subject,]$Length),]
# sort by left coordinate of longest hit
  ycoord <- ycoord[order(ycoord$maxLengthLeft),]
  ycoord$Y <- 1:nrow(ycoord)
  plot.new()
  par(mar = c(5,10,4,3) + 0.1)
  height <- 1
  id_height <- 2
  sep <- 0.2
  xmax <- max(end(x[[1]]))
  xlim <- c(-0.07*xmax, 1.07*xmax)
  ylim <- c(-height/2, length(unique(x[[3]])) * (height + sep) + id_height)
  x[[1]][x[[1]]==""] <- "None"
  totalnt = max(end(x[[1]]))
  plasmidnt <- length(IRanges::which(x[[1]]=="Plasmid"))
  chromosoment <- length(IRanges::which(x[[1]]=="Chromosome"))
  nonent <- length(IRanges::which(x[[2]]==0))
  hue <- c("Plasmid" = 0, "Chromosome" = 0.5)
  plot.window(xlim, ylim)
  for(i in 1:length(segments)) {
    type <- unique(x[[1]][start(segments)[i]:end(segments)[i]])
    id <- unique(x[[2]][start(segments)[i]:end(segments)[i]])
    subject <- unique(x[[3]][start(segments)[i]:end(segments)[i]])
    if (id > 0 & subject != "None") {
      if (id > 1) { id <- id/100 }
      if (id < minid) { id <- minid }
      ybottom <- ycoord[subject,"Y"] * (height + sep) - height + id_height
      rect(start(segments)[i]-0.5, ybottom, end(segments)[i]+0.5, ybottom + height, col=hsv(h=hue[type], s=(id-minid)/(1-minid), v=1), border=border, ...)
    }
  }
  barcolor <- c("Plasmid" = hsv(h=hue["Plasmid"], s=1, v=1),
                "Chromosome" = hsv(h=hue["Chromosome"], s=1, v=1),
                "None" = hsv(h=0, s=0, v=0))
  for(i in 1:length(start(x[[1]]))) {
    rect(start(x[[1]])[i], -height/2, end(x[[1]])[i], -sep, col=barcolor[x[[1]]@values[i]], border=NA)
  }
  lines(x[[2]]/100 * id_height)
  xaxis_at <- axTicks(1)
  xaxis_at <- xaxis_at[which(xaxis_at < max(end(x[[1]])))]
  xaxis_at <- c(xaxis_at, max(end(x[[1]])))
  axis(1, at=xaxis_at)
  axis(2, labels=F, at=c(id_height, (ycoord$Y) * (height + sep) + id_height), tick=T)
  axis(2, labels=c("Source call", "Best Hit ID", ycoord$Subject), at=c(-(height/2+sep)/2, id_height/2, (ycoord$Y - 0.5) * (height + sep) + id_height), tick=F, las=1, cex.axis = 0.8) 
  if(!is.null(classification)) {
    rownames(classification) <- paste(classification$Query, classification$Contig, sep="__")
    ycoord$Total <- ycoord$Total / classification[ycoord$Subject,"Length"]
    for(i in 1:nrow(ycoord)) {
      if(ycoord[i,"Subject"] %in% rownames(classification)) {
        ybottom <- ycoord[i,"Y"] * (height + sep) - height + id_height
        rect(-0.07*xmax, ybottom, -0.052*xmax, ybottom + height, col=hsv(h=hue["Plasmid"], s=classification[ycoord$Subject[i],"Plasmid"], v=1), border=border, ...)
        rect(-0.05*xmax, ybottom, -0.032*xmax, ybottom + height, col=hsv(h=hue["Chromosome"], s=classification[ycoord$Subject[i],"Chromosome"], v=1), border=border, ...)
        rect(-0.03*xmax, ybottom, -0.012*xmax, ybottom + height, col=hsv(h=0, s=0, v=1-(classification[ycoord$Subject[i],"None"])), border=border, ...)
        rect(1.02*xmax, ybottom, (1.02 + 0.05*ycoord$Total[i])*xmax, ybottom + height, col="black", border=border, ...)
      }
    }
    segments(1.07*xmax, id_height, 1.07*xmax, max(ycoord$Y) * (height + sep) + id_height, lty="dotted", lwd=0.5)
    axis(3, labels=c("Plasmid", "Chromosome", "None"), at=c(-0.061*xmax, -0.041*xmax, -0.021*xmax), las=2, cex.axis=0.8)
    axis(1, labels="% rep", at=1.045*xmax, las=2, cex.axis=0.8)
  }
  axis(4, at=c(0, id_height), labels=range(x[[2]]/100), tick=T, cex.axis=0.6, las=1)
  title(main=paste(sep="\n", "Contig Detail plot", x[[4]]), xlab="Contig coordinate")
  legendvector <- c(sprintf("Plasmid (%.2f%%)", 100*plasmidnt/totalnt),
                    sprintf("Chromosome (%.2f%%)", 100*chromosoment/totalnt),
                    sprintf("None (%.2f%%)", 100*nonent/totalnt))
  title(ylab="Subject", line=6)
  legend("bottomright", legend=legendvector, fill=barcolor, cex=0.6, inset=c(0.1, 0.1))
  return(c("Total" = totalnt, "Plasmid" = plasmidnt, "Chromosome" = chromosoment, "None" = nonent))
}
