#' Blastn scatter plot with windowing
#'
#' Plot slcp-formatted blastn data frame, can do windowing.
#'
#'
#'
#'
#'
#' @export
#'
blastScatterWindow <- function (x, s = 1, e = 10000, border=NA, min.id = 0.9, ...) {
  # scatter plot for the blasts of the data frame
  # concatenate the subject refs
  # this currently also relies on the db data frame
  #
  subjects <- unique(x$Subject)
  refs <- db[subjects,]
  refs <- refs[order(refs$Type, -refs$Length),]
  refs$Y <- 1:nrow(refs)
  x <- subset(x, Qend >= s)
  x <- subset(x, Qstart <= e)
  x$SstartAdj <- x$Sstart
  x$SendAdj <- x$Send
  par(mar = c(5,8,4,3) + 0.1)
  plot.new()
  xlim <- c(0, 100)
  ylim <- c(0, max(refs$Y))
  height <- 0.8
  plot.window(xlim, ylim, main="Blast Scatter plot")
  axis(1, labels=T)
  axis(2, labels=F, at=0:max(refs$Y), tick=T)
  axis(2, labels=refs$Subject, at=refs$Y - 0.6, tick=F, las=1, cex.axis = 0.8) 
  title(main=paste(sep="\n", "Blast Scatter Plot", floor((s+e)/2)), xlab = "% Position on Reference sequence")
  title(ylab="Subject", line=6)
  if (nrow(x) > 1) {
    for(i in 1:nrow(x)) {
      if (x$Qstart[i] < s) {
        x$SstartAdj[i] <- x$Sstart[i] + round((s-x$Qstart[i]+1)/(x$Qend[i]-x$Qstart[i]+1) * (x$Send[i]-x$Sstart[i]+1))
      }
      if (x$Qend[i] > e) {
        x$SendAdj[i] <- x$Send[i] - round((x$Qend[i]-e+1)/(x$Qend[i]-x$Qstart[i]+1) * (x$Send[i]-x$Sstart[i]+1))
      }
    }
    x$X1 <- pmin(x$SstartAdj, x$SendAdj)
    x$X2 <- pmax(x$SstartAdj, x$SendAdj)
    x$X1 <- x$X1/x$Slen * 100
    x$X2 <- x$X2/x$Slen * 100
    x$Y <- refs[x$Subject, "Y"]
    x$ID <- x$ID/100
    x$Sat <- x$ID
    x$Sat[which(x$Sat < min.id)] <- min.id
    x$Sat <- (x$Sat - min.id) / (1 - min.id)
    rect(x$X1, x$Y-height, x$X2, x$Y, col=hsv(h=x$Y/max(x$Y), s=x$Sat, v=1), border=border, ...)
    points((x$X1+x$X2)/2, x$Y-height/2, col=hsv(h=x$Y/max(x$Y), s=x$Sat, v=1), cex=abs(x$X1-x$X2)/max(abs(x$X1-x$X2)), pch=16, ...)
  }
}
