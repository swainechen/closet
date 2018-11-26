#' Summarize blast hits for a given contig (query)
#'
#' Needs slcp-formatted data frame with data from only one contig from one
#' query genome.
#'
#' Internal procedure, generally used in conjunction with ddply.
#'
#' @export
#'
contigSummary <- function(x, meta=ref_meta) {
  # assume we only have one query genome and contig in the input data frame
  # expecting to get blast data, so columns should be:
  # c("Query", "Contig", "Subject", "ID", "AlignLen", "Mismatch", "Gaps", "Qstart", "Qend", "Sstart", "Send", "Evalue", "Bitscore", "Qlen", "Slen")
  # i.e. there needs to be an extra "Query" column
  # for example Query might be the strain ID (WBB1160), then Contig is the
  # actual contig from the assembly, of which there may be many for that strain
  # this needs some metadata defined, which has info on the sequences in the
  # blast reference database, for example:
  # NC_003385.1     106516  Plasmid
  # these should be in the rowname, Length, and Type columns of a data frame
  #
  total <- x$Qlen[1]
  start <- pmin(x$Qstart, x$Qend)
  end <- pmax(x$Qstart, x$Qend)
  length <- end - start + 1
  v_Chrom <- Rle(rep(0, total))
  v_Plasm <- Rle(rep(0, total))
  for(i in 1:nrow(x)) {
    if (meta[x$Subject[i], "Type"] == "Chromosome") {
      v_Chrom <- pmax(v_Chrom, Rle(c(0, x$ID[i], 0),
                                   c(start[i]-1, length[i], total-end[i])))
    } else {
      v_Plasm <- pmax(v_Plasm, Rle(c(0, x$ID[i], 0),
                                   c(start[i]-1, length[i], total-end[i])))
    }
  }
  ir_Chrom <- IRanges(v_Chrom > 0)
  ir_Plasm <- IRanges(v_Plasm > 0)
  final <- Rle(rep("", total))
  if (length(IRanges::setdiff(ir_Chrom, ir_Plasm)) > 0) {
    for(i in 1:length(IRanges::setdiff(ir_Chrom, ir_Plasm))) {
      s <- start(IRanges::setdiff(ir_Chrom, ir_Plasm))
      e <- end(IRanges::setdiff(ir_Chrom, ir_Plasm))
      final[s[i]:e[i]] <- "Chromosome"
    }
  }
  if (length(IRanges::setdiff(ir_Plasm, ir_Chrom)) > 0) {
    for(i in 1:length(IRanges::setdiff(ir_Plasm, ir_Chrom))) {
      s <- start(IRanges::setdiff(ir_Plasm, ir_Chrom))
      e <- end(IRanges::setdiff(ir_Plasm, ir_Chrom))
      final[s[i]:e[i]] <- "Plasmid"
    }
  }
  if (length(IRanges::intersect(ir_Chrom, ir_Plasm)) > 0) {
    for(i in 1:length(IRanges::intersect(ir_Chrom, ir_Plasm))) {
      # set these to chromosome first
      # where plasmid is higher then switch to plasmid
      s <- start(IRanges::intersect(ir_Chrom, ir_Plasm))
      e <- end(IRanges::intersect(ir_Chrom, ir_Plasm))
      final[s[i]:e[i]] <- "Chromosome"
      final[s[i]:e[i]][v_Plasm[s[i]:e[i]] > v_Chrom[s[i]:e[i]]] <- "Plasmid"
    }
  }
  perc_chrom <- length(IRanges::which(final == "Chromosome"))/total
  perc_plasm <- length(IRanges::which(final == "Plasmid"))/total
  perc_none <- length(IRanges::which(final == ""))/total
  best_df <- ddply(x, .(Subject), summarize, TopBitscore = max(Bitscore), TotalBitscore=sum(Bitscore), Longest = max(AlignLen*ID/100), Total = sum(AlignLen*ID/100))
  best_df$Type <- meta[best_df$Subject, "Type"]
  if(perc_chrom >= perc_plasm) {
    if (nrow(subset(best_df, Type == "Chromosome")) > 0) {
      best_df <- subset(best_df, Type == "Chromosome")
    }
  } else {
    if (nrow(subset(best_df, Type == "Plasmid")) > 0) {
      best_df <- subset(best_df, Type == "Plasmid")
    }
  }
  return(data.frame(Chromosome=perc_chrom,
           Plasmid=perc_plasm,
           None=perc_none,
           Length=total,
           BestSingle=best_df[which.max(best_df$TopBitscore), "Subject"],
           BestTotal=best_df[which.max(best_df$TotalBitscore), "Subject"]))
}
