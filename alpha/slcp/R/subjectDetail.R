#' Return details of blast results for a single subject (database/reference)
#'
#' Requires a data frame of slcp-formatted blastn results for only one subject
#' that was hit. This might be useful for example for checking coverage of a
#' known reference plasmid from Genbank
#'
#' @export
#'

# TODO: if no info it assumes all hits are Plasmids
# We don't always have any info, should fall back to hit or not
# This doesn't deal well with multiple contigs on the subject side
# can pull in coordinate plots as in how SVRE does it

subject_detail <- function(x) {
  # switch the place of subject and query - i.e. now subject based
  # assume we only have one subject in the input data frame
  # expecting to get blast data, so columns should be:
  # c("Query", "Contig", "Subject", "ID", "AlignLen", "Mismatch", "Gaps", "Qstart", "Qend", "Sstart", "Send", "Evalue", "Bitscore", "Qlen", "Slen")
  # this also needs a queries variable to get info on the query sequences
  #
  queries <- ddply(x, .(Query, Contig), summarize, Length=head(Qlen, 1))
  rownames(queries) <- paste(queries$Query, queries$Contig, sep="__")
  queries$Type <- ""
  total <- x$Slen[1]
  start <- pmin(x$Sstart, x$Send)
  end <- pmax(x$Sstart, x$Send)
  length <- end - start + 1
  # these store the highest ID for a chrom/plasmid hit for each query coordinate
  v_Chrom <- Rle(rep(0, total))
  v_Plasm <- Rle(rep(0, total))
  v_Chrom_Query <- Rle(rep("", total))
  v_Plasm_Query <- Rle(rep("", total))
  for(i in 1:nrow(x)) {
    if (queries[paste(x$Query[i], x$Contig[i], sep="__"), "Type"] == "Chromosome") {
      this_hit <- Rle(c(0, x$ID[i], 0), c(start[i]-1, length[i], total-end[i]))
      v_Chrom_Query[this_hit > v_Chrom] <- paste(x$Query[i], x$Contig[i], sep="__")
      v_Chrom <- pmax(v_Chrom, this_hit)
    } else {
      this_hit <- Rle(c(0, x$ID[i], 0), c(start[i]-1, length[i], total-end[i]))
      v_Plasm_Query[this_hit > v_Plasm] <- paste(x$Query[i], x$Contig[i], sep="__")
      v_Plasm <- pmax(v_Plasm, this_hit)
    }
  }
  # convert these to IRanges to do figure out where there's hit to both
  # chromosome and plasmid for the same coordinates
  ir_Chrom <- IRanges(v_Chrom > 0)
  ir_Plasm <- IRanges(v_Plasm > 0)
  # final stores the final classification of whether that coordinate best hit
  # to a chromosome or plasmid sequence in the database
  final <- Rle(rep("", total))
  final_id <- Rle(rep(0, total))
  final_query <- Rle(rep("", total))
  perc_chrom <- length(which(final == "Chromosome"))/total
  # if only a chromosomal hit, those are by definition chromosomal
  if (length(setdiff(ir_Chrom, ir_Plasm)) > 0) {
    for(i in 1:length(setdiff(ir_Chrom, ir_Plasm))) {
      s <- start(setdiff(ir_Chrom, ir_Plasm))
      e <- end(setdiff(ir_Chrom, ir_Plasm))
      final[s[i]:e[i]] <- "Chromosome"
      final_id[s[i]:e[i]] <- v_Chrom[s[i]:e[i]]
      final_query[s[i]:e[i]] <- v_Chrom_Query[s[i]:e[i]]
    }
  }
  # if only a plasmid hit, those are by definition plasmid
  if (length(setdiff(ir_Plasm, ir_Chrom)) > 0) {
    for(i in 1:length(setdiff(ir_Plasm, ir_Chrom))) {
      s <- start(setdiff(ir_Plasm, ir_Chrom))
      e <- end(setdiff(ir_Plasm, ir_Chrom))
      final[s[i]:e[i]] <- "Plasmid"
      final_id[s[i]:e[i]] <- v_Plasm[s[i]:e[i]]
      final_query[s[i]:e[i]] <- v_Plasm_Query[s[i]:e[i]]
    }
  }
  # check the identity of the hit to chromosome and plasmid when there are both
  if (length(intersect(ir_Chrom, ir_Plasm)) > 0) {
    for(i in 1:length(intersect(ir_Chrom, ir_Plasm))) {
      # set these to chromosome first
      # where plasmid is higher then switch to plasmid
      s <- start(intersect(ir_Chrom, ir_Plasm))
      e <- end(intersect(ir_Chrom, ir_Plasm))
      final[s[i]:e[i]] <- "Chromosome"
      final_id[s[i]:e[i]] <- v_Chrom[s[i]:e[i]]
      final_query[s[i]:e[i]] <- v_Chrom_Query[s[i]:e[i]]
      final[s[i]:e[i]][v_Plasm[s[i]:e[i]] > v_Chrom[s[i]:e[i]]] <- "Plasmid"
      final_id[s[i]:e[i]][v_Plasm[s[i]:e[i]] > v_Chrom[s[i]:e[i]]] <- v_Plasm[s[i]:e[i]][v_Plasm[s[i]:e[i]] > v_Chrom[s[i]:e[i]]]
      final_query[s[i]:e[i]][v_Plasm[s[i]:e[i]] > v_Chrom[s[i]:e[i]]] <- v_Plasm_Query[s[i]:e[i]][v_Plasm[s[i]:e[i]] > v_Chrom[s[i]:e[i]]]
    }
  }
  return(list(final, final_id, final_query, x$Subject[1]))
}
