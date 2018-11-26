#' Return details of blast results for a single contig (query)
#'
#' Requires a data frame of slcp-formatted blastn results for only one contig
#' from one query genome. This might be useful for example for known complete
#' plasmid or chromosome on a single contig.
#'
#'
#'
#' @export
#'
contigDetail <- function(x, meta=ref_meta) {
  # this is more for whole analysis of a chromosome
  # assume we only have one query genome and contig in the input data frame
  # expecting to get blast data, so columns should be:
  # c("Query", "Contig", "Subject", "ID", "AlignLen", "Mismatch", "Gaps", "Qstart", "Qend", "Sstart", "Send", "Evalue", "Bitscore", "Qlen", "Slen")
  # this also needs a db variable defined to get info on sequences in the
  # reference blast database
  #
  total <- x$Qlen[1]
  start <- pmin(x$Qstart, x$Qend)
  end <- pmax(x$Qstart, x$Qend)
  length <- end - start + 1
  # these store the highest ID for a chrom/plasmid hit for each query coordinate
  v_Chrom <- Rle(rep(0, total))
  v_Plasm <- Rle(rep(0, total))
  v_Chrom_Subject <- Rle(rep("", total))
  v_Plasm_Subject <- Rle(rep("", total))
  for(i in 1:nrow(x)) {
    if (meta[x$Subject[i], "Type"] == "Chromosome") {
      this_hit <- Rle(c(0, x$ID[i], 0), c(start[i]-1, length[i], total-end[i]))
      v_Chrom_Subject[this_hit > v_Chrom] <- x$Subject[i]
      v_Chrom <- pmax(v_Chrom, this_hit)
    } else {
      this_hit <- Rle(c(0, x$ID[i], 0), c(start[i]-1, length[i], total-end[i]))
      v_Plasm_Subject[this_hit > v_Plasm] <- x$Subject[i]
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
  final_subject <- Rle(rep("", total))
  # if only a chromosomal hit, those are by definition chromosomal
  if (length(setdiff(ir_Chrom, ir_Plasm)) > 0) {
    for(i in 1:length(setdiff(ir_Chrom, ir_Plasm))) {
      s <- start(setdiff(ir_Chrom, ir_Plasm))
      e <- end(setdiff(ir_Chrom, ir_Plasm))
      final[s[i]:e[i]] <- "Chromosome"
      final_id[s[i]:e[i]] <- v_Chrom[s[i]:e[i]]
      final_subject[s[i]:e[i]] <- v_Chrom_Subject[s[i]:e[i]]
    }
  }
  # if only a plasmid hit, those are by definition plasmid
  if (length(setdiff(ir_Plasm, ir_Chrom)) > 0) {
    for(i in 1:length(setdiff(ir_Plasm, ir_Chrom))) {
      s <- start(setdiff(ir_Plasm, ir_Chrom))
      e <- end(setdiff(ir_Plasm, ir_Chrom))
      final[s[i]:e[i]] <- "Plasmid"
      final_id[s[i]:e[i]] <- v_Plasm[s[i]:e[i]]
      final_subject[s[i]:e[i]] <- v_Plasm_Subject[s[i]:e[i]]
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
      final_subject[s[i]:e[i]] <- v_Chrom_Subject[s[i]:e[i]]
      final[s[i]:e[i]][v_Plasm[s[i]:e[i]] > v_Chrom[s[i]:e[i]]] <- "Plasmid"
      final_id[s[i]:e[i]][v_Plasm[s[i]:e[i]] > v_Chrom[s[i]:e[i]]] <- v_Plasm[s[i]:e[i]][v_Plasm[s[i]:e[i]] > v_Chrom[s[i]:e[i]]]
      final_subject[s[i]:e[i]][v_Plasm[s[i]:e[i]] > v_Chrom[s[i]:e[i]]] <- v_Plasm_Subject[s[i]:e[i]][v_Plasm[s[i]:e[i]] > v_Chrom[s[i]:e[i]]]
    }
  }
  perc_chrom <- length(which(final == "Chromosome"))/total
  return(list(final, final_id, final_subject, paste(sep=" ", x$Query[1], x$Contig[1])))
}
