#' Blast against an slcp-formatted fasta database
#'
#' Manages temp files and making a blast database
#'
#'
#' @export
#'
doBlast <- function(file, name="", database=ref_fna_file) {
  # figure out what we have
  # database could have been passed a filename or an accession number
  if(!file.exists(database)) {
    # the only other option right now is we assume it's a Genbank accession
    tempfilename <- paste(download_tempdir, "/", database, ".fasta", sep="")
    write.dna(read.GenBank(database), file=tempfilename, format="fasta", nbcol=1, colw=60)
    if(!file.exists(tempfilename)) {
      return(NA)
    } else {
      db_file <- tempfilename
    }
  } else {
    # we have a file that was passed in
    db_file <- database
  }
  # check for blast indices
  if(!file.exists(paste(db_file, ".nin", sep="")) &
     !file.exists(paste(db_file, ".00.nin", sep=""))) {
    system(paste("makeblastdb -dbtype nucl -in", db_file, sep=" "))
  }
  evalue_cutoff <- 1e-10
  command <- paste(sep=" ",
                   "blastn -db", db_file,
                   "-query", file,
                   "-out - -evalue", evalue_cutoff,
                   "-dust no",
                   '-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"',
                   "-max_target_seqs 20")
  blastcon <- pipe(command, open="r")
  blastdata <- read.table(blastcon, header=F, stringsAsFactors=F, sep="\t", quote="")
  close(blastcon)
  colnames(blastdata) <- c("Contig", "Subject", "ID", "AlignLen", "Mismatch", "Gaps", "Qstart", "Qend", "Sstart", "Send", "Evalue", "Bitscore", "Qlen", "Slen")
  if (name == "") {
    blastdata$Query <- basename(file)
  } else {
    blastdata$Query <- name
  }
  blastdata <- blastdata[,c(ncol(blastdata),1:(ncol(blastdata)-1))]
  return(blastdata)
}
