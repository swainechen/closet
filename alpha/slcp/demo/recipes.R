library(IRanges)
library(plyr)
library(ggplot2)
library(ape)
# need to figure out which actually are still needed
library(Biostrings)
library(biomartr)
library(slcp)

ref_fna_file = "/home/slchen/Documents/Projects/Plasmids/RC1/full/Gamma_plasmids.fna"
ref_meta_file = "/home/slchen/Documents/Projects/Plasmids/RC1/full/Gamma_plasmids_meta.txt"
#ref_fna_file = "/home/slchen/Documents/Projects/Plasmids/RC1/Gamma_plasmids.fna"
#ref_meta_file = "/home/slchen/Documents/Projects/Plasmids/RC1/Gamma_plasmids_meta.txt"
ref_meta <- read.table(ref_meta_file, sep="\t", quote="", header=T, comment.char="", stringsAsFactors=F)
rownames(ref_meta) <- ref_meta$Accession

# input assembly file, should have contigs
name <- "WBB857"
assembly <- "WBB857.assembly"
known_plasmid <- "/mnt/projects/slchen/recalibration/pacbio/GridION-181112/TTSH16-unicycler.fasta"
download_tempdir <- tempdir()

# predict plasmid contigs and their best hit
query_fasta <- readDNAStringSet(assembly)
blastn_data <- do_blast(assembly, name=name, database=ref_fna_file)
# TODO: problem here because need to capture sequence names during db generation
summary_class <- ddply(blastn_data, .(Query, Contig), contig_summary, .progress="text")
#summary_classification(blastn_data, ref_meta)
stackplot(summary_class, col=c("gray", "red", "blue"))
p_contigs <- subset(summary_class[order(summary_class$Length, decreasing=T),], Plasmid > 0.5)$Contig
plotDetail(contig_detail(subset(blastn_data, Contig==p_contigs[1])))

# analyze against a single plasmid or assembly
blastn <- do_blast(assembly, name=assembly, database=known_plasmid)
subject_meta <- ddply(blastn, .(Subject), summarize, Length=max(Slen), Type="None")
rownames(subject_meta) <- subject_meta$Subject
# go fix up the real types
subject_meta$Type <- c("Chromosome", "Plasmid", "None", "None", "None")
ref_meta <- subject_meta
summary_class <- ddply(blastn, .(Query, Contig), contig_summary, .progress="text")
plotDetail(subject_detail(subset(blastn, Subject=="2_plasmid")), classification=summary_class)

# analyze against a plasmid in Genbank
# for ex CP024856.1
accession <- "CP024856.1"
plasmid_blastn <- do_blast(assembly, name=name, database=accession)
plotDetail(subject_detail(plasmid_blastn))

# need to add gene annotation, from gff file, ? prokka to automate?
# need to add mouseover information
# dynamic zoom? shiny, D3
