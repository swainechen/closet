library(IRanges)
library(plyr)
library(ggplot2)
library(ape)
# need to figure out which actually are still needed
library(Biostrings)
library(biomartr)
library(slcp)

ref_fna_file <- "/tmp/plasmids/Gamma_plasmids.fna"
ref_meta_file <- "/tmp/plasmids/Gamma_plasmids_meta.txt"
ref_meta <- read.table(ref_meta_file, sep="\t", quote="", header=T, comment.char="", stringsAsFactors=F)
rownames(ref_meta) <- ref_meta$Accession

# input assembly file, should have contigs
name <- "WBB857"
assembly <- "WBB857.assembly"
known_plasmid <- "/mnt/projects/slchen/recalibration/pacbio/GridION-181112/TTSH16-unicycler.fasta"
download_tempdir <- tempdir()

# predict plasmid contigs and their best hit
query_fasta <- readDNAStringSet(assembly)
blastn_data <- doBlast(assembly, name=name, database=ref_fna_file)
# TODO: problem here because need to capture sequence names during db generation
summary_class <- ddply(blastn_data, .(Query, Contig), contigSummary, .progress="text")
#summary_classification(blastn_data, ref_meta)
stackPlot(summary_class, col=c("gray", "red", "blue"))
p_contigs <- subset(summary_class[order(summary_class$Length, decreasing=T),], Plasmid > 0.5)$Contig
plotDetail(contigDetail(subset(blastn_data, Contig==p_contigs[1])))
# logit_width stretches out the extremes more the closer it is to 1
# logit_width must be less than 1
plotDetail(contigDetail(subset(blastn_data, Contig==p_contigs[1])), logit_width=0.99)

# analyze against a single plasmid or assembly
blastn <- doBlast(assembly, name=assembly, database=known_plasmid)
subject_meta <- ddply(blastn, .(Subject), summarize, Length=max(Slen), Type="None")
rownames(subject_meta) <- subject_meta$Subject
# go fix up the real types
subject_meta$Type <- c("Chromosome", "Plasmid", "None", "None", "None")
ref_meta <- subject_meta
summary_class <- ddply(blastn, .(Query, Contig), contigSummary, .progress="text")
plotDetail(subjectDetail(subset(blastn, Subject=="2_plasmid")), classification=summary_class)

# analyze against a plasmid in Genbank
# for ex CP024856.1
accession <- "CP024856.1"
plasmid_blastn <- doBlast(assembly, name=name, database=accession)
plotDetail(subjectDetail(plasmid_blastn))

# need to add gene annotation, from gff file, ? prokka to automate?
# need to add mouseover information
# dynamic zoom? shiny, D3
