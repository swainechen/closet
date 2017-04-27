#!/usr/bin/perl -w
#
# run srst2 for E. coli specific stuff
# for MLST just keep the MLST line
# others take the full gene outputs
#
use warnings;
use strict;
use File::Basename;
use File::Temp;
use File::Spec;
use Cwd;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $name = "";
my $q1 = "";
my $q2 = "";
my $combined_fasta = "/home/slchen/Documents/Projects/SRST2/Ecoli-srst2.fasta";
my $mlst_fasta = "/home/slchen/Documents/Projects/SRST2/MLST/Escherichia_coli#1.fasta";
my $mlst_def = "/home/slchen/Documents/Projects/SRST2/MLST/ecoli.txt";
my $mlst_delimiter = "-";
my $samtools = "/mnt/software/stow/samtools-0.1.18/bin/samtools";
my $bowtie2 = "/mnt/software/stow/bowtie2-2.2.4/bin/bowtie2";
my $bowtie2_build = "/mnt/software/stow/bowtie2-2.2.4/bin/bowtie2-build";

GetOptions (
  'name=s' => \$name,
  'q1=s' => \$q1,
  'q2=s' => \$q2,
  'mlst_db=s' => \$mlst_fasta,
  'mlst_def=s' => \$mlst_def,
  'mlst_delimiter=s' => \$mlst_delimiter
);

if (!length($name) || !length($q1) || !length($q2) || !-f $q1 || !-f $q2) {
  print "Usage: $0 -q1 <fastq1> -q2 <fastq2> -name <sample name>\n";
  print "Runs SRST2 with E. coli specific databases\n";
  print "fastq1 and fastq2 can be gzipped or not\n";
  exit;
}

if (!File::Spec->file_name_is_absolute($q1)) {
  $q1 = File::Spec->rel2abs($q1);
}
if (!File::Spec->file_name_is_absolute($q2)) {
  $q2 = File::Spec->rel2abs($q2);
}

my $currentdir = getcwd;
my $tempdir = File::Temp::tempdir( CLEANUP => 1 );
#my $oldpath = $ENV{'PATH'};
#$ENV{'PATH'} = "/usr/local/bin:/usr/bin:$oldpath";
# SRST2 0.1.8 enables environment variables to set programs
$ENV{'SRST2_SAMTOOLS'} = $samtools;
$ENV{'SRST2_BOWTIE2'} = $bowtie2;
$ENV{'SRST2_BOWTIE2_BUILD'} = $bowtie2_build;
chdir($tempdir);

# these filenames will follow the SRST2 naming system
my ($clean_q1, $clean_q2);
if ($q1 =~ /\.gz$/) {
  $clean_q1 = $name . "_1.fastq.gz";
} else {
  $clean_q1 = $name . "_1.fastq";
}
if ($q2 =~ /\.gz$/) {
  $clean_q2 = $name . "_2.fastq.gz";
} else {
  $clean_q2 = $name . "_2.fastq";
}

# we are using a local srst2 due to library differences - and we need to mangle the path for it so it gets the right bowtie2 also
#my $resistance = "/home/slchen/Documents/Projects/SRST2/resistance/ARGannot.fasta";
#my $VF = "/home/slchen/Documents/Projects/SRST2/VFDB/Escherichia_VF_clustered.fasta";
#my $serotype = "/home/slchen/Documents/Projects/SRST2/serotype/OH.fasta";

`ln -s $q1 $tempdir/$clean_q1`;
`ln -s $q2 $tempdir/$clean_q2`;
my $output;
#$output = `srst2 --input_pe $clean_q1 $clean_q2 --log --output $name-resistance --gene_db $resistance`;
#$output = `srst2 --input_pe $clean_q1 $clean_q2 --log --output $name-VFDB --gene_db $VF`;
#$output = `srst2 --input_pe $clean_q1 $clean_q2 --log --output $name-serotype --gene_db $serotype`;
$output = `srst2 --input_pe $clean_q1 $clean_q2 --log --output $name-combined --gene_db $combined_fasta 2>&1`;
$output = `srst2 --input_pe $clean_q1 $clean_q2 --log --output $name-mlst --mlst_db $mlst_fasta --mlst_definitions $mlst_def --mlst_delimiter "$mlst_delimiter" 2>&1`;

my $mlst_final = `tail -n +2 -q $name-mlst__*__results.txt | sed -e 's/^/# /'`;
my $fullgenes_final = `tail -n +2 -q $name-combined__fullgenes__*__results.txt`;

#$ENV{'PATH'} = $oldpath;
chdir($currentdir);
print $mlst_final, $fullgenes_final;
