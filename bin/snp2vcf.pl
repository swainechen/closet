#!/usr/bin/perl -w
#
# convert mummer show-snps -T to vcf file
#
use slchen;
use File::Basename;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $dp = 1000;
my $qual = 1000;
my $all = 0;	# if all is 0, then only print out SNPs
my $original_command = "$0 " . join (" ", @ARGV);

GetOptions (
  "depth=i" => \$dp,
  "quality=i" => \$qual,
  "all!" => \$all
);
if (defined $ARGV[0] && -f $ARGV[0]) {
  $file = $ARGV[0];
} else {
  $file = "STDIN";
}
$sample = basename($file);
$sample =~ s/\.snps$//;
$sample =~ s/\.snp$//;

my $ref_fasta = "";

while (<>) {
  if (!length $ref_fasta || !-f $ref_fasta) {
    @f = split /\s+/, $_;
    $ref_fasta = $f[0] if -f $f[0];
  }
  last if /^NUCMER/;
}

my $sequence;
open F, $ref_fasta;
my @a = <F>;
$sequence = slchen::fasta2hash(@a);
close F;

while (<>) {
  last if /^\[P1\]/;
}
print "##fileformat=VCFv4.0\n";
print "##fileDate=", `date "+%Y%m%d"`;
print "##source=$original_command\n";
print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n";
print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
print "##FORMAT=<ID=FT,Number=2,Type=String,Description=\"Filter\">\n";
print "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$sample\n";

my $lastP1 = -1;
my $lastP2 = -1;
my $pos = -1;
my $ref = "";
my $alt = "";
my $chrom = "";
while (<>) {
  chomp;
  @f = split /\t/, $_;
  # if we're past an indel, print out the data
  if ($f[0] != $lastP1 && $f[3] != $lastP2 && $ref ne "" && $alt ne "" && $pos != -1 && $chrom ne "") {
    print join ("\t", $chrom, $pos, ".", $ref, $alt, $qual, "PASS", "DP=$dp", "GT:DP:FT", "1/1:$dp:PASS"), "\n";
    $chrom = "";
    $pos = -1;
    $ref = "";
    $alt = "";
  }
  if (!$all) {
    next if $f[1] eq ".";
    next if $f[2] eq ".";
    next if length $f[1] != 1;
    next if length $f[2] != 1;
  }
  $lastP1 = $f[0];
  $lastP2 = $f[3];
  $chrom = $f[8];
  if ($f[1] ne "." && $f[2] ne ".") {
    # snp, should be pretty safe
    $pos = $f[0];
    $ref = $f[1];
    $alt = $f[2];
  } elsif ($f[1] eq ".") {
    # insertion relative to reference - file gives 1-based position
    # vcf rules say to put in the base just before
    # position isn't going to change so we're safe here
    if (!defined $sequence->{$chrom}) {
      die "Can't find sequence in $ref_fasta for $chrom, position $f[0]\n";
    }
    if ($ref eq "" && $alt eq "") {
      $ref = substr($sequence->{$chrom}, $f[0]-1, 1);
      $alt = $ref;
      $pos = $f[0];
    }
    $alt .= $f[2];
  } else {
    # deletion relative to reference
    if (!defined $sequence->{$chrom}) {
      die "Can't find sequence in $ref_fasta for $chrom, position $f[0]\n";
    }
    if ($ref eq "" && $alt eq "") {
      # for the first one have to add the position before
      # and we need to catch the position here before it increments
      # take this chance to do a sanity check
      if ($f[1] ne substr($sequence->{$chrom}, $f[0]-1, 1)) {
        die "Discrepancy in $chrom ($ref_fasta), position $f[0]; I see " . substr($sequence->{$chrom}, $f[0]-1, 1) . " but file reports $f[1]\n";
      }
      $ref = substr($sequence->{$chrom}, $f[0]-2, 1);
      $alt = $ref;
      $pos = $f[0] - 1;
    }
    $ref .= $f[1];
  }
}
if ($lastP1 != -1 && $lastP2 && $ref ne "" && $alt ne "" && $pos != -1 && $chrom ne "") {
  print join ("\t", $chrom, $pos, ".", $ref, $alt, $qual, "PASS", "DP=$dp", "GT:DP:FT", "1/1:$dp:PASS"), "\n";
}
