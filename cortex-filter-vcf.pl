#!/usr/bin/perl -w
#
# filter FINAL vcf from cortex
# i.e. take only PASS lines
# and eliminate anything that doesn't have variable SNP calls across samples
# VCF format should be:
# 0 - CHROM
# 1 - POS
# 2 - ID
# 3 - REF
# 4 - ALT
# 5 - QUAL
# 6 - FILTER
# 7 - INFO
# 8 - FORMAT
# 9 and on - samples
#
my @f;
my @g;
my $i;
my $different;
my $compare;

while (<>) {
  if (/^#/) { print; next; }
  chomp;
  @f = split /\t/, $_;
  next if $f[6] ne 'PASS';
  $different = 0;
  if ($#f <= 9) {
    $different = 1;
  } else {
    @g = split /:/, $f[9];
    $compare = $g[0];
    foreach $i (10..$#f) {
      @g = split /:/, $f[$i];
      if ($g[0] ne $compare) {
        $different = 1;
        last;
      }
    }
  }
  if ($different) {
    print $_, "\n";
  }
}
