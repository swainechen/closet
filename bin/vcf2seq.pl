#!/usr/bin/perl -w
#
# take vcf file
# give only reduced sequences indicated
# insert - for blanks
# options for no indels and biallelic only (default)
# useful for making SNP-only alignments
#
my $biallelic_only = 1;
my $indels = 0;
my @header = ();
my @alleles = ();
my $is_indel;
my $line = 0;
my $sequence;
while (<>) {
  $line++;
  chomp;
  if (/^#CHROM/) {
    @f = split /\t/, $_;
    @header = @f;
    foreach $i (9..$#f) {
      $sequence->{$i} = "";
    }
  }
  next if /^#/;
  # example line: 1 8577 . C T 49314 PASS AC=2;AN=1980 GT:DP:FT 0/0:50:PASS
  # ref is field (0-based) 3, alt is 4, genotypes start at 9
  @f = split /\t/, $_;
  @g = split /,/, $f[4];
  next if length($f[3]) > 1 && !$indels;
  next if $#g != 0 && $biallelic_only;
  $is_indel = 0;
  foreach $i (0..$#g) {
    $is_indel++ if length($g[$i]) > 1;
  }
  next if $is_indel && !$indels;
  @alleles = ($f[3], @g);
  foreach $i (9..$#f) {
    @h = split /[:\/]/, $f[$i];
    if ($#h < 1 || $h[0] ne $h[1]) {
      die "Error in genotype on line $line, field $i (0-based)\n";
    }
    if ($h[0] eq ".") {
      $sequence->{$i} .= "-" x length($f[3]);
    } else {
      $sequence->{$i} .= $alleles[$h[0]];
    }
  }
}
foreach $i (9..$#f) {
  print ">$header[$i]\n$sequence->{$i}\n";
}
