#!/usr/bin/perl -w
#
# Take vcf file
# fix remapping if it's there
# convert genotype calls to the actual SNP
#
use warnings;
use strict;
use Getopt::Long;
Getopt::Long::Configure("pass_through");
my $mincov = 10;
GetOptions (
  'mincov=i' => \$mincov
);


my $remap = ();
my $line = 0;
my $i;
my @f;
my @g;
my $ref;
my @alt;
my $allele;

while (<>) {
  $line++;
  if (/^## Name change/) {
    /^## Name change: Field \d+ has name longer than \d+; (\S+) changed to (Remap_\d+)$/;
    $remap->{$2} = $1;
    next;
  }
  if (/^##/) {
    print;
    next;
  }
  if (/^#CHROM/) {
    chomp;
    @f = split /\t/, $_;
    foreach $i (9..$#f) {
      if ($f[$i] =~ /^Remap_\d+$/ && defined $remap->{$f[$i]}) {
        $f[$i] = $remap->{$f[$i]};
      }
    }
    print join ("\t", @f), "\n";
    next;
  }
  chomp;
  @f = split /\t/, $_;
  $ref = $f[3];
  @alt = split /,/, $f[4];
  foreach $i (9..$#f) {
    @g = split /:/, $f[$i];
    $allele = substr($g[0], 0, 1);
    if ($g[$#g] ne "PASS" || $g[1] < $mincov || $allele eq "." || $allele !~ /\d+/) {
      $f[$i] = "-";
    } elsif ($allele == 0) {
      $f[$i] = $ref;
    } elsif ($allele-1 <= $#alt) {
      $f[$i] = $alt[$allele-1];
    } else {
      print STDERR "Had problem with assigning allele on line $line:\n$_\n";
    }
  }
  print join ("\t", @f), "\n";
}
