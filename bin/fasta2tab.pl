#!/usr/bin/perl -w
#
# take fasta
# spit out tab delimited
#
use warnings;
use strict;

my $name = "";
my @out = ();
while (<>) {
  if (/^#/) {
    print;
    next;
  }
  chomp;
  next if /^$/;
  if (s/^>//) {
    if (scalar @out) {
      print join ("\t", @out), "\n";
    }
    @out = ($_);
  } else {
    push @out, $_;
  }
}
if (scalar @out) {
  print join ("\t", @out), "\n";
}
