#!/usr/bin/perl
#
# take tab delimited, spit out fasta
#
use warnings;
use strict;

my @f;
while (<>) {
  if (/^#/) {
    print;
    next;
  }
  chomp;
  next if /^$/;
  @f = split /\t/, $_;
  if (scalar @f == 2) {
    print ">$f[0]\n$f[1]\n";
  } else {
    print "$_\n";
  }
}
