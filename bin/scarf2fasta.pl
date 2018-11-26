#!/usr/bin/perl -w
#
# take solexa scarf format
# throw away quality and just get fasta
#
while (<>) {
  chomp;
  @f = split /:/, $_;
  print ">", join (":", @f[0..4]), "\n";
  print $f[5], "\n";
}
