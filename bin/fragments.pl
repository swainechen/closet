#!/usr/bin/perl -w
#
# take list of positions
# give fragment sizes
#
use slchen;

@pos = ();
while (<>) {
  chomp;
  next if /^$/;
  next if /^#/;
  next if /^>/;
  push @pos, $_;
}

@pos = sortu @pos;
foreach $i (0..$#pos-1) {
  print $pos[$i+1] - $pos[$i], "\n";
}
