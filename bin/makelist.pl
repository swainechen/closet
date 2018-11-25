#!/usr/bin/perl -w
#
# take list of numbers one per line
# take it to comma separated range we can use on command line
#
use slchen;
@a = ();
while (<>) {
  chomp;
  next if /^#/ || /^$/;
  push @a, $_;
}
print make_list(@a);
