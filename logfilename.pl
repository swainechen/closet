#!/usr/bin/perl
#
# log command history
#
use warnings;
use strict;

my $date = `date +%F`;
chomp $date;

my $i = 1;
while (-f "$date"."_$i.log") {
  $i++;
}

print "$date", "_$i.log";
