#!/usr/bin/perl -w
#
# delete fasta sequences less than a certain length
#

use Getopt::Long;

my $length = 0;
my $name = "";
my $seq = "";

GetOptions (
  "length=i" => \$length
);

while (<>) {
  next if /^#/;
  chomp;
  if (/^>/) {
    if ($name ne "" || $seq ne "") {
      print "$name\n$seq\n" if length($seq) >= $length;
    }
    $name = $_;
    $seq = "";
  } else {
    $seq .= $_;
  }
}
if ($name ne "" || $seq ne "") {
  print "$name\n$seq\n" if length($seq) >= $length;
}
