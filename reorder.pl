#!/usr/bin/perl -w
#
# take a file and a column number
# reorder that file based on another file
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $order_file = "";
my $col = 0;
my %in;
my @f;
my @input;
my $delimiter = "\t";
my @comments = ();

GetOptions (
  'col=i' => \$col,
  'order=s' => \$order_file,
  'd=s' => \$delimiter
);

if (!-f $order_file) {
  print "Usage: $0 <file to sort> -order <file to get order from> [ -col <column number to sort on> ] [ -d <delimiter> ]\n";
  exit;
}

while (<>) {
  chomp;
  next if /^$/;
  if (/^#/) {
    push @comments, $_;
  } else {
    @f = split /$delimiter/, $_;
    $in{$f[$col]} = $_;
  }
}
close F;

open F, $order_file;
while (<F>) {
  chomp;
  next if /^$/;
  next if /^#/;
  if (defined $in{$_}) {
    print $in{$_}, "\n";
  } else {
    print "$_\t0\n";
  }
}
close F;
