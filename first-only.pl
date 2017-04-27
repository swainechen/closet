#!/usr/bin/perl -w
#
# print only first occurrence of a given field
#
if ((defined $ARGV[0]) && (($ARGV[0] eq '-h') || ($ARGV[0] eq '--help'))) {
  print "Usage: $0 [ -d DELIMITER ] -c # INPUT\n";
  print "Prints only lines matching the first occurrence of values in a given column\n";
  print "  -d specifies delimiter, default is tab.\n";
  print "  -c specifies column to look at, 0-based.  Default is column 0 (the first).\n";
  exit (-1);
}

use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$delimiter = "\t";
$col = 0;
GetOptions ('d=s' => \$delimiter,
            'c=s' => \$col);

my $seen;
my @f;
while (<>) {
  if (/^#/ || /^$/) {
    print; next;
  }
  chomp;
  @f = split /$delimiter/, $_;
  if (!defined $seen->{$f[$col]}) {
    print "$_\n";
    $seen->{$f[$col]} = 1;
  }
}
