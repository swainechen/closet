#!/usr/bin/perl -w
#
# take a column of numbers, get a histogram
# input n = number of bins, col = column to read (default n = 100, col = 0)
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$n = 50;
$col = 0;
$percent = 0;
$help = 0;
GetOptions ('n=i' => \$n, 'col=i' => \$col, 'percent' => \$percent, 'help' => \$help);
if ($help) {
  print "Usage: hist-graph.pl [ -n bins ] [ -col column ] [ -percent ] <infile>\n";
  exit;
}
use PDL;
$PDL::IO::Misc::colsep = "\t";
if (defined $ARGV[0]) { $fh = $ARGV[0]; }
else { $fh = *STDIN; }
$x = rcols ($fh, $col, { EXCLUDE => '/^#/' });
($min, $max) = (stats($x))[3,4];
($bins, $hist) = hist($x, $min, $max, ($max-$min)/$n);
if ($percent) { $hist /= sum($hist); }
wcols $bins, $hist;
