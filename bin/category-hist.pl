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
$cum = 0;
$help = 0;
$category = -1;
undef $max;
undef $min;
GetOptions (
  'max=f' => \$max,
  'min=f' => \$min,
  'n=i' => \$n,
  'col=i' => \$col,
  'percent' => \$percent,
  'cumulative' => \$cum,
  'category=i' => \$category,
  'help' => \$help
);
if ($help || $category < 0) {
  print "Usage: $0 [ -n bins ] [ -col column ] [ -percent ] [ -cumulative ] [ -max <max value> ] [ -min <min value> ] [ -category <category column> ] <infile>\n";
  print "Data will be partitioned according to category column.\n";
  print "Can only handle one column for -col option.\n";
  exit;
}

use PDL;
$PDL::IO::Misc::colsep = "\t";

my @data = <>;
my %label = ();
my $j = 0;

foreach $i (0..$#data) {
  chomp $data[$i];
  next if $data[$i] =~ /^#/;	# we ignore labels since we have categories
  @f = split /\t/, $data[$i];
  if (!defined $label{$f[$category]}) {
    $label{$f[$category]} = $j;
    $j++;
  }
  next if !defined $f[$col] || $f[$col] !~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/;
  next if (defined $max && $f[$col] > $max);
  next if (defined $min && $f[$col] < $min);
  push @{"x".$label{$f[$category]}}, $f[$col];
}
foreach $j (keys %label) {
  $x[$label{$j}] = pdl @{"x".$label{$j}};
}

$p = $x[0]->copy;
foreach $i (1..$#x) {
  $p = append ($p, $x[$i]);
}
($min, $max) = (stats($p))[3,4];

@label = sort { $label{$a} <=> $label{$b} } keys %label;
foreach $i (0 .. $#x) {
  ($bins, $y[$i]) = hist($x[$i], $min, $max, ($max-$min)/$n);
  if ($percent) {
    $y[$i] /= sum ($y[$i]);
    $label[$i] .= "(" . nelem($x[$i]) . ")";
  }
  if ($cum) { $y[$i] = cumusumover $y[$i]; }
}
print join ("\t", "# bin", @label), "\n";
wcols $bins, @y;
