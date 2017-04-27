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
undef $max;
undef $min;
GetOptions (
  'max=f' => \$max,
  'min=f' => \$min,
  'n=i' => \$n,
  'col=s' => \$col,
  'percent' => \$percent,
  'cumulative' => \$cum,
  'help' => \$help
);
if ($help) {
  print "Usage: dual-hist.pl [ -n bins ] [ -col column ] [ -percent ] [ -cumulative ] [ -max <max value> ] [ -min <min value> ] <infile>\n";
  exit;
}

use PDL;
$PDL::IO::Misc::colsep = "\t";

$multifile = 0;
$stdin = 0;
if (@ARGV) {
  if ($#ARGV > 0) { $multifile = 1; }
  else { $fh = $ARGV[0]; }
} else {
  $fh = *STDIN;
  $stdin = 1;
}

if ($multifile) {
  foreach $i (0..$#ARGV) {
    $fh = $ARGV[$i];
    $label[$i] = $fh;
    $x[$i] = rcols ($fh, $col, { EXCLUDE => '/^#/' });
    if( defined $min ){
      $temp = $x[$i]->copy;
      $tmp = $temp->where($temp >= $min); # get values > min
      $x[$i] = $tmp->copy;
    }
    if( defined $max){
      $temp = $x[$i]->copy;
      $tmp = $temp->where($temp <= $max); # get values > max 
      $x[$i] = $tmp->copy;
    }
  }
} else {
  @cols = split /,/, $col;
  foreach $i (reverse 0..$#cols) {
    if ($cols[$i] =~ /(\d+)[-.]+(\d+)/) {
      splice @cols, $i, 1, ($1..$2);
    }
  }
#  @x = rcols ($fh, @cols, { EXCLUDE => '/^#/' });
  if (!$stdin) {
    open F, $fh;
    @data = <F>;
    close F;
  } else {
    @data = <>;
  }
  foreach $i (0..$#data) {
    chomp $data[$i];
    @f = split /\t/, $data[$i];
    if ($data[$i] =~ s/^#\s+//) {
      @label = @f[@cols];
      next;
    }
    foreach $j (0..$#cols) { 
      next if !defined $f[$cols[$j]] || $f[$cols[$j]] !~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/;	# we need numbers, throw away blanks also
      next if (defined $max && $f[$cols[$j]] >= $max);
      next if (defined $min && $f[$cols[$j]] <= $min);
      push @{"x$j"}, $f[$cols[$j]];
    }
  }
  foreach $j (0..$#cols) {
    $x[$j] = pdl @{"x$j"};
  }
}

$p = $x[0]->copy;
foreach $i (1..$#x) {
  $p = append ($p, $x[$i]);
}
($min, $max) = (stats($p))[3,4];

foreach $i (0 .. $#x) {
  ($bins, $y[$i]) = hist($x[$i], $min, $max, ($max-$min)/$n);
  if ($percent) {
    $y[$i] /= sum ($y[$i]);
    $label[$i] .= "(" . nelem($x[$i]) . ")";
  }
  if ($cum) { $y[$i] = cumusumover $y[$i]; }
}
if (scalar @label) {
  print join ("\t", "# bin", @label), "\n";
}
wcols $bins, @y;
