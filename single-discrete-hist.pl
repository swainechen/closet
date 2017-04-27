#!/usr/bin/perl -w
#
# make discrete histogram - i.e. count number of times each label appears
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$col = 0;
$percent = 0;
$cum = 0;
$delim = "\t";
$help = 0;
GetOptions ('col=s' => \$col, 'percent' => \$percent, 'cumulative' => \$cum, 'delimiter=s' => \$delim, 'help' => \$help);

if ($help) {
  print "Usage: $0 [ -col column(s) ] [ -percent ] [ -cumulative ] [ -d delimiter ] <infile>\n";
  exit;
}

while (<>) {
  next if /^#/;
  chomp;
  @f = split /$delim/, $_;
  next if !defined $f[$col];
  if (defined $count{$f[$col]}) {
    $count{$f[$col]}++;
  } else {
    $count{$f[$col]} = 1;
  }
}

if ($percent) {
  $total = 0;
  foreach $k (keys %count) {
    $total += $count{$k};
  }
  if ($total) {
    foreach $k (keys %count) {
      $count{$k} /= $total;
    }
  }
}

if ($cum) {
  $running = 0;
  foreach $k (sort { $a <=> $b } keys %count) {
    $running += $count{$k};
    $count{$k} = $running;
  }
}

$numeric = 1;
foreach $k (keys %count) {
  if ($k !~ /^\s*([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?\s*$/) {
    $numeric = 0;
  }
}
if ($numeric) {
  foreach $k (sort { $a <=> $b } keys %count) {
    print "$k\t$count{$k}\n";
  }
} else {
  foreach $k (sort { $a cmp $b } keys %count) {
    print "$k\t$count{$k}\n";
  }
}
