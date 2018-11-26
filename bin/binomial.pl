#!/usr/bin/perl -w
use Statistics::Distributions qw(uprob);

$leux_expect = 21493/(6180+16823+78357+17942+21351+21493);
$leuz_expect = 21351/(6180+16823+78357+17942+21351+21493);
$serx_expect = 24327/(24327+14312+12032+13580+13567+13449);
print join ("\t", "# systematic", "leuX cumu prob", "total leucine", "serX cumu prob", "total serine", "trivial", "annotation"), "\n";
while (<>) {
  chomp;
  next if /^#/;
  @f = split /\t/, $_;
  print $f[0], "\t";
  print binomial($leux_expect, $f[3], $f[1]*$f[3]), "\t$f[3]\t";
  print binomial($leuz_expect, $f[3], $f[2]*$f[3]), "\t$f[3]\t";
  print binomial($serx_expect, $f[5], $f[4]*$f[5]), "\t$f[5]\t";
  print join ("\t", $f[5], $f[6]), "\n";
}
#while (<>) {
#  chomp;
#  @f = split /\s+/, $_;
#  print binomial($f[0], $f[1], $f[2]), "\n";
#}

sub binomial {
  # take expected success rate, number of trials, number of successes, return
  # prob of getting that many successes or less
  my ($rate, $trial, $success) = @_;
  my $cum = 0;
  if ($trial > 50) {
    if ($trial * $rate >= 6 && $trial * (1-$rate) >= 10) {
      # use normal approxmation
      my $mean = $trial * $rate;
      my $stdev = sqrt ($trial * $rate * (1-$rate));
      if ($success < $trial) {
        # use continuity correction
        my $z = ($success + 0.5 - $mean) / $stdev;
        return (uprob($z));
      } else {
        return 1;
      }
    } elsif ($trial * $rate < 6) {
      # use poisson approximation
      return (poisson($trial*$rate, $success));
    } else {
      # get probability of number of failures, using poisson approximation
      return (1-poisson($trial*(1-$rate), $trial-$success));
    }
  }
  foreach my $i (0..$success) {
    $cum += choose($trial, $i) * $rate**$i * (1-$rate)**($trial-$i);
  }
  return ($cum);
}

sub poisson {
  # $l is lambda, i.e. expected occurences
  # $n is actual occurrences
  my ($p, $n) = @_;
  my $poisson;
  my $cum = 0;
  foreach my $i (0..$n) {
    $cum += exp(-$p) * ($p**$i) / fact($i);
  }
  if ($cum >= 1) { return 1; }
  else { return ($cum); }
}

sub choose {
  my ($n, $m) = @_;
  return (fact($n)/fact($m)/fact($n-$m));
}

sub fact {
  my ($n) = @_;
  if ($n <= 1) { return 1; }
  my $cum = 1;
  foreach my $i (2..$n) {
    $cum *= $i;
  }
  return $cum;
}
