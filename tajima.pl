#!/usr/bin/perl -w
#
# take a fasta alignment
# calculate some parameters
# give Tajima's D value
# ref: Tajima, F. (1989) Genetics 123:585-95
#
use Math::GSL::SF qw(:gamma);
my $name;
my $seq;
my $debug = 0;
while (<>) {
  next if /^#/;
  chomp;
  if (/^>/) {
    $name = $_;
  } elsif (defined $name) {
    $seq->{$name} = $_;
    undef $name;
  }
}

$numsamples = scalar(keys %$seq);
$segsites = segregating_sites($seq);
$avgdiff = average_difference($seq);
($tajimad, $p) = tajimas_d($numsamples, $segsites, $avgdiff);

print "num samples: $numsamples\nseg sites: $segsites\naverage diff: $avgdiff\ntajima's D: $tajimad\np-value: $p\n";

sub tajimas_d {
  my ($n, $S, $k) = @_;
  return 0 if !$S;
  my ($a1, $a2, $b1, $b2, $c1, $c2, $e1, $e2, $Dmin, $Dmax, $D);
  my $i;

  # a1 is sum of 1/i for i from 1 to n-1
  # a2 is sum of 1/i^2 for i from 1 to n-1
  $a1 = 0;
  $a2 = 0;
  foreach $i (1..$n-1) {
    $a1 += 1/$i;
    $a2 += 1/($i*$i);
  }
  $debug && print "a1 = $a1\na2 = $a2\n";

  # b1 is (n+1)/3(n-1)
  $b1 = ($n + 1) / (3 * ($n-1));
  # b2 is 2(n^2 + n + 3)/9n(n-1)
  $b2 = 2 * ($n*$n + $n + 3) / (9 * $n * ($n-1));
  $debug && print "b1 = $b1\nb2 = $b2\n";

  # c1 is b1 - 1/a1
  $c1 = $b1 - 1/$a1;
  # c2 is b2 - (n+2)/a1*n + a2/a1^2
  $c2 = $b2 - ($n+2)/($a1*$n) + $a2/($a1*$a1);
  $debug && print "c1 = $c1\nc2 = $c2\n";

  # e1 is c1/a1
  $e1 = $c1/$a1;
  # e2 is c2/(a1^2 + a2)
  $e2 = $c2/($a1*$a1 + $a2);
  $debug && print "e1 = $e1\ne2 = $e2\n";

  # min is (2/n - 1/a1) / sqrt(e2)
  $Dmin = (2/$n - 1/$a1) / sqrt($e2);
  # max is ( n+1/2n - 1/a1 ) / sqrt(e2)
  $Dmax = (($n+1)/(2*$n) - 1/$a1) / sqrt($e2);
  $debug && print "Dmin = $Dmin\nDmax = $Dmax\n";

  # now we can calculate D
  $D = ($k - $S/$a1) / sqrt($e1*$S + $e2*$S*($S-1));

  # get p value, use 1000 samples of distribution
  $p_less = tajima_beta($Dmin, $Dmax, 1000, $D);

  return ($D, $p_less);
}

sub segregating_sites {
  my ($seq) = @_;
  my $nt;
  my $length;
  my $i;
  my $sites = 0;
  foreach $name (keys %$seq) {
    @{$nt->{$name}} = split //, $seq->{$name};
    $length = length $seq->{$name} if !defined $length;
    die "sequence length ", length $seq->{name}, " but expected $length\n" if $length != length $seq->{$name};
  }
  foreach $i (0..$length-1) {
    $consensus = "";
    foreach $name (keys %$nt) {
      $consensus = $nt->{$name}->[$i] if $consensus eq "";
      if ($nt->{$name}->[$i] ne $consensus) {
        $sites++;
        last;
      }
    }
  }
  return $sites;
}

sub average_difference {
  my ($seq) = @_;
  my @name = keys %$seq;
  my $i;
  my $j;
  my $n = 0;
  my $total = 0;
  foreach $i (0..$#name-1) {
    foreach $j ($i+1..$#name) {
      $n++;
      $total += differences($seq->{$name[$i]}, $seq->{$name[$j]});
    }
  }
  return $total/$n if $n;
}

sub differences {
  my ($a, $b) = @_;
  # take 2 sequences
  # give # of differences
  # they need to be the same length and aligned already
  my $diff = 0;
  my $i;
  my @a = split //, $a;
  my @b = split //, $b;
  die "sub differences, lengths different:\n$a\n$b\n" if scalar @a != scalar @b;
  foreach $i (0..$#a) {
    $diff++ if $a[$i] ne $b[$i];
  }
  return $diff;
}

sub tajima_beta {
  # given Dmin, Dmax, and sampling resolution, approximate a p-value for
  # a given D value
  my ($min, $max, $sample, $D) = @_;
  my $alpha = - (1 + $min*$max)*$max / ($max - $min);
  my $beta = (1 + $min*$max)*$min / ($max - $min);
  my $constant = gsl_sf_gamma($alpha + $beta)/gsl_sf_gamma($alpha)/gsl_sf_gamma($beta)/($max-$min)**($alpha+$beta-1);
  my $prob_less = 0;
  my $prob_more = 0;
  my $probability;
  foreach $i (0..$sample) {
    $d = $min + $i*($max - $min)/$sample;
    $d = $max if $d > $max;
    $probability = $constant * ($max - $d)**($alpha-1) * ($d - $min)**($beta-1);
    if ($d <= $D) {
      $prob_less += $probability/$sample*($max-$min);
    } else {
      $prob_more += $probability/$sample*($max-$min);
    }
  }
  $norm = $prob_less + $prob_more;
  $debug && print "beta normalization: $norm\n";
  $prob_less /= $norm;
  $prob_more /= $norm;
  return ($prob_less);
}
