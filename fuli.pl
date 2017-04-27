#!/usr/bin/perl -w
#
# calculate Fu and Li's D* and F*
# take aligned fasta sequence as input
#
my $name;
my $seq;
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

my $n;
my $eta;
my $eta_s;
my $Dstar;
my $pi_n;
my $Fstar;

$n = scalar (keys %$seq);
($eta, $eta_s) = eta($seq);
$pi_n = average_difference($seq);

$Dstar = $n / ($n-1) * $eta - a($n) * $eta_s;
$Dstar /= sqrt(u_D($n) * $eta + v_D($n) * $eta * $eta);

$Fstar = $pi_n - ($n - 1) / $n * $eta_s;
$Fstar /= sqrt(u_F($n) * $eta + v_F($n) * $eta * $eta);

print "Number of sequences: $n\n";
print "Average differences: $pi_n\n";
print "eta: $eta\n";
print "eta_s: $eta_s\n";
print "D* = $Dstar\n";
print "F* = $Fstar\n";

sub a {
  my ($n) = @_;
  my $i;
  my $a = 0;
  foreach $i (1..$n-1) {
    $a += 1/$i;
  }
  return $a;
}

sub b {
  my ($n) = @_;
  my $i;
  my $b = 0;
  foreach $i (1..$n-1) {
    $b += 1/($i*$i);
  }
  return $b;
}

sub c {
  my ($n) = @_;
  my $c;
  if ($n == 2) {
    $c = 1;
  } elsif ($n > 2) {
    $c = 2 * ($n * a($n) - 2*($n-1)) / ($n-1) / ($n-2);
  }
  return $c;
}

sub d {
  my ($n) = @_;
  my $d = c($n);
  $d += ($n - 2) / ($n - 1) / ($n - 1);
  $d += 2 / ($n - 1) * (3/2 - (2*a($n+1) - 3)/($n-2) - 1/$n);
  return $d;
}

sub v_D {
  my ($n) = @_;
  my $v = $n*$n/($n-1)/($n-1) * b($n);
  $v += a($n) * a($n) * d($n);
  $v -= 2 * $n * a($n) * (a($n)+1) / ($n - 1) / ($n - 1);
  $v /= a($n) * a($n) + b($n);
  return $v;
}

sub u_D {
  my ($n) = @_;
  my $u = $n / ($n-1);
  $u *= a($n) - $n/($n-1);
  $u -= v_D($n);
  return $u;
}

sub v_F {
  my ($n) = @_;
  my $v = d($n);
  $v += 2 * ($n*$n + $n + 3) / 9 / $n / ($n-1);
  $v -= 2 / ($n-1) * (4 * b($n) - 6 + 8/$n);
  $v /= a($n) * a($n) + b($n);

  $v = 2*$n*$n*$n + 110*$n*$n - 255*$n + 153;
  $v /= 9 * $n * $n * ($n-1);
  $v += 2 * ($n-1) * a($n) / $n / $n;
  $v -= 8 * b($n) / $n;
  $v /= a($n) * a($n) + b($n);
  return $v;
}

sub u_F {
  my ($n) = @_;
  my $u = $n / ($n - 1);
  $u += ($n + 1) / 3 / ($n - 1);
  $u -= 2 * 2 / $n / ($n-1);
  $u += 2 * ($n + 1) / ($n - 1) / ($n - 1) * (a($n+1) - 2 * $n / ($n + 1));
  $u /= a($n);
  $u -= v_F($n);

  $u = 4*$n*$n + 19*$n + 3 - 12*($n+1)*a($n+1);
  $u /= 3 * $n * ($n - 1);
  $u /= a($n);
  $u -= v_F($n);
  return $u;
}

sub eta {
  my ($seq) = @_;
  my $name;
  my $i;
  my $data;
  my $eta = 0;
  my $eta_s = 0;
  foreach $name (keys %$seq) {
    foreach $i (0..length($seq->{$name})-1) {
      $data->{$i}->{uc substr($seq->{$name}, $i, 1)}++;
    }
  }
  foreach $i (keys %$data) {
    foreach $name (keys %{$data->{$i}}) {
      $eta_s++ if $data->{$i}->{$name} == 1;
      $eta++;
    }
    $eta-- if keys %{$data->{$i}};	# eta is number of nt at each position minus 1
  }
  return ($eta, $eta_s);
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
