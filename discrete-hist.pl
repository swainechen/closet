#!/usr/bin/perl -w
#
# make discrete histogram - i.e. count number of times each label appears
#
use slchen;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$col = 0;
$percent = 0;
$cum = 0;
$delim = "\t";
$help = 0;
$blank = 0;	# whether to filter blanks
GetOptions (
  'col=s' => \$col,
  'percent' => \$percent,
  'cumulative' => \$cum,
  'blank!' => \$blank,
  'delimiter=s' => \$delim,
  'help' => \$help
);

if ($help) {
  print "Usage: $0 [ -col column(s) ] [ -percent ] [ -cumulative ] [ -d delimiter ] <infile>\n";
  exit;
}

$multifile = 0;
$stdin = 0;
if (@ARGV) {
  if ($#ARGV > 0) { $multifile = 1; }
  else { $fh = $ARGV[0]; }
} else {
  $fh = *STDIN;
  $stdin = 1;
}

$count = {};

if ($multifile) {
  foreach $i (0..$#ARGV) {
    $fh = $ARGV[$i];
    if (-f $fh) {
      open I, $fh;
      while (<I>) {
        chomp;
        next if /^#/ || /^$/;
        @f = split /\t/, $_;
        next if !defined $f[$col];
        next if $blank && $f[$col] eq "";
        if (defined $count->{$f[$col]}->[$i]) {
          $count->{$f[$col]}->[$i]++;
        } else {
          $count->{$f[$col]}->[$i] = 1;
        }
      }
      close I;
      $label[$i] = $fh;
    }
  }
} else {
  @cols = parse_list($col);
  if (!$stdin) {
    open I, $fh;
    $input = *I;
  } else {
    $input = *STDIN;
  }
  while (<$input>) {
    chomp;
    @f = split /\t/, $_;
    if ($_ =~ s/^#\s+//) {
      @label = @f[@cols];
      next;
    }
    foreach $j (0..$#cols) {
      next if !defined $f[$cols[$j]];
      next if $blank && $f[$cols[$j]] eq "";
      if (defined $count->{$f[$cols[$j]]}->[$j]) {
        $count->{$f[$cols[$j]]}->[$j]++;
      } else {
        $count->{$f[$cols[$j]]}->[$j] = 1;
      }
    }
  }
}

# data cleaning, put in 0's where there are undefineds
my $maxi = 0;
foreach $k (keys %$count) {
  $maxi = $#{$count->{$k}} if $#{$count->{$k}} > $maxi;
}
foreach $k (keys %$count) {
  foreach $i (0..$maxi) {
    if (!defined $count->{$k}->[$i]) {
      $count->{$k}->[$i] = 0;
    }
  }
}

if ($percent) {
  $total = [];
  foreach $k (keys %$count) {
    foreach $i (0..$#{$count->{$k}}) {
      if (!defined $total->[$i]) {
        $total->[$i] = $count->{$k}->[$i];
      } else {
        $total->[$i] += $count->{$k}->[$i];
      }
    }
  }
  foreach $i (0..$#$total) {
    if ($total->[$i]) {
      foreach $k (keys %$count) {
        $count->{$k}->[$i] /= $total->[$i] if $count->{$k}->[$i];
      }
    }
  }
}

if ($cum) {
  $running = [];
  foreach $k (sort { $a <=> $b } keys %$count) {
    foreach $i (0..$#{$count->{$k}}) {
      if (defined $count->{$k}->[$i]) {
        $running->[$i] += $count->{$k}->[$i];
        $count->{$k}->[$i] = $running->[$i];
      }
    }
  }
}

$numeric = 1;
foreach $k (keys %$count) {
  $numeric = 0 if !isfloat($k);
#  if ($k !~ /^\s*([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?\s*$/) {
#    $numeric = 0;
#  }
}

if (@label) {
  foreach $k (0..$#label) {
    if (defined $total->[$k]) {
      $label[$k] .= "(" . $total->[$k] . ")";
    }
  }
  print join ("\t", "# value", @label), "\n";
}
if ($numeric) {
  foreach $k (sort { $a <=> $b } keys %$count) {
    print join ("\t", $k, @{$count->{$k}}), "\n";
  }
} else {
  foreach $k (sort { $a cmp $b } keys %$count) {
    print join ("\t", $k, @{$count->{$k}}), "\n";
  }
}
