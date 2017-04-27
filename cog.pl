#!/usr/bin/perl -w
#
# take orgcode
# take list of genes on STDIN
# print out genome-wide cog makeup, then makeup of list of genes
# do simple binomial statistics for significant differences
#
$verbose = 0;	# if verbose, print out significantly enriched/depleted cogs
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
use Statistics::Distributions qw(uprob);
use Orgmap qw(:DEFAULT $pttfile $genefield);
&read_orgmap;

$exclude = "";
GetOptions (
  'v' => \$verbose,
  'exclude=s' => \$exclude
);

%exclude = ();
if (-f $exclude) {
  open E, $exclude;
  while (<E>) {
    chomp;
    next if /^$/ || /^#/;
    $exclude{$_} = 1;
  }
}

open PTT, $pttfile;
my @f;
my ($gid, $syst);
my %cog;
my @all = ();
my $cutoff = 0.05;
while (<PTT>) {
  if (/^\s*\d+\.\.\d+\s+/) {
    @f = split /\t/, $_;
    $gid = $f[3];
    $syst = $f[$genefield];
    next if $gid <= 0;
    next if $exclude{$gid};
    next if $exclude{$syst};
    if ($f[6]) {
      $cog{$gid} = $f[6];
      $cog{$syst} = $f[6];
    }
    push @all, $gid;
  }
}

@some = ();
while (<>) {
  chomp;
  next if /^#/;
  @f = split /\s+/;
  push @some, $f[0];
}

my %cogcats = ();
set_cog_names(\%cogcats);
my %cogall = ();
my %cogsome = ();
my %p = ();
&count_cog(\@all, \%cogall);
&count_cog(\@some, \%cogsome);
foreach my $code (keys %cogcats) {
  my $rate = $cogall{$code}/$cogall{'total'};
  my $trial = $cogsome{'total'};
  my $success = $cogsome{$code};
  if ($success > $rate * $trial) {	# want prob of this or greater
					# so get prob of success-1 or fewer
					# subtract from 1
    $p{$code} = 1 - binomial($rate, $trial, $success-1);
  } else {				# want prob of this or less
    $p{$code} = binomial($rate, $trial, $success);
  }
#  $p{$code} = binomial($cogall{$code}/$cogall{'total'}, $cogsome{'total'}, $cogsome{$code});
  $p{$code} = 1 - $p{$code} if $p{$code} > 0.5;

  printf "%s\t%d\t%d\t%.5f\t%d\t%d\t%.5f\t%.5f", $code, $cogall{$code}, $cogall{'total'}, $cogall{$code}/$cogall{'total'}, $cogsome{$code}, $cogsome{'total'}, $cogsome{$code}/$cogsome{'total'}, $p{$code};
  if ($verbose) {
    if ($cogall{$code}/$cogall{'total'} < $cogsome{$code}/$cogsome{'total'}) {
      if ($p{$code} < $cutoff) {
        print "\tEnriched";
      }
    } else {
      if ($p{$code} < $cutoff) {
        print "\tDepleted";
      }
    }
  }
  print "\n";
  if ($verbose && $p{$code} < $cutoff) {
    foreach my $item (@some) {
      if ($cog{$item} eq $code) {
        print "  $item\t$cogcats{$code}\n";
      }
    }
  }
}

sub count_cog {
  my ($a, $h) = @_;
  foreach my $code (keys %cogcats) {
    $h->{$code} = 0;
    $h->{'total'} = 0;
  }
  foreach my $item (@$a) {
    if (defined $cog{$item}) {
      $h->{$cog{$item}}++;
      $h->{'total'}++;
    }
  }
}

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

sub set_cog_names {
  # COG Categories at http://www.ncbi.nlm.nih.gov/sutils/coxik.cgi?gi=115
  # JAKLBDYVTMNZWUOCGEFHIPQRS-
  my ($ref) = @_;
  %$ref = (
    "J" => "Translation",
    "A" => "RNA processing and modification",
    "K" => "Transcription",
    "L" => "Replication, recombination and repair",
    "B" => "Chromatin structure and dynamics",
    "D" => "Cell cycle control, mitosis and meiosis",
    "Y" => "Nuclear structure",
    "V" => "Defense mechanisms",
    "T" => "Signal transduction mechanisms",
    "M" => "Cell wall/membrane biogenesis",
    "N" => "Cell motility",
    "Z" => "Cytoskeleton",
    "W" => "Extracellular structures",
    "U" => "Intracellular trafficking and secretion",
    "O" => "Posttranslational modification, protein turnover, chaperones",
    "C" => "Energy production and conversion",
    "G" => "Carbohydrate transport and metabolism",
    "E" => "Amino acid transport and metabolism",
    "F" => "Nucleotide transport and metabolism",
    "H" => "Coenzyme transport and metabolism",
    "I" => "Lipid transport and metabolism",
    "P" => "Inorganic ion transport and metabolism",
    "Q" => "Secondary metabolites biosynthesis, transport and catabolism",
    "R" => "General function prediction only",
    "S" => "Function unknown",
    "-" => "not in COGs",
  );
}
