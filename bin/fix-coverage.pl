#!/usr/bin/perl -w
#
# use gcovs.gz or gcov.gz files to get coverage
# take sample names from vcf file
# replace no-call positions with the right info
# default PASS
# change chromosome to numbers
#
use Getopt::Long;
Getopt::Long::Configure("pass_through");
my $verbose = 0;
my $do_remap = 1;	# whether to mangle long sample names
my $do_gcov = 1;	# whether to look up gcov files
			# if no, then just fix the header
my $data_dir = ".";	# where to look for gcov files
my $min_cov = 10;	# based on lofreq default
my $number_chromosomes = 1;	# change chromsome names to numbers for SNPhylo
my $target = "";
my $chrom = ();
my $chrom_counter = 1;
my $max_name_length = 10;
my @map = ();

GetOptions (
  "datadir=s" => \$data_dir,
  "verbose" => \$verbose,
  "remap!" => \$do_remap,
  "number_chromosomes!" => \$number_chromosomes,
  "gcov!" => \$do_gcov
);

while (<>) {
  if (/^##/) {
    print; next;
  }
  if (/^#CHROM/) {
    chomp;
    @map = split /\t/, $_;
    @out = @map;
    foreach $i (9..$#map) {
      if ($do_remap && length $map[$i] > $max_name_length) {
        print "##Name_change=Field $i has name longer than $max_name_length; $map[$i] changed to Remap_$i\n";
        $out[$i] = "Remap_$i";
      }
    }
    print join ("\t", @out), "\n";
    last;
  }
}
if (!$do_gcov) {
  if ($number_chromosomes) {
    while (<>) {
      chomp;
      @f = split /\t/, $_;
      if (!defined $chrom->{$f[0]}) {
        $chrom->{$f[0]} = $chrom_counter;
        $chrom_counter++;
      }
      $f[0] = $chrom->{$f[0]};
      print join ("\t", @f), "\n";
    }
  } else {
    while (<>) { print; }
  }
  exit;
}
@in = <>;
my $out;
my @f;
foreach $i (0..$#in) {
  chomp $in[$i];
  @f = split /\t/, $in[$i];
  if (!defined $chrom->{$f[0]}) {
    $chrom->{$f[0]} = $chrom_counter;
    $chrom_counter++;
  }
  if ($number_chromosomes) {
    $f[0] = $chrom->{$f[0]};
  }
  @{$out->[$i]} = @f;
}
foreach $i (9..$#map) {
  next if length($target) && $map[$i] ne $target;
  if (-f "$data_dir/$map[$i].gcovs.gz") {
    open F, "zcat $data_dir/$map[$i].gcovs.gz |";
  } elsif (-f "$data_dir/$map[$i].gcov.gz") {
    open F, "zcat $data_dir/$map[$i].gcov.gz |";
  } elsif (-f "$data_dir/$map[$i].gcovs") {
    open F, "$data_dir/$map[$i].gcovs";
  } elsif (-f "$data_dir/$map[$i].gcov") {
    open F, "$data_dir/$map[$i].gcov";
  } else {
    print STDERR "Can't find $map[$i].gcovs file in $data_dir, skipping...\n";
    next;
  }
  $verbose && print STDERR "Starting $i...\n";
  $cov = [];
  while (<F>) {
    next if /^#/;
    next if /^$/;
    chomp;
    @g = split /\t/, $_;
    if ($g[2] >= $min_cov) {
      $cov->[$g[1]] = $g[2];
    } else {
      $cov->[$g[1]] = 0;
    }
  }
  close F;
  foreach $j (0..$#$out) {
    if ($out->[$j]->[$i] eq "./.:.:." && $cov->[$out->[$j]->[1]]) {
      $out->[$j]->[$i] = "0/0:$cov->[$out->[$j]->[1]]:PASS";
    }
  }
  $verbose && print STDERR "Finished $i\n";
}
foreach $i (0..$#$out) {
  print join ("\t", @{$out->[$i]}), "\n";
}
if ($number_chromosomes) {
  print STDERR "Chromosome\tNumber\n";
  foreach $i (sort {$chrom->{$a} <=> $chrom->{$b}} keys %$chrom) {
    print STDERR join ("\t", $i, $chrom->{$i}), "\n";
  }
}
