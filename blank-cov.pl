#!/usr/bin/perl -w
#
#
use slchen;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $gcov;
my $mincov = 10;
GetOptions (
  'gcov=s' => \$gcov,
  'mincov=i' => \$mincov
);

if (!defined $gcov || !-f $gcov) {
  print "Usage: $0 -gcov <gcov.gz file> <reference sequence>\n";
  exit;
}

if ($gcov =~ /gz$/) {
  open G, "zcat $gcov |";
} else {
  open G, $gcov;
}

@a = <>;
$a = slchen::fasta2hash(@a);
$error = {};
foreach $i (keys %$a) {
  @{$edit->{$i}} = split //, $a->{$i};
}

while (<G>) {
  next if /^#/;
  next if /^$/;
  chomp;
  @f = split /\t/, $_;
  next if $f[2] >= $mincov;
  if (!defined $a->{$f[0]} && !defined $error->{$f[0]}) {
    $error->{$f[0]} = 1;
    print STDERR "Error on chromosome $f[0] in gcov, not in sequence\n";
  }
  if ($f[1]-1 > $#{$edit->{$f[0]}}) {
    print STDERR "Error coordinate $f[1] out of range for chromosome $f[0], max expected $#{$edit->{$f[0]}}\n";
  }
  $edit->{$f[0]}->[$f[1]-1] = "-";
}

foreach $i (keys %$a) {
  $a->{$i} = join ("", @{$edit->{$i}});
}
print join ("\n", slchen::hash2fasta($a)), "\n";
