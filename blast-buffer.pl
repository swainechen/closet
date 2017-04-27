#!/usr/bin/perl -w
#
# only feed 200 sequences at a time to blast
#
use Getopt::Long;
use Orgmap;
&Getopt::Long::Configure("pass_through");

if (!@ARGV) {
  print "Usage: $0 <orgcode> [ -F T|F ] [ -m # ] [ -b # ] [ -v # ] [ -e <float> ] [ -n <number of fasta lines to feed blast at a time> ] <fasta sequence file>\n";
  exit;
}

if ($ARGV[0] eq 'nt') {
  $db = "/usr/local/lib/Genomes/nr/nt";
  shift;
} else {
  use Orgmap qw(:DEFAULT $fileprefix);
  &read_orgmap;
  $db = $fileprefix.".fna";
}

$eValue = 0.0000000001;
$F = 'F';		# no filtering if this is 'F'
$m = 8;
$b = 1;
$v = 0;
$n = 200;		# number of sequences to give at a time

GetOptions (
  'db=s' => \$db,
  'F=s' => \$F,
  'b=i' => \$b,
  'v=i' => \$v,
  'm=i' => \$m,
  'e=f' => \$eValue,
  'n=i' => \$n
);

@out = ();
$i = 0;
while (<>) {
  $i++ if /^>/;
  if ($i >= $n) {
    open F, "| blastn $orgcode -F $F -b $b -v $v -m $m -e $eValue";
    print F @out;
    close F;
    @out = ();
    $i = 1;
  }
  push @out, $_;
}
if ($i) {
  open F, "| blastn $orgcode -F $F -b $b -v $v -m $m -e $eValue";
  print F @out;
  close F;
}
