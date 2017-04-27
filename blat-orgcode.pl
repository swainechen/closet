#!/usr/bin/perl -w

# blat-orgcode.pl
# This is just a frontend for blat
# assumes we're only doing DNA to DNA blat

use Getopt::Long;
use File::Temp;
&Getopt::Long::Configure("pass_through");
if ($ARGV[0] eq 'nt') {
  $db = "/usr/local/lib/Genomes/nr/nt";
  shift;
} else {
  use Orgmap qw(:DEFAULT $fileprefix);
  &read_orgmap;
  $db = $fileprefix.".fna";
}

my $tempdir = File::Temp::tempdir( CLEANUP => 1);
my $tempin = "blat-in";
my $tempout = "blat-output";

$dbtype = "dna";
$qtype = "dna";
$minIdentity = 90;
$m = 'psl';
%output = (
  0 => 'blast',
  8 => 'blast8',
  9 => 'blast9'
);

GetOptions (
  'db=s' => \$db,
  'm=i' => \$m,
);

if (!@ARGV) {
  @in = <>;
  open F, ">$tempdir/$tempin";
  print F @in;
  close F;
  @file = ("$tempdir/$tempin");
} else {
  @file = @ARGV;
}

if (defined $output{$m}) {
  $output = $output{$m};
} else {
  $output = $m;
}

@out = ();
foreach $inputFile (@file) {
  system "blat -t=$dbtype -q=$qtype -minIdentity=$minIdentity -out=$output $db $inputFile $tempdir/$tempout";
  open F, "$tempdir/$tempout";
  push @out, <F>;
  close F;
  unlink "$tempdir/$tempout";
}
print @out;

if (-f "$tempdir/$tempin") {
  unlink "$tempdir/$tempin";
}
