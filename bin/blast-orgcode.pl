#!/usr/bin/perl -w

# blast-orgcode.pl
# This is just a frontend for blastall so that you can do a blastp or
# blastn against a given organism more easily.
# This assumes you only want the top hit for each query sequence
# and uses 2 CPUs by default, no alignments are shown.
# Only allows blastp and blastn searches.

use Getopt::Long;
&Getopt::Long::Configure("pass_through");
use Orgmap qw(:DEFAULT $fileprefix);
&read_orgmap;

$program = '';
$inputFile = '';
$eValue = 100.0;
$F = 'T';		# no filtering if this is 'F'

GetOptions ('input=s' => \$inputFile,
            'program=s' => \$program,
            'F=s' => \$F,
            'e=f' => \$eValue);

if ($inputFile eq '-') {
  @in = <>;
  $tempfile = "blast-orgcode-".rand().".temp";
  open F, ">$tempfile";
  print F @in;
  close F;
  $inputFile = $tempfile;
}

if (!$inputFile || !-e $inputFile || ($program ne 'blastp' && $program ne 'blastn')) {
  print "Usage: blast-orgcode.pl ORG-CODE -program [blastp|blastn] -input [query file | - ] -e [Expectation Value]\n";
  exit(-1);
}

if ($program eq 'blastn') {
#  $db = $fileprefix.".ffn";
  $db = $fileprefix.".fna";
} else {
  $db = $fileprefix.".faa";
}

system "blastall -p $program -d $db -i $inputFile -e $eValue -v 0 -b 1 -a 2 -F $F";

if ($tempfile) {
  unlink $tempfile;
}
