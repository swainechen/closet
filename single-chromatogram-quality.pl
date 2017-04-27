#!/usr/bin/perl -w
#
# get high quality sequence using phred from each chromatogram
# import it back into mysql
#
# 050518: add blast information as well
# 060328: mods to use mysql databases
# 080910: mod to take single chromatograms from the command line
#
use slchen;
use File::Temp qw(tempdir);
my $debug = 0;

use Getopt::Long;
&Getopt::Long::Configure("pass_through");
GetOptions (
  'debug=i' => \$debug,
);

# expect file names on command line
if (!@ARGV) {
  print "Usage: $0 <chromatogram1> [ <chromatogram2> [ ... ] ]\n";
  print "High quality sequence for each chromatogram will be output in fasta\n";
  print "format to standard output, with sequence names given by\n";
  print "chromatogram file names\n";
  exit;
}

# make a temp directory
$working_dir = tempdir( CLEANUP => 1 );
mkdir "$working_dir/chromat";

# set up database connection, pull out chromatograms
# give them St. Louis convention filenames
# i.e. just need to get primer direction
# then b = +, g = -
# readname.[bg]anything_else.ab1

foreach $file (@ARGV) {

  next if !-f $file;
  my $i = 0;
  # default type to abi
  $extension = '.abi';
  if ($file =~ /scf$/) {
    $extension = '.scf';
  }
  $filename = "chromatogram.bchromatogram" . $extension;

  system "cp $file $working_dir/chromat/$filename";
  $debug && print ("$file\t$filename\n");
  system "phred $working_dir/chromat/$filename -sa $working_dir/seq -qa $working_dir/qual -trim_fasta -trim_alt \"\" 2>&1 > /dev/null";
  open Q, "$working_dir/qual";
  while (<Q>) {
    if (/^>/) {
      my @f = split;
      while ($f[$#f] =~ /SCF/ || $f[$#f] =~ /trim/ || $f[$#f] =~ /ABI/) { pop @f; }
      $highqual = pop @f;
    }
  }
  close Q;
  $debug && print ("highquality: $highqual\n");

  my $hqseq = "";
  open S, "$working_dir/seq";
  while (<S>) {
    next if /^>/;
    chomp;
    $hqseq .= $_;
  }
  close S;

  # print out fasta format
  print ">$file\n";
  print "$hqseq\n";

  # cleanup
  !$debug && system "rm -f $working_dir/chromat/*";
  $debug && exit;
}
