#!/usr/bin/perl -w
#
# make a simple one-line graphic to show direction and location of genes
#

if ((!defined ($ARGV[0])) || ($ARGV[0] eq '-h') || ($ARGV[0] eq '--help')) { &PrintUsage; }

use Orgmap qw(:DEFAULT read_sequence $sequence $pttfile $genefield);
&read_orgmap;

use Getopt::Long;
&Getopt::Long::Configure("pass_through");

$position = -1;
$gene = '';
$length = -1;

open PTT, $pttfile;
@ptt = <PTT>;
close PTT;
&read_sequence;
$genome_length = length ($sequence);
if ($length > $genome_length) { print "$length greater than genome length ($genome_length).  Exiting.\n"; exit (1); }
$fullscale = '.' x $genome_length;

while (defined ($in = <>)) {
  $position = -1; $gene = ''; $length = -1;
  chomp $in;
  ($input, $length) = split /\s+/, $in;
  if ((!defined $input) || (!defined $length)) { next; }
  if ($input =~ m/^\d+$/) { $position = $input; } else { $gene = $input; }
  if ($position == -1) { $found = 0; } else { $found = 1; }
  foreach $line (@ptt) {
    if ($line =~ m/^\s*?(\d+)\.\.(\d+)/) { 
      $s = $1; $e = $2;
      @line = split /\t/, $line;
      if ($line[1] =~ m/\+/) { $char = '>'; } else { $char = '<'; }
      substr($fullscale, $s-1, $e-$s+1) = $char x ($e-$s+1);
      if ($gene ne '') {
        if ($line[$genefield] eq $gene) {
          $position = int (($s + $e)/2);
          ++$found;
          $char =~ tr/<>/-+/;
          substr($fullscale, $s-1, $e-$s+1) = $char x ($e-$s+1);
        }
      }
    }
  }
  
  if ($found > 1) { print "Found multiple matches to $gene.  Showing last one.\n"; }
  if ($found >= 1) {
    if (length(int($length/80)) == 1) { $pixel = int($length/8 + 0.5)/10; }
    else { $pixel = substr(int($length/80), 0, 2) . ('0'x(length(int($length/80))-2)); }
    print "$orgcode genome, $length bp centered at $position.  Each char about $pixel bp.\n";
    $display_region = substr($fullscale, int($position - $length/2), $length);
    #@display = split //, $display_region;
    foreach $i (1 .. 80) {
      $index = int(-$length/160 + $i*$length/80 + 0.5);
      if ($i == 40 && $gene eq '') { print '|'; }
      else { print substr($display_region, $index, 1); }
    }
    print "\n";
  }
  elsif ($found == 0) { print "Couldn't find $gene.  Exiting.\n"; exit (2); }
}


sub PrintUsage {
  print "Usage: graphic.pl ORGCODE\n";
  print "Then it takes on STDIN:\n[ position | gene ] length\n";
  print "You should only specify one of position or gene.  You must specify length.\n";
  print "'>' means gene going to the right.  '<' means gene going to the left.\n";
  print "'+' means selected gene, and it goes to the right, '-' it goes to the left.\n";
  print "'.' means intergenic region.\n";
  exit (-1);
}
