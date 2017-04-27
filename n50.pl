#!/usr/bin/perl -w
#
# calculate N50 and N50 length from fasta sequences
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$output_table = 0;

GetOptions (
  'table!' => \$output_table
);

$total = 0;
$contig_num = 0;
$large_contig_num = 0;
$contig_cutoff = 1000;
@length = ();
$file = "<stdin>";
$file = $ARGV[0] if defined $ARGV[0] && -f $ARGV[0];
while (<>) {
  chomp;
  next if /^>/;
  $total += length ($_);
  $contig_num++;
  $large_contig_num++ if length($_) >= $contig_cutoff;
  push @length, length ($_);
}

@length = sort { $a <=> $b } @length;
$sub = 0;
foreach $i (reverse 0..$#length) {
  $sub += $length[$i];
  if ($sub >= $total/2) {
    if ($output_table) {
      print join ("\t", $file, $total, $#length - $i + 1, $length[$i], $contig_num, $large_contig_num), "\n";
    } else {
      print "total contig size: $total\n";
      print "N50: ", $#length - $i + 1, "\n";
      print "N50 length: $length[$i]\n";
      print "Largest contig: $length[$#length]\n";
      print "Number of contigs: $contig_num\n";
      print "Number of contigs > $contig_cutoff: $large_contig_num\n";
    }
    last;
  }
}
