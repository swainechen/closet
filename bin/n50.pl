#!/usr/bin/perl -w
#
# calculate N50 and N50 length from fasta sequences
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$output_table = 0;
$contig_cutoff = 1000;

GetOptions (
  'table!' => \$output_table,
  'contig_cutoff=i' => \$contig_cutoff
);

$total = 0;
$contig_num = 0;
$large_contig_num = 0;
@length = ();
$file = "<stdin>";
$file = $ARGV[0] if defined $ARGV[0] && -f $ARGV[0];
$name = "__n50_notset";
$subtotal = 0;
while (<>) {
  chomp;
  if (/^>(.*)/) {
    if ($name ne "__n50_notset") {
      $contig_num++;
      $large_contig_num++ if $subtotal >= $contig_cutoff;
      $total += $subtotal;
      push @length, $subtotal;
      $subtotal = 0;
    }
    $name = $1;
    next;
  }
  if ($name eq "__n50_notset") {
    $name = "__n50_firstcontig";
  }
  $subtotal += length ($_);
}
if ($name ne "__n50_notset") {
  $contig_num++;
  $large_contig_num++ if $subtotal >= $contig_cutoff;
  $total += $subtotal;
  push @length, $subtotal;
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
