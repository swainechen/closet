#!/usr/bin/perl
# modified from https://github.com/bioperl/bioperl-live/blob/master/examples/subsequence.cgi
# take accession number, pull genbank sequence, spit out fasta formatted sequence
#
if (!@ARGV) {
  print "Usage: $0 <accession> [ accession [ ... ] ]\n";
  print "Will grab sequence from GenBank and print fasta sequence to STDOUT\n";
  print "Accession can be accession number (J00522) or GI number (405830)\n";
  exit;
}

use Bio::DB::GenBank;
use slchen;
$gb = new Bio::DB::GenBank;
foreach $acc (@ARGV) {
  if (isint($acc)) {
    eval {
      $seq = $gb->get_Seq_by_gi($acc);
    };
    if ($@) {
      print ">Error: GI number $acc not found\n";
      next;
    }
  } else {
    eval {
      $seq = $gb->get_Seq_by_acc($acc);
    };
    if ($@) {
      print ">Error: Accession $acc not found\n";
      next;
    }
  }
  print ">$acc\n";
  print $seq->seq, "\n";
}
