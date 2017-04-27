#!/usr/bin/perl -w
#
# take fasta sequence, probably protein
# filter out redundant sequences
# align it
# spit out entropy at each position
# also spit out the percent consensus, i.e. percent that the most frequent one accounts for
#
use slchen;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $do_align = 1;

GetOptions (
  'align!' => \$do_align	# whether to do an alignment, i.e. we don't if
				# input is aligned already
);

my @a = <>;
my $seq = slchen::fasta2hash(@a);
$seq = slchen::filter_redundant($seq);
@a = slchen::hash2fasta($seq);
@a = slchen::align(@a) if $do_align;
print join ("\n", @a), "\n";
$seq = slchen::fasta2hash(@a);

my $i;
my $name;
my $length;
my $seqarray;
my $top_frequency;
my $top_allele;
my %count;
my $entropy;
foreach $name (keys %$seq) {
  $length = length $seq->{$name};
  @{$seqarray->{$name}} = split //, $seq->{$name};
}

print join ("\t", "Position", "Entropy", "Top Freq", "Alleles..."), "\n";
foreach $i (0..$length-1) {
  @a = ();
  foreach $name (keys %$seqarray) {
    push @a, $seqarray->{$name}->[$i] if $seqarray->{$name}->[$i] ne '-' && $seqarray->{$name}->[$i] ne ' ';
  }
  $entropy = slchen::shannon(@a);
  $top_frequency = 0;
  $top_allele = "";
  %count = ();
  @a = sort @a;
  foreach $a (@a) {
    $count{$a}++;
  }
  foreach $a (keys %count) {
    $top_frequency = $count{$a} if $count{$a} > $top_frequency;
  }
  if (scalar @a) {
    $top_frequency /= scalar(@a);
  }
  print join ("\t", $i+1, $entropy, $top_frequency, sort {$count{$a} <=> $count{$b}} keys %count), "\n";
}
