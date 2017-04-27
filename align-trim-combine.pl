#!/usr/bin/perl -w
#
# take fasta sequence
# align using clustal
# rename sequences and take out identical sequences
# deal with blank fasta names too
#
use slchen;
use File::Temp;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my @seq;	# fasta
my $seq;	# hash ref
my $name_hash;
my $tempdir = File::Temp::tempdir( CLEANUP => 1 );
my $type = "protein";	# where to do alignment, in dna or protein
my $z;		# holds translation errors from paln2daln
my @temp;
my $special = "";	# holds special name to be put first and to reference numbers off of
my $i;
my $special_seq;
my @spec;	# holds special sequence as array
my @index;
my (@left, @right);
my $j;
my $mode_l;
my $mode_r;
my $codon;
my $name;
my $seq_ref;
my $sequence;
my $name_map;

GetOptions (
  'type=s' => \$type,
  'special=s' => \$special
);

while (<>) {
  chomp;
  next if /^#/; next if /^$/;
  push @seq, $_;
}

# this back and forth should get rid of blank fasta lines
$seq = slchen::fasta2hash(@seq);
@seq = slchen::hash2fasta($seq);

# align
$name_hash = slchen::to_phylip_names(\@seq);
open T, "| translate.pl -noseq -nosix -notri > $tempdir/translate.tmp";
print T join ("\n", @seq), "\n";
close T;
open T, "$tempdir/translate.tmp";
@pseq = ();
while (<T>) {
  chomp; push @pseq, $_;
}
unlink "$tempdir/translate.tmp";

if ($type eq 'dna') {
  # for dna alignment
  @seq = slchen::dnaalign(@seq);
} else {
  # for protein alignment, then convert back to dna
  @pseq = slchen::align(@pseq);
  $z = 0;
  @temp = slchen::paln2daln(\@pseq, \@seq, \$z);
  @seq = @temp;
}

&slchen::from_phylip_names(\@seq, $name_hash);

# find special and trim
if (length($special)) {
  foreach $i (0..$#seq) {
    if ($seq[$i] =~ /^>$special/i) {
      $special_seq = $seq[$i+1];
      last;
    }
  }
  @spec = split //, $special_seq if length $special_seq;
}
if (!@spec) {
  @spec = split //, $seq[1];
}
@index = ();
foreach $i (0..$#spec) {
  push @index, $i if $spec[$i] ne '-';
}

# trim to special sequence
# collect information on gaps at left and right end
foreach $i (0..$#seq) {
  if ($seq[$i] =~ /^>/) { next; }
  @temp = split //, $seq[$i];
  $seq[$i] = join "", @temp[@index];
  foreach $j (0..$#index) {
    next if $temp[$index[$j]] eq '-';
    push @left, $j;
    last;
  }
  foreach $j (reverse 0..$#index) {
    next if $temp [$index[$j]] eq '-';
    push @right, $j;
    last;
  }
}

# use the mode of the left and right gaps and trim down to the nearest codon
$mode_l = slchen::array_mode(@left);
$mode_r = slchen::array_mode(@right);
# actually, hope that sequences have been selected with min homology length
# take minimal sequence with no gaps on either side
# we want 0 to go to 0, but 1 and 2 to go to 3
@left = slchen::sortu @left;
@right = slchen::sortu @right;
$mode_l = pop @left;
$mode_r = shift @right;
$mode_l += 2;
$mode_l -= $mode_l % 3;
# this time we want 30 and 31 to go to 29 but 32 to go to 32
$mode_r += 1;
$mode_r -= $mode_r % 3;
$mode_r -= 1;

foreach $i (0..$#seq) {
  next if ($seq[$i] =~ /^>/);
  $seq[$i] = substr($seq[$i], $mode_l, $mode_r - $mode_l + 1);
  # blank out stop codons
  $j = 0;
  while ($j < length $seq[$i]) {
    $codon = substr ($seq[$i], $j, 3);
    if ($codon =~ /taa/i || $codon =~ /tag/i || $codon =~ /tga/i) {
      substr ($seq[$i], $j, 3) = '---';
    }
    $j += 3;
  }
}

# make new sequence names for unique sequences
# spit out table to map back from new sequence names to old names
foreach $i (0..$#seq) {
  if ($seq[$i] =~ /^>/) {
    $name = $seq[$i];
    $name =~ s/^>//;
  } else {
    push @{$seq_ref->{$seq[$i]}}, $name;
  }
}

$name_map = 'name_map.txt';
open N, ">$name_map";
print N "# All sequences trimmed to nucleotides ", $mode_l+1, " to ", $mode_r+1, " of $special sequence\n" if length($special);
print N "# Sequence\tOriginal Name\n";
$i = 0;
foreach $sequence (keys %$seq_ref) {
  $i++;
  print ">Seq_$i\n$sequence\n";
  foreach $j (0..$#{$seq_ref->{$sequence}}) {
    print N "$i\t", $seq_ref->{$sequence}->[$j], "\n";
  }
}
