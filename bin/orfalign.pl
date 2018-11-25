#!/usr/bin/perl -w
#
# take fasta DNA sequence
# check against orgcode
# pick out best reading frame based on blast
# trim sequence to give best reading frame
# align as protein
# trim gaps and stop codons
# print out DNA alignment
#
use slchen;
use File::Temp;
use Getopt::Long;
use Orgmap;
&Getopt::Long::Configure("pass_through");
if ($ARGV[0] eq 'nr') {
  $orgcode = $ARGV[0];
} else {
  &read_orgmap;
}

my $unique = 1;	# whether to only use unique sequences
my $trim_stop = 1;
my $db = "";
GetOptions (
  'unique!' => \$unique,
  'db=s' => \$db,
  'stop!' => \$trim_stop
);

my ($fh, $fn);
my @dna;
my @protein;
my $i;

@dna = <>;
foreach $i (0..$#dna) {
  chomp $dna[$i];
}

my $name;
my $seq_ref;
my $name_hash;
my $key;
if ($unique) {
  foreach $i (0..$#dna) {
    if ($dna[$i] =~ /^>/) {
      $name = $dna[$i];
      $name =~ s/^>//;
    } else {
      push @{$seq_ref->{$dna[$i]}}, $name;
    }
  }
  undef $name_hash;
  undef @dna;
} else {
  $name_hash = to_phylip_names(\@dna);
  undef $seq_ref;
}

-f "name-map.txt" && die "name-map.txt exists...\n";
open N, ">name-map.txt";
if (defined $name_hash) {
  foreach $i (sort {$a<=>$b} keys %$name_hash) {
    print N "$i\t", $name_hash->{$i}, "\n";
  }
} else {
  $i = 0;
  @dna = ();
  foreach $key (keys %$seq_ref) {
    $i++;
    push @dna, ">$i", $key;
    foreach $j (0..$#{$seq_ref->{$key}}) {
      print N "Seq_$i\t", $seq_ref->{$key}->[$j], "\n";
    }
  }
}
close N;

# should have sequences in @dna
# make a $dna_ref, then clobber @dna, then remake it with trimmed sequneces
# get frames using blast
# trim sequence at front and end appropriately
# make $dna_ref first
my $dna_ref = {};
foreach $i (0..$#dna) {
  if ($dna[$i] =~ /^>/) {
    $name = $dna[$i];
    $name =~ s/^>//;
  } else {
    $dna_ref->{$name} = $dna[$i];
  }
}
my %frame = ();
my @blast;
my $tempseq;
my $end_trim;
@dna = ();
foreach $i (keys %$dna_ref) {
  $tempseq = $dna_ref->{$i};
  if (length $db) {
    @blast = `echo $tempseq | blastx $orgcode -m 0 -db $db | grep 'Frame'`;
  } else {
    @blast = `echo $tempseq | blastx $orgcode -m 0 | grep 'Frame'`;
  }
  if (@blast) {
    $blast[0] =~ /Frame\s+=\s+([-+]\d)/;
    $frame{$i} = $1;
  } else {
    $frame{$i} = 1;
  }

  # trim sequence, reverse if needed
  if ($frame{$i} < 0) {
    $dna_ref->{$i} = Orgmap::revcomp($name_hash->{$i});
  }
  substr ($dna_ref->{$i}, 0, abs($frame{$i}) - 1) = "";
  $end_trim = length($dna_ref->{$i});
  $end_trim = $end_trim - 3*(int($end_trim/3));
  substr ($dna_ref->{$i}, -$end_trim) = "" if $end_trim;

  push @dna, ">$i", $dna_ref->{$i};
}

($fh, $fn) = File::Temp::tempfile( UNLINK => 1 );
open T, "| translate.pl -nosix -notri -noseq > $fn";
print T join ("\n", @dna), "\n";
close T;
open T, $fn;
@protein = <T>;
close T;
foreach $i (0..$#protein) {
  chomp $protein[$i];
}
@protein = align(@protein);
@dna = paln2daln(\@protein, \@dna, undef);
@dna = align_stripgaps(@dna);

# take out ambiguous codon sequences
foreach $i (0..$#dna) {
  next if $dna[$i] =~ /^>/;
  for ($j = length($dna[$i]) - 3; $j >= 0; $j -= 3) {
    if (substr($dna[$i], $j, 3) =~ /[^GATCgatc]/) {
      foreach $k (0..$#dna) {
        next if $dna[$k] =~ /^>/;
        substr($dna[$k], $j, 3) = "";
      }
    }
  }
}

#&from_phylip_names(\@dna, $name_hash);

print join ("\n", @dna), "\n";
