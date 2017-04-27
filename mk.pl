#!/usr/bin/perl -w
#
# do a McDonald-Kreitman test
# tabulate synonymous and nonsynonymous fixed and polymorphic sites
# need 2 sequence files, assume fasta and already aligned
# also assume we are already in frame
# reference: McDonald JH and Kreitman M. (1991) Nature 351:652-4
#
use slchen;

if (!$ARGV[0] || !-f $ARGV[0] || !$ARGV[1] || !-f $ARGV[1]) {
  print "Usage: $0 <file1> <file2>\n";
  exit;
}

my ($group1, $group2) = @ARGV[0,1];

my ($seq1, $nt1) = get_seq($group1);
my ($seq2, $nt2) = get_seq($group2);
my ($aa1, $aa2);
my ($length, $aalength);
my ($fixed_s, $fixed_n) = (0, 0);
my ($poly_s, $poly_n) = (0, 0);
my $name;
my $i;

# check lengths
$length = 0;
foreach $name (keys %$seq1) {
  $length = length $seq1->{$name} if !$length;
  die "Wrong length for $name, expected $length but got ", length $seq1->{$name}, "\n" if $length != length $seq1->{$name};
}
foreach $name (keys %$seq2) {
  $length = length $seq2->{$name} if !$length;
  die "Wrong length for $name, expected $length but got ", length $seq2->{$name}, "\n" if $length != length $seq2->{$name};
}

# translate to aa
$aalength = 0;
foreach $name (keys %$seq1) {
  @{$aa1->{$name}} = ();
  for ($i = 0; $i < length $seq1->{$name}; $i += 3) {
    push @{$aa1->{$name}}, aa(substr($seq1->{$name}, $i, 3));
  }
  $aalength = length $aa1->{$name} if !$aalength;
}

foreach $name (keys %$seq2) {
  @{$aa2->{$name}} = ();
  for ($i = 0; $i < length $seq2->{$name}; $i += 3) {
    push @{$aa2->{$name}}, aa(substr($seq2->{$name}, $i, 3));
  }
}

# go through sequence
foreach $i (0..$length-1) {
  $j = int ($i/3);	# amino acid index
  # check if it's all the same
  $data = {};
  foreach $name (keys %$seq1) {
    $data->{1}->{uc $nt1->{$name}->[$i]} = 1;
  }
  foreach $name (keys %$seq2) {
    $data->{2}->{uc $nt2->{$name}->[$i]} = 1;
  }

  $data->{COUNT}->{1} = 0;
  foreach $base (keys %{$data->{1}}) {
    $data->{COUNT}->{1} += $data->{1}->{$base};
  }
  $data->{COUNT}->{2} = 0;
  foreach $base (keys %{$data->{2}}) {
    $data->{COUNT}->{2} += $data->{2}->{$base};
  }
  $data->{COUNT}->{ALL} = 0;
  foreach $base qw(G A T C) {
    $data->{COUNT}->{ALL}++ if $data->{1}->{$base} || $data->{2}->{$base};
  }

  # short circuit if no polymorphism anywhere
  next if $data->{COUNT}->{ALL} <= 1;

  $consensus_aa = "";
  $nonsyn = 0;
  foreach $name (keys %$aa1) {
    last if $nonsyn;
    $consensus_aa = $aa1->{$name}->[$j] if $consensus_aa eq "";
    if ($aa1->{$name}->[$j] ne $consensus_aa) {
      $nonsyn = 1;
    }
  }
  foreach $name (keys %$aa2) {
    last if $nonsyn;
    if ($aa2->{$name}->[$j] ne $consensus_aa) {
      $nonsyn = 1;
    }
  }

  if ($nonsyn) {
    if ($data->{COUNT}->{1} == 1 && $data->{COUNT}->{2} == 1) {
      $fixed_n++;
    } else {
      $poly_n++;
    }
  } else {
    if ($data->{COUNT}->{1} == 1 && $data->{COUNT}->{2} == 1) {
      $fixed_s++;
    } else {
      $poly_s++;
    }
  }
}

print "Fixed:\n    syn: $fixed_s\n    non: $fixed_n\n";
print "Polymorphic:\n    syn: $poly_s\n    non: $poly_n\n";

sub get_seq {
  # take a file name, get the sequences into a hash ref
  my ($file) = @_;
  my $seq;
  my $nt;
  my $name;
  undef $name;
  open F, $file;
  while (<F>) {
    chomp;
    next if /^#/;
    next if /^$/;
    if (/^>/) {
      $name = $_;
    } else {
      $seq->{$name} = $_ if defined $name;
      @{$nt->{$name}} = split //, $_ if defined $name;
      undef $name;
    }
  }
  close F;
  return $seq, $nt;
}
