#!/usr/bin/perl -w
#
# given an orgcode
# and a mutation (like G352A)
# print whether it's synonymous or nonsynonymous or intergenic
#
use slchen;
use Orgmap;
&read_orgmap;
&Orgmap::read_sequence;

my $ptt;
my ($s, $e);
open PTT, $Orgmap::pttfile;
while (<PTT>) {
  if (/^\s*(\d+)\.\.(\d+)\s/) {
    ($s, $e) = ($1, $2);
    @f = split /\t/, $_;
    $strand = $f[1];
    $ptt->{$s}->{END} = $e;
    $ptt->{$s}->{STRAND} = $strand;
  }
}

while (<>) {
  if (/^#/ || /^>/) {
    print;
    next;
  }
  chomp;
  if (uc $_ =~ /([GATC])(\d+)([GATC])/) {
    $orig = $1;
    $pos = $2;
    $mut = $3;
    next if $orig eq $mut;
    if (uc substr($Orgmap::sequence, $pos-1, 1) eq $orig) {
      ($orig_codon, $new_codon) = codon($pos, $mut, $ptt);
      if ($orig_codon eq "INTERGENIC") {
        print "$_\tI\n";
      } else {
        if (slchen::aa($orig_codon) eq slchen::aa($new_codon)) {
          print "$_\tS\n";
        } else {
          print "$_\tN\n";
        }
      }
    }
  }
}

sub codon {
  # given a position, mutation, and ptt ref
  # return the appropriate codon and its mutation
  my ($p, $m, $ptt) = @_;
  my $i;
  my $codon_num;
  my $orig;
  my $new;
  my $s;
  foreach $s (sort {$a <=> $b} keys %$ptt) {
    last if $s > $p;
    if ($p <= $ptt->{$s}->{END}) {
      if ($ptt->{$s}->{STRAND} eq "-") {
        # muck with rev complement
        $i = $ptt->{$s}->{END};
        while ($i - 3 >= $p) {
          $i -= 3;
        }
        $orig = substr($Orgmap::sequence, $i - 3, 3);
        $new = $orig;
        substr($new, $p - $i + 2, 1) = $m;
        $orig =~ tr/GATCgatc/CTAGctag/;
        $new =~ tr/GATCgatc/CTAGctag/;
        $orig = reverse $orig;
        $new = reverse $new;
      } else {
        # we're on the positive strand, easier
        $i = $s;
        while ($i + 3 <= $p) {
          $i += 3;
        }
        $orig = substr($Orgmap::sequence, $i - 1, 3);
        $new = $orig;
        substr($new, $p - $i, 1) = $m;
      }
      return ($orig, $new);
    }
  }
  return ("INTERGENIC", "INTERGENIC");
}
