#!/usr/bin/perl -w
#
# take in a position and a nucleotide change
# spit out information
# use GERMS library
# generally will need all the NCBI Refseq information - ptt, faa, fna
#
use strict;
use lib '/home/slchen/bin/germs/lib';
use GERMS;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
my $gbkfile = "";
my $ignore_accession = 1;
GetOptions (
  'gbkfile=s' => \$gbkfile,
  'ignore_sequence_name!' => $ignore_accession
);
if (!-f $gbkfile) {
   print "Usage: $0 -gbkfile <gbk file> vcffile\n";
   print "vcf file can come on stdin or as a filename\n";
   exit;
}

my @f;
my $chrom;
my $pos;
my $ref;
my $alt;
my $data;
my @out;
my $accession;
my @a;
my ($ptt_ref, $sequence, $gbk_ref) = GERMS::ns_setup_gbk($gbkfile);
my $sorted_ptt;
my $topology = 0;
my $locus;
my @alts;	# in case of multiple alleles at the same position

@$sorted_ptt = ();
foreach $locus (sort { $ptt_ref->{$a}->{START} <=> $ptt_ref->{$b}->{START} ||
                     $ptt_ref->{$a}->{END} <=> $ptt_ref->{$b}->{END} }
                keys %$ptt_ref) {
    $ptt_ref->{$locus}->{LOCUS} = $locus;
    push @$sorted_ptt, $ptt_ref->{$locus};
}
foreach $locus (keys %{$gbk_ref->{__HEADER__}}) {
  $topology = 1 if $gbk_ref->{__HEADER__}->{$locus} =~ /(^|\s)linear($|\s)/;
}

if (defined $ptt_ref && defined $sequence) {
  # get accession / version
  open G, $gbkfile;
  while (<G>) {
    if (/^VERSION/) {
      chomp;
      @f = split /\s+/, $_;
      if (defined $f[1]) {
        $accession = $f[1];
      } else {
        $accession = "";
      }
      last;
    }
  }
#  print join ("\t", "# Chromosome", "Position", "VCFID", "Ref", "Alt", "Mutation", "Location", "Impact/Left", "GID/Right", "Systematic", "Trivial", "Annotation(s)"), "\n";
  while (<>) {
    if (/^$/) { print; next; }
    if (/^##/) { print; next; }
    if (/^#CHROM/) { chomp; print join ("\t", $_, "Mutation", "Location", "Impact/Left", "GID/Right", "Systematic", "Trivial", "Annotation(s)"), "\n"; next; }
    chomp;
    # VCF input
    @f = split /\t/, $_;
    $chrom = $f[0];
    $pos = $f[1];
    $ref = uc $f[3];
    $alt = uc $f[4];
    @alts = split /,/, $alt;
    foreach $alt (@alts) {
      undef $data;
      if ($chrom eq $accession || $accession eq "" || $ignore_accession) {
        if ($ref eq substr($sequence, $pos-1, length($ref))) {
          $data = GERMS::ns($pos, $alt, $ptt_ref, $sequence);
        } else {
          print STDERR "Error, $chrom $pos is noted as $ref but I see ", substr($sequence, $pos-1, length($ref)), " in $gbkfile...skipping\n";
        }
      }
#      @out = @f[0..3];
      @out = @f;
      $out[4] = $alt;
      if (defined $data) {
        if ($data->{INTERGENIC}) {
          if ($f[1] eq '.') {
            push @out, "Insertion";
          } elsif ($alt eq '.') {
            push @out, "Deletion";
          } else {
            push @out, "Substitution";
          }
          @a = nearestorf($pos, $sorted_ptt, $topology);
          if (defined $a[0] && defined $a[1]) {
            chomp $a[0];
            chomp $a[1];
          } else {
            $a[0] = "";
            $a[1] = "";
          }
          push @out, "Intergenic", $a[0], $a[1], $ptt_ref->{$a[0]}->{TRIVIAL}, $ptt_ref->{$a[1]}->{TRIVIAL}, "$ptt_ref->{$a[0]}->{ANNOTATION}; $ptt_ref->{$a[1]}->{ANNOTATION}";
        } else {
          if ($f[1] eq '.') {
            push @out, "Insertion";
          } elsif ($alt eq '.') {
            push @out, "Deletion";
          } else {
            push @out, "Substitution";
            if ($ptt_ref->{$data->{GID}}->{TYPE} eq "CDS") {
              if ($ptt_ref->{$data->{GID}}->{COMPLETE}) {
                if ($data->{SYNONYMOUS}) {
                  push @out, "Synonymous";
                  push @out, "-";
                } else {
                  push @out, "Nonsynonymous";
                  push @out, $data->{ORIGINALAA} . $data->{AAPOSITION} . $data->{NEWAA};
                }
              } else {
                push @out, "Incomplete CDS", "Unknown";
              }
            } else {
              push @out, $ptt_ref->{$data->{GID}}->{TYPE}, "-";
            }
            push @out, $data->{GID}, $data->{SYSTEMATIC}, $ptt_ref->{$data->{GID}}->{TRIVIAL}, $ptt_ref->{$data->{GID}}->{ANNOTATION};
          }
        }
        print join ("\t", @out), "\n"; 
      }
    }
  }
}

sub nearestorf {
  my ($n, $ptt, $topology) = @_;

  my $low = 0;
  my $high = scalar (@$ptt) - 1;
  my ($lowest_s, $lowest_e);
  my ($highest_s, $highest_e);
  my ($s, $e);
  my $now;

# first boundary conditions
  ($lowest_s, $lowest_e) = ($ptt->[$low]->{START}, $ptt->[$low]->{END});
  ($highest_s, $highest_e) = ($ptt->[$high]->{START}, $ptt->[$high]->{END});
  if ((($n < $lowest_s) && ($n < $lowest_e)) ||
     (($n > $highest_s) && ($n > $highest_e))) {
    if ($topology == 0) {       # circular topology, do some wrap-around
      $high = $low;
      $low = scalar (@$ptt) - 1;
      $now = $high;
    } else {                    # linear topology, no wrap-around
      if (($n < $lowest_s) && ($n < $lowest_e)) {
        $high = $low;
        $low = -1;
        $now = $high;
      } else {
        $low = $high;
        $high = -1;
        $now = $low;
      }
    }
  }

# Binary search.  If inside an orf set $high = $low at that orf.  If
# not in an orf then $high = $low + 1, $high and $low should be nearest orfs.
  while (($high - $low > 1) && ($low >= 0)) {
    $now = int (($low + $high)/2);
    ($s, $e) = ($ptt->[$now]->{START}, $ptt->[$now]->{END});
    if (($n < $s) && ($n < $e)) {               # we're not far enough
      $high = $now;
    }
    elsif (($n > $s) && ($n > $e)) {            # we're past it
      $low = $now;
    }
    elsif ((($n >= $s) && ($n <= $e)) || (($n <= $s) && ($n >= $e))) {
      $low = $now;                              # we're inside the orf
      $high = $now;
    }
  }
  # check if we missed checking one
  if (($high != $low) && ($high >= 0) && ($low >= 0)) {
    if ($low == $now) {
      ($s, $e) = ($ptt->[$high]->{START}, $ptt->[$high]->{END});
    } else {
      ($s, $e) = ($ptt->[$low]->{START}, $ptt->[$low]->{END});
    }
    if ((($n >= $s) && ($n <= $e)) || (($n <= $s) && ($n >= $e))) {
      if ($low == $now) { $low = $high; }
      else { $high = $low; }
    }
  }
  if ($high == $low) {
    return($ptt->[$high]->{LOCUS}, "__inside__");
  } else {
    return($ptt->[$low]->{LOCUS}, $ptt->[$high]->{LOCUS});
  }
}

sub wraparound {
  my ($i, $max) = @_;
  # assume have an array with $max+1 elements ($max is the largest index)
  # perl already takes care of negative $i
  # if $i is greater than $max, we're going to wrap that around to the beginning
  # just to be complete we'll take care of negative $i also
  if ($i < 0) {
    return ($max + 1 + $i);     # want to start at $max if $i = -1
  }
  if ($i > $max) {
    return ($i - 1 - $max);     # want to start at 0 if $i = $max + 1
  }
  return $i;
}

