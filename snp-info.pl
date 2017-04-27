#!/usr/bin/perl -w
#
# take in a position and a nucleotide change
# spit out information
# use GERMS library
# generally will need all the NCBI Refseq information - ptt, faa, fna
#
use strict;
use GERMS;
use Orgmap;
&read_orgmap;

my @f;
my $chrom;
my $pos;
my $ref;
my $alt;
my $data;
my @out;
my $accession;
my $gbkfile;
my @a;
my $fna_header;
my ($ptt_ref, $sequence) = GERMS::ns_setup($Orgmap::pttfile, $Orgmap::fnafile);
$gbkfile = $Orgmap::fnafile;
$gbkfile =~ s/\.fna$/.gbk/;
open F, $Orgmap::fnafile;
$fna_header = <F>;
if ($fna_header =~ s/^>//) {
  @f = split /\s+/, $fna_header;
  $fna_header = $f[0];
} else {
  $fna_header = "";
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
  while (<>) {
    if (/^$/) { print; next; }
    if (/^#/) { print; next; }
    chomp;
    # VCF input
    @f = split /\t/, $_;
    $chrom = $f[0];
    $pos = $f[1];
    $ref = uc $f[3];
    $alt = uc $f[4];
    undef $data;
    if ($chrom eq $accession || $chrom eq $fna_header || $accession eq "") {
      if ($ref eq substr($sequence, $pos-1, length($ref))) {
        $data = GERMS::ns($pos, $alt, $ptt_ref, $sequence);
      } else {
        print STDERR "Error, $chrom $pos is noted as $ref but I see ", substr($sequence, $pos-1, length($ref)), " in $Orgmap::fnafile...skipping\n";
      }
    }
    @out = @f;
    if (defined $data) {
      if ($data->{INTERGENIC}) {
        if ($f[1] eq '.') {
          push @out, "Insertion";
        } elsif ($alt eq '.') {
          push @out, "Deletion";
        } else {
          push @out, "Substitution";
        }
        @a = `echo $pos | nearestorf.pl $orgcode -ptt $Orgmap::pttfile -fna $Orgmap::fnafile | tail -n 2 | cut -f5- -d' '`;
        if (defined $a[0] && defined $a[1]) {
          chomp $a[0];
          chomp $a[1];
        } else {
          $a[0] = "";
          $a[1] = "";
        }
        push @out, "Intergenic", $a[0], $a[1];
      } else {
        if ($f[1] eq '.') {
          push @out, "Insertion";
        } elsif ($alt eq '.') {
          push @out, "Deletion";
        } else {
          push @out, "Substitution";
          if ($data->{SYNONYMOUS}) {
            push @out, "Synonymous";
            push @out, "";
          } else {
            push @out, "Nonsynonymous";
            push @out, $data->{ORIGINALAA} . $data->{AAPOSITION} . $data->{NEWAA};
          }
        }
        push @out, $data->{GID}, $data->{SYSTEMATIC}, $ptt_ref->{$data->{GID}}->{TRIVIAL}, $ptt_ref->{$data->{GID}}->{ANNOTATION};
      }
      print join ("\t", @out), "\n"; 
    }
  }
}
