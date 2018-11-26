#!/usr/bin/perl
#
use warnings;
use strict;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
use File::Temp;
use LWP::Simple;

my $help = 0;
my $proks = "";
my $max_decoys = 100;
my $output_fna = "Gamma_plasmids.fna";
my $output_meta = "Gamma_plasmids_meta.txt";
my $tempdir = File::Temp::tempdir( CLEANUP => 1 );
my $line = 0;
my @header = ();
my $species = "";
my $species_decoy = ();
my $plasmid_included = ();
my $printed_header = ();
my $outputted_header;
my $secondary_header;
my @seq;
my $seq;
my $fna;
my $url;
my $name;
my $type;
my $acc1;
my $acc2;
my $length;
my $fields;
my @f;
my @g;
my @h;
my $i;
my @out;

GetOptions (
  'help' => \$help,
  'proks=s' => \$proks,
  'decoy=i' => \$max_decoys
);

if ($help) {
  print "Usage: $0 [ -proks <downloaded prokaryotes.txt file> ] [ -decoy <int> ]\n";
  print "Will download get plasmids indicated in the following files:\n";
  print "  ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt\n";
  print "  ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/plasmids.txt\n";
  print "Will include some non-plasmid contigs/chromosomes, maximum number per species set by -decoy (default 100).\n";
  print "Set -decoy 0 to not include them.\n";
  exit;
}

if (-f $output_fna) {
  die "$output_fna already exists, refusing to overwrite...\n";
}

my $return;
if ($proks eq "" || ! -f $proks) {
  $proks = "$tempdir/prokaryotes.txt";
  $return = getstore("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt", $proks);
}

open O, ">$output_fna";
open M, ">$output_meta";
print M join ("\t", qw(Accession Accession2 Length Type Name Species Strain AssemblyAccession Status FTP OrganismName)), "\n";

# ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt first
# expecting these in the header line - * are what we want
#  *0 - Organism/Name
#  1 - TaxID
#  2 - BioProject Accession
#  3 - BioProject ID
#  4 - Group
#  *5 - SubGroup
#  6 - Size (Mb)
#  7 - GC%
#  *8 - Replicons
#  9 - WGS
#  10 - Scaffolds
#  11 - Genes
#  12 - Proteins
#  13 - Release Date
#  14 - Modify Date
#  *15 - Status
#  16 - Center
#  17 - BioSample Accession
#  *18 - Assembly Accession
#  19 - Reference
#  *20 - FTP Path
#  21 - Pubmed ID
#  22 - Strain
#
#open P, "head -n 20 $proks |";
open P, $proks;

while (<P>) {
  $line++;
  chomp;
  s/\r$//;      # DOS newlines
  if ($line == 1 && s/^#//) {
    @header = split /\t/, $_;
    foreach $i (0..$#header) {
      $fields->{$header[$i]} = $i;
    }
    next;
  }
  next if /^#/;
  next if /^$/;
  @f = split /\t/, $_;
  next if $f[$fields->{"SubGroup"}] ne "Gammaproteobacteria";
  next if $f[$fields->{"Status"}] ne "Complete Genome";
  if ($f[$fields->{"Replicons"}] =~ /plasmid/) {
    @g = split /\s+/, $f[$fields->{"Organism/Name"}];
    $species = join (" ", $g[0], $g[1]);
    if (!defined $species_decoy->{$species}) {
      $species_decoy->{$species} = 0;
    }
    @out = ($species, $f[$fields->{"Strain"}], $f[$fields->{"Assembly Accession"}]);
    $url = $f[$fields->{"FTP Path"}];
    @h = split /\//, $url;
    $url .= "/$h[$#h]_genomic.fna.gz";
    next if $url eq "" || $url !~ /^ftp:\/\//;
    $fna = "$tempdir/$h[$#h]_genomic.fna.gz";
    push @out, $f[$fields->{"Status"}];
    push @out, $url;
    $return = getstore($url, $fna);
    if ($fna =~ /\.gz$/) {
      open F, "zcat $fna |";
    } else {
      open F, $fna;
    }
    @seq = <F>;
    foreach $i (0..$#seq) {
      if ($seq[$i] =~ /^>/) {
        $seq[$i] =~ s/\s.*$//;
      }
    }
    $seq = fasta2hash(@seq);
    close F;
    @g = split /; /, $f[$fields->{"Replicons"}];
    foreach $i (0..$#g) {
#      next if $g[$i] !~ /^plasmid/;
      if ($g[$i] !~ /:/) {
        $g[$i+1] = $g[$i] . "; " . $g[$i+1];
        next;
      }
      $type = "";
      $name = "";
      $acc1 = "";
      $acc2 = "";
      if ($g[$i] =~ /^plasmid:(.+?)\/(.+)/) {
        $type = "Plasmid";
        $name = "";
        $acc1 = $1;
        $acc2 = $2;
      } elsif ($g[$i] =~ /^plasmid:(.+)/) {
        $type = "Plasmid";
        $name = "";
        $acc1 = $1;
        $acc2 = "";
      } elsif ($g[$i] =~ /^plasmid\s+(.+?):(.+?)\/(.+)/) {
        $type = "Plasmid";
        $name = $1;
        $acc1 = $2;
        $acc2 = $3;
      } elsif ($g[$i] =~ /^plasmid\s+(.+?):(.+)/) {
        $type = "Plasmid";
        $name = $1;
        $acc1 = $2;
        $acc2 = "";
      } elsif ($g[$i] =~ /^chromosome:(.+?)\/(.+)/) {
        $type = "Chromosome";
        $name = "";
        $acc1 = $1;
        $acc2 = $2;
      } elsif ($g[$i] =~ /^chromosome:(.+)/) {
        $type = "Chromosome";
        $name = "";
        $acc1 = $1;
        $acc2 = "";
      } elsif ($g[$i] =~ /^chromosome\s+(.+?):(.+?)\/(.+)/) {
        $type = "Chromosome";
        $name = $1;
        $acc1 = $2;
        $acc2 = $3;
      } elsif ($g[$i] =~ /^chromosome\s+(.+?):(.+)/) {
        $type = "Chromosome";
        $name = $1;
        $acc1 = $2;
        $acc2 = "";
      }
      if ($type ne "") {
        $length = 0;
        $outputted_header = "";
        $secondary_header = "";
        if ($type eq "Plasmid") {
          if (length $acc1 && defined $seq->{$acc1} && !$plasmid_included->{$acc1}) {
            if (!$printed_header->{$acc1}) {
              print O ">$acc1\n$seq->{$acc1}\n";
              $length = length($seq->{$acc1});
              $outputted_header = $acc1;
              $secondary_header = $acc2;
            }
            $printed_header->{$acc1} = 1;
          } elsif (length $acc2 && defined $seq->{$acc2} && !$plasmid_included->{$acc2}) {
            if (!$printed_header->{$acc2}) {
              print O ">$acc2\n$seq->{$acc2}\n";
              $length = length($seq->{$acc2});
              $outputted_header = $acc2;
              $secondary_header = $acc1;
            }
            $printed_header->{$acc2} = 1;
          }
          if ($acc1 !~ /^\s*$/) {
            $plasmid_included->{$acc1} = 1;
          }
          if ($acc2 !~ /^\s*$/) {
            $plasmid_included->{$acc2} = 1;
          }
        } elsif ($species_decoy->{$species} < $max_decoys) {
          if (length $acc1 && defined $seq->{$acc1} && !$plasmid_included->{$acc1}) {
            if (!$printed_header->{$acc1}) {
              print O ">$acc1\n$seq->{$acc1}\n";
              $length = length($seq->{$acc1});
              $outputted_header = $acc1;
              $secondary_header = $acc2;
            }
            $printed_header->{$acc1} = 1;
          } elsif (length $acc2 && defined $seq->{$acc2} && !$plasmid_included->{$acc1}) {
            if (!$printed_header->{$acc2}) {
              print O ">$acc2\n$seq->{$acc2}\n";
              $length = length($seq->{$acc2});
              $outputted_header = $acc2;
              $secondary_header = $acc1;
            }
            $printed_header->{$acc2} = 1;
          }
          $species_decoy->{$species}++;
          if ($acc1 !~ /^\s*$/) {
            $plasmid_included->{$acc1} = 1;
          }
          if ($acc2 !~ /^\s*$/) {
            $plasmid_included->{$acc2} = 1;
          }
        }
        if ($outputted_header ne "") {
          print M join ("\t", $outputted_header, $secondary_header, $length, $type, $name, @out, $f[$fields->{"Organism/Name"}]), "\n";
        }
      }
    }
  }
}

close P;

close O;
close M;

sub fasta2hash {
  # try to chomp newlines so we can take any format
  # try to be smart about strings vs. arrays
  my @a = @_;
  my @b = ();
  my $seq;
  my $name;
  my $a;
  my $i;
  my $test;
  foreach $a (@a) {
    push @b, split (/\n/, $a);
  }
  foreach $a (@b) {
    next if $a =~ /^#/;
    chomp $a;
    if ($a =~ s/^>//) {
      $i = 1;
      $name = $a;
      $test = $name;
      while (defined $seq->{$test}) {
        $test = $name . $i;
        $i++;
      }
      $name = $test;
      $seq->{$name} = "";
    } elsif (defined $name) {
      $seq->{$name} .= $a;
    }
  }
  return $seq;
}

