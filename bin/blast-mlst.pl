#!/usr/bin/perl -w
#
# use blastn against MLST database
# pull in definitions from a text file
# assign MLST
# use files from SRST2 - this is for full genomes
# set 80% of allele must be covered, then look for highest % identity
#
use XML::Simple;
use FindBin;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

# We expect most files to be in somewhere like /usr/local/lib/SRST2/MLST
# This is where we will look in the procedure set_data at the bottom

my $mlst_fasta = "";
my $mlst_def = "";
my $delimiter = "-";
my $species = "";
my $meta = &set_data;

my $min_coverage = 0.9;	# SRST2 default
my $min_identity = 0.9;	# SRST2 default
my $debug = 0;
my $out_name = "";
GetOptions (
  'mlst_fasta=s' => \$mlst_fasta,
  'mlst_def=s' => \$mlst_def,
  'delimiter=s' => \$delimiter,
  'species=s' => \$species,
  'name=s' => \$out_name,
  'debug!' => \$debug
);

if (length $species) {
  if (defined $meta->{$species}) {
    $mlst_fasta = $meta->{$species}->{FASTA};
    $mlst_def = $meta->{$species}->{DEFINITIONS};
    $delimiter = $meta->{$species}->{DELIMITER};
  } else {
    print "Unknown -species argument $species\n";
    undef($ARGV[0]);
  }
}

if (!-f $mlst_fasta || !-f $mlst_def || !defined $ARGV[0] || !-f $ARGV[0]) {
  print "Usage: $0 <assembled genome> -mlst_fasta <mlst fasta file> -mlst_def <mlst definition file> -delimiter <char>\n";
  print "       $0 <assembled genome> -species <species>\n";
  print "expects mlst_fasta and mlst_def to be in same format that SRST2 uses\n";
  print "need to know the allele delimiter, i.e \"-\" for E. coli or \"_\" for GBS\n";
  print "Alternatively, can specify species from these:\n";
  foreach (sort keys %$meta) {
    print "  $_\n";
  }
  exit;
}

my $xs = XML::Simple->new;
my @data = `slc-blastn -db $ARGV[0] $mlst_fasta -m 5`;
foreach $i (0..$#data) {
  chomp $data[$i];
}
my $data = join ("", @data);
my $ref = $xs->XMLin($data);
#use Data::Dumper;
#print Dumper $ref;
#exit;

my $h;
my $query;
my $qlen;
my $gene;
my $allele;
my $hits;
my $length;
my $t;
my $hitlength;
my $hitidentity;
my $hitgaps;
my $score;
if ($out_name eq "") {
  $out_name = $ARGV[0];
  $out_name =~ s/\.fna$//;
  $out_name =~ s/\.fa$//;
  $out_name =~ s/\.fasta$//;
}

if (defined $ref->{"BlastOutput_iterations"}->{"Iteration"}) {
  $h = $ref->{"BlastOutput_iterations"}->{"Iteration"};
  foreach $i (0..$#$h) {
    $query = $h->[$i]->{"Iteration_query-def"};
    $qlen = $h->[$i]->{"Iteration_query-len"};
    if ($query =~ /(.+)$delimiter(\d+)$/) {
      ($gene, $allele) = ($1, $2);
      $length = $h->[$i]->{"Iteration_query-len"};
      $hitlength = 0;
      $hitidentity = 0;
      $hitgaps = 0;
      if (defined $h->[$i]->{"Iteration_hits"}->{"Hit"}->{"Hit_hsps"}->{"Hsp"}) {
        if (ref($h->[$i]->{"Iteration_hits"}->{"Hit"}->{"Hit_hsps"}->{"Hsp"}) eq "HASH") {
          $t = $h->[$i]->{"Iteration_hits"}->{"Hit"}->{"Hit_hsps"}->{"Hsp"};
        } elsif (ref($h->[$i]->{"Iteration_hits"}->{"Hit"}->{"Hit_hsps"}->{"Hsp"}) eq "ARRAY") {
          $t = $h->[$i]->{"Iteration_hits"}->{"Hit"}->{"Hit_hsps"}->{"Hsp"}->[0];
        }
        $hitlength = $t->{"Hsp_align-len"};
        $hitidentity = $t->{"Hsp_identity"};
        $hitgaps = $t->{"Hsp_gaps"};

        if (!$hitlength || $hitlength - $hitgaps < $min_coverage * $qlen || $hitidentity/$hitlength < $min_identity) {
          delete($hits->{$gene}->{$allele});
          next;
        }
        $hits->{$gene}->{$allele}->{IDENTITY} = $hitidentity/$hitlength;
        $hits->{$gene}->{$allele}->{COVERAGE} = ($hitlength - $hitgaps)/$qlen;
        $hits->{$gene}->{$allele}->{RANK} = -1;
        $hits->{$gene}->{$allele}->{DIFF} = "";
        if ($hitidentity < $hitlength - $hitgaps) {
          $hits->{$gene}->{$allele}->{DIFF} .= ($hitlength - $hitidentity - $hitgaps) . "snp";
        }
        if ($hitgaps) {
          $hits->{$gene}->{$allele}->{DIFF} .= $hitgaps . "indel";
        }
      }
    }
  }
}

# rank everything
my $rank;
my $best;
foreach $gene (keys %$hits) {
  $rank = 0;
  foreach $allele (sort
      { $hits->{$gene}->{$b}->{IDENTITY} <=> $hits->{$gene}->{$a}->{IDENTITY} ||
        $hits->{$gene}->{$b}->{COVERAGE} <=> $hits->{$gene}->{$a}->{COVERAGE} }
      keys %{$hits->{$gene}}) {
    $hits->{$gene}->{$allele}->{RANK} = $rank;
    if ($rank == 0) {
      $best->{$gene} = $allele;
    }
    $rank++;
  }
}

#my $min_sum1 = -1;	# one mismatched allele
#my $min_sum2 = -1;	# one missing allele
#my $min_sum3 = -1;	# two mismatched alleles
#my $min_sum4 = -1;	# two missing alleles
my $sum = 0;
my $cc = "";
my $close = ();
my @nonzero = ();
my @noallele = ();
my $note = "";
my $header = "# Sample\tST";
my @out1;
my @out2;
my @allele_out = ();
my @diff = ();
my $st = 0;
my $missing = 0;
my %mlst_genes = ();
open M, $mlst_fasta;
while (<M>) {
  chomp;
  if (s/^>//) {
    @f = split /$delimiter/, $_;
    if ($f[$#f] =~ /^\d+$/) {
      $mlst_genes->{join ($delimiter, @f[0..$#f-1])} = 1;
    }
  }
}
close M;
open M, $mlst_def;
while (<M>) {
  chomp;
  @f = split /\t/, $_;
  if (/^ST/) {
    foreach $i (1..$#f) {
#      next if $f[$i] eq "clonal_complex";
      next if !defined $mlst_genes->{$f[$i]};
      $map->[$i] = $f[$i];
      $header .= "\t$f[$i]";
      if (defined $hits->{$f[$i]} && defined $best->{$f[$i]}) {
        if ($hits->{$f[$i]}->{$best->{$f[$i]}}->{DIFF} eq "") {
          push @allele_out, $best->{$f[$i]};
        } else {
          push @allele_out, $best->{$f[$i]} . "*";
          push @diff, $hits->{$f[$i]}->{$best->{$f[$i]}}->{DIFF};
        }
      } else {
        push @allele_out, "-";
        $missing++;
      }
    }
    $header .= "\tmismatches\tuncertainty\tdepth\tmaxMAF";
    if (scalar @diff) {
      push @allele_out, join (";", @diff);
    } else {
      push @allele_out, "";
    }
    push @allele_out, "", 1000, 0;
    next;
  }

  @noallele = ();
  @nonzero = ();
  $sum = 0;
  foreach $i (1..$#$map) {
    if (defined $hits->{$map->[$i]}->{$f[$i]}->{RANK} && $hits->{$map->[$i]}->{$f[$i]}->{RANK} >= 0) {
      $debug && print "$map->[$i] $f[$i] Rank:$hits->{$map->[$i]}->{$f[$i]}->{RANK}\n";
      if ($hits->{$map->[$i]}->{$f[$i]}->{RANK} > 0) {
        push @nonzero, $i;
      }
      $sum += $hits->{$map->[$i]}->{$f[$i]}->{RANK};
    } else {
      push @noallele, $i;
    }
  }
  $debug && print "ST$f[0] Sum:$sum; nonzero:", join (",", @nonzero), "; noallele:", join (",", @noallele), "\n";

  if (scalar @nonzero == 0 && scalar @noallele == 0 && $sum == 0) {
    # a perfect allele
    $st = $f[0];
    last;
  }
#  } elsif (scalar @nonzero == 1 && scalar @noallele == 0 && $sum > 0) {
#    if ($sum < $min_sum1 || $min_sum1 == -1) {
#      $min_sum1 = $sum;
#      my $diff = $map->[$nonzero[0]] . "-" . $f[$nonzero[0]] . "/";
#      $diff .= $hits->{$map->[$nonzero[0]]}->{$f[$nonzero[0]]}->{DIFF};
#      $f[0] .= "*";
#      $f[$nonzero[0]] .= "*";
#      @out1 = ($out_name, @f[0..$#$map], $diff, "", 1000, 0);
#      if ($#f > $#$map) {
#        push @out1, @f[$#$map+1..$#f];
#      }
#    }
#  } elsif (scalar @noallele == 1) {
#    if ($sum < $min_sum2 || $min_sum2 == -1) {
#      $min_sum2 = $sum;
#      $note = "No allele for $map->[$noallele[0]]";
#      $f[0] .= "*";
#      $f[$noallele[0]] = "-";
#      @out2 = ($out_name, @f, $note);
#    }
#  } elsif (scalar @nonzero == 2 && scalar @noallele == 0 && $sum > 0) {
#    if ($sum < $min_sum3 || $min_sum3 == -1) {
#      $min_sum3 = $sum;
#      $note = "Has $map->[$nonzero[0]]$delimiter$best->{$map->[$nonzero[0]]}, $map->[$nonzero[1]]$delimiter$best->{$map->[$nonzero[1]]}";
#      if (defined $hits->{$map->[$nonzero[0]]}->{$f[$nonzero[0]]}->{NOTES}) {
#        $note .= "; $hits->{$map->[$nonzero[0]]}->{$f[$nonzero[0]]}->{NOTES}";
#      }
#      $f[0] .= "*";
#      $f[$nonzero[0]] .= "*";
#      $f[$nonzero[1]] .= "*";
#      @out3 = ($out_name, @f, $note);
#    }
#  } elsif (scalar @noallele == 2) {
#    if ($sum < $min_sum4 || $min_sum4 == -1) {
#      $min_sum4 = $sum;
#      $note = "No allele for $map->[$noallele[0]] and $map->[$noallele[1]]";
#      $f[0] .= "*";
#      $f[$noallele[0]] = "-";
#      $f[$noallele[1]] = "-";
#      @out4 = ($out_name, @f, $note);
#    }
#  }
}
close M;

# SRST2 outputs this:
#   Sample
#   ST
#   adk
#   fumC
#   gyrB
#   icd
#   mdh
#   purA
#   recA
#   mismatches
#   uncertainty (should be nothing for us)
#   depth (defined as 1000 for us)
#   maxMAF (defined as 0 for us)
print "$header\n";
if ($st) {
  # check for differences, these are summarized in $allele_out[$#allele_out-3]
  if ($allele_out[$#allele_out-3] eq "") {
    print join ("\t", $out_name, $st, @allele_out), "\n";
  } else {
    print join ("\t", $out_name, $st . "*", @allele_out), "\n";
  }
} elsif ($missing) {
  print join ("\t", $out_name, "NF*", @allele_out), "\n";
} else {
  print join ("\t", $out_name, "NF", @allele_out), "\n";
}

sub set_data {
  my $meta;
  my $mlst_base = $FindBin::Bin . "/../lib/SRST2/MLST";

  $meta->{"Abaumannii"}->{FASTA} = "$mlst_base/Acinetobacter_baumannii#1.fasta";
  $meta->{"Abaumannii"}->{DEFINITIONS} = "$mlst_base/abaumannii.txt";
  $meta->{"Abaumannii"}->{DELIMITER} = "_";

  $meta->{"Ecoli"}->{FASTA} = "$mlst_base/Escherichia_coli#1.fasta";
  $meta->{"Ecoli"}->{DEFINITIONS} = "$mlst_base/ecoli.txt";
  $meta->{"Ecoli"}->{DELIMITER} = "_";

  $meta->{"GBS"}->{FASTA} = "$mlst_base/Streptococcus_agalactiae.fasta";
  $meta->{"GBS"}->{DEFINITIONS} = "$mlst_base/sagalactiae.txt";
  $meta->{"GBS"}->{DELIMITER} = "_";

  $meta->{"Sagalactiae"}->{FASTA} = "$mlst_base/Streptococcus_agalactiae.fasta";
  $meta->{"Sagalactiae"}->{DEFINITIONS} = "$mlst_base/sagalactiae.txt";
  $meta->{"Sagalactiae"}->{DELIMITER} = "_";

  $meta->{"Kpneumoniae"}->{FASTA} = "$mlst_base/Klebsiella_pneumoniae.fasta";
  $meta->{"Kpneumoniae"}->{DEFINITIONS} = "$mlst_base/kpneumoniae.txt";
  $meta->{"Kpneumoniae"}->{DELIMITER} = "_";

  $meta->{"Spneumoniae"}->{FASTA} = "$mlst_base/Klebsiella_pneumoniae.fasta";
  $meta->{"Spneumoniae"}->{DEFINITIONS} = "$mlst_base/kpneumoniae.txt";
  $meta->{"Spneumoniae"}->{DELIMITER} = "_";

  $meta->{"Cjejuni"}->{FASTA} = "$mlst_base/Campylobacter_jejuni.fasta";
  $meta->{"Cjejuni"}->{DEFINITIONS} = "$mlst_base/campylobacter.txt";
  $meta->{"Cjejuni"}->{DELIMITER} = "_";

  $meta->{"Lmonocytogenes"}->{FASTA} = "$mlst_base/Listeria_monocytogenes.fasta";
  $meta->{"Lmonocytogenes"}->{DEFINITIONS} = "$mlst_base/lmonocytogenes.txt";
  $meta->{"Lmonocytogenes"}->{DELIMITER} = "_";

  $meta->{"GAS"}->{FASTA} = "$mlst_base/Streptococcus_pyogenes.fasta";
  $meta->{"GAS"}->{DEFINITIONS} = "$mlst_base/spyogenes.txt";
  $meta->{"GAS"}->{DELIMITER} = "_";

  $meta->{"Spyogenes"}->{FASTA} = "$mlst_base/Streptococcus_pyogenes.fasta";
  $meta->{"Spyogenes"}->{DEFINITIONS} = "$mlst_base/spyogenes.txt";
  $meta->{"Spyogenes"}->{DELIMITER} = "_";

  $meta->{"Senterica"}->{FASTA} = "$mlst_base/Salmonella_enterica.fasta";
  $meta->{"Senterica"}->{DEFINITIONS} = "$mlst_base/senterica.txt";
  $meta->{"Senterica"}->{DELIMITER} = "_";

  $meta->{"Efaecalis"}->{FASTA} = "$mlst_base/Enterococcus_faecalis.fasta";
  $meta->{"Efaecalis"}->{DEFINITIONS} = "$mlst_base/efaecalis.txt";
  $meta->{"Efaecalis"}->{DELIMITER} = "_";

  $meta->{"Efaecium"}->{FASTA} = "$mlst_base/Enterococcus_faecium.fasta";
  $meta->{"Efaecium"}->{DEFINITIONS} = "$mlst_base/efaecium.txt";
  $meta->{"Efaecium"}->{DELIMITER} = "_";

  $meta->{"Ecloacae"}->{FASTA} = "$mlst_base/Enterobacter_cloacae.fasta";
  $meta->{"Ecloacae"}->{DEFINITIONS} = "$mlst_base/ecloacae.txt";
  $meta->{"Ecloacae"}->{DELIMITER} = "_";

  $meta->{"Paeruginosa"}->{FASTA} = "$mlst_base/Pseudomonas_aeruginosa.fasta";
  $meta->{"Paeruginosa"}->{DEFINITIONS} = "$mlst_base/paeruginosa.txt";
  $meta->{"Paeruginosa"}->{DELIMITER} = "_";

  $meta->{"Saureus"}->{FASTA} = "$mlst_base/Staphylococcus_aureus.fasta";
  $meta->{"Saureus"}->{DEFINITIONS} = "$mlst_base/saureus.txt";
  $meta->{"Saureus"}->{DELIMITER} = "_";

  $meta->{"IncAC"}->{FASTA} = "$mlst_base/Plasmid/IncAC.fasta";
  $meta->{"IncAC"}->{DEFINITIONS} = "$mlst_base/Plasmid/IncAC.txt";
  $meta->{"IncAC"}->{DELIMITER} = "_";

  $meta->{"IncF"}->{FASTA} = "$mlst_base/Plasmid/IncF.fasta";
  $meta->{"IncF"}->{DEFINITIONS} = "$mlst_base/Plasmid/IncF.txt";
  $meta->{"IncF"}->{DELIMITER} = "_";

  $meta->{"IncHI1"}->{FASTA} = "$mlst_base/Plasmid/IncHI1.fasta";
  $meta->{"IncHI1"}->{DEFINITIONS} = "$mlst_base/Plasmid/IncHI1.txt";
  $meta->{"IncHI1"}->{DELIMITER} = "_";

  $meta->{"IncHI2"}->{FASTA} = "$mlst_base/Plasmid/IncHI2.fasta";
  $meta->{"IncHI2"}->{DEFINITIONS} = "$mlst_base/Plasmid/IncHI2.txt";
  $meta->{"IncHI2"}->{DELIMITER} = "_";

  $meta->{"IncI1"}->{FASTA} = "$mlst_base/Plasmid/IncI1.fasta";
  $meta->{"IncI1"}->{DEFINITIONS} = "$mlst_base/Plasmid/IncI1.txt";
  $meta->{"IncI1"}->{DELIMITER} = "_";

  $meta->{"IncN"}->{FASTA} = "$mlst_base/Plasmid/IncN.fasta";
  $meta->{"IncN"}->{DEFINITIONS} = "$mlst_base/Plasmid/IncN.txt";
  $meta->{"IncN"}->{DELIMITER} = "_";

  $meta->{"Bpertussis"}->{FASTA} = "$mlst_base/Bordetella.fasta";
  $meta->{"Bpertussis"}->{DEFINITIONS} = "$mlst_base/bordetella.txt";
  $meta->{"Bpertussis"}->{DELIMITER} = "_";

  foreach $i (keys %$meta) {
    if (!-f $meta->{$i}->{FASTA} || !-f $meta->{$i}->{DEFINITIONS}) {
      delete($meta->{$i});
    }
  }

  return $meta;
}
