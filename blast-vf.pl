#!/usr/bin/perl -w
#
# use blastn against any fasta database
# assign best allele for groups
# use files from SRST2 - this is for full genomes
# set 90% of allele must be covered, 90% identity minimum
# (defaults for SRST2)
# then pick best blast score
#
use XML::Simple;
use File::Basename;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
my $fasta = "/home/slchen/Documents/Projects/SRST2/Ecoli-srst2.fasta";
my $min_coverage = 0.9;	# SRST2 defaults
my $min_identity = 0.9;	# SRST2 defaults
my $out_name = "";
my $report_alleles = 0;
GetOptions (
  'fasta=s' => \$fasta,
  'coverage=f' => \$min_coverage,
  'identity=f' => \$min_identity,
  'report_all_alleles' => \$report_alleles,
  'name=s' => \$out_name
);

if (!-f $fasta || !defined $ARGV[0] || !-f $ARGV[0]) {
  print "Usage: $0 <assembled genome> -fasta <SRST2 fasta file> [ -report_all_alleles ]\n";
  print "expects fasta to be in same format that SRST2 uses so genes can be grouped\n";
  exit;
}

my $xs = XML::Simple->new;
my @data = `slc-blastn -db $ARGV[0] $fasta -m 5`;
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
my $clusterid;
my $clustersym;
my $allelesym;
my $alleleid;
my $annotation;
my $cmap;
my $amap;
my $hits;
my $t;
my $hitlength;
my $hitidentity;
my $hitgaps;
my $qsequence;
my $score;
my $full_allele;
my $dbname = basename($fasta);
$dbname =~ s/\.fa$//;
$dbname =~ s/\.fna$//;
$dbname =~ s/\.fasta$//;
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
    @f = split /\s+/, $query;
    @g = split /__/, $f[0];
    if (defined $g[3]) {
      ($clusterid, $clustersym, $allelesym, $alleleid) = @g;
      $full_allele->{$clusterid}->{$alleleid} = $f[0];
      if (!defined $cmap->{$alleleid}->{$clusterid}) {
        $cmap->{$alleleid}->{$clusterid} = $clustersym;
      }
      if ($cmap->{$alleleid}->{$clusterid} ne $clustersym) {
        print STDERR "Database error, saw cluster ID $clusterid with symbols $cmap->{$alleleid}->{$clusterid} and $clustersym\n";
      }
      if (!defined $amap->{$clusterid}->{$alleleid}) {
        $amap->{$clusterid}->{$alleleid} = $allelesym;
      }
      if ($amap->{$clusterid}->{$alleleid} ne $allelesym) {
        print STDERR "Database error, saw allele ID $alleleid with symbols $amap->{$clusterid}->{$alleleid} and $allelesym\n";
      }
      if (defined $f[1]) {
        $annotation = join (" ", @f[1..$#f]);
      } else {
        $annotation = "";
      }
      if (defined $h->[$i]->{"Iteration_hits"}->{"Hit"}->{"Hit_hsps"}->{"Hsp"}) {
        if (ref($h->[$i]->{"Iteration_hits"}->{"Hit"}->{"Hit_hsps"}->{"Hsp"}) eq "HASH") {
          $t = $h->[$i]->{"Iteration_hits"}->{"Hit"}->{"Hit_hsps"}->{"Hsp"};
        } elsif (ref($h->[$i]->{"Iteration_hits"}->{"Hit"}->{"Hit_hsps"}->{"Hsp"}) eq "ARRAY") {
          # take the top hit
          $t = $h->[$i]->{"Iteration_hits"}->{"Hit"}->{"Hit_hsps"}->{"Hsp"}->[0];
        }
        $hitlength = $t->{"Hsp_align-len"};
        $hitidentity = $t->{"Hsp_identity"};
        $hitgaps = $t->{"Hsp_gaps"};
        $qsequence = $t->{"Hsp_qseq"};
        if (!$hitlength || $hitlength - $hitgaps < $min_coverage * $qlen || $hitidentity/$hitlength < $min_identity) {
          delete($hits->{$clusterid}->{$alleleid});
          next;
        }
        $hits->{$clusterid}->{$alleleid}->{QLENGTH} = $qlen;
        $hits->{$clusterid}->{$alleleid}->{QSEQ} = $qsequence;
        $hits->{$clusterid}->{$alleleid}->{IDENTITY} = $hitidentity/$hitlength;
        $hits->{$clusterid}->{$alleleid}->{COVERAGE} = ($hitlength - $hitgaps)/$qlen;
  if ($hits->{$clusterid}->{$alleleid}->{COVERAGE}> 1) {
    print "high cov query $query length $qlen hit $hitlength gaps $hitgaps id $hitidentity\n";
  }
        $hits->{$clusterid}->{$alleleid}->{ANNOTATION} = $annotation;
        $hits->{$clusterid}->{$alleleid}->{RANK} = -1;
        $hits->{$clusterid}->{$alleleid}->{DIFF} = "";
        if ($hitidentity < $hitlength - $hitgaps) {
          $hits->{$clusterid}->{$alleleid}->{DIFF} .= ($hitlength - $hitgaps - $hitidentity)."snp";
        }
        if ($hitgaps) {
          $hits->{$clusterid}->{$alleleid}->{DIFF} .= $hitgaps."indel";
        }
      }
    }
  }
}

# output should be:
#   Sample
#   DB
#   gene (seems to be clusterSymbol)
#   allele (usually alleleSymbol_alleleUniqueIdentifier)
#   coverage (in %)
#   depth (will be by definition 1000 for us)
#   diffs
#   uncertainty (blank for us)
#   divergence
#   length (seems to be of the query)
#   maxMAF (will be by definition 0 for us)
#   clusterid
#   seqid
#   annotation
foreach $i (keys %$hits) {
  foreach $j (sort
            { $hits->{$i}->{$b}->{IDENTITY} <=> $hits->{$i}->{$a}->{IDENTITY} ||
              $hits->{$i}->{$b}->{COVERAGE} <=> $hits->{$i}->{$a}->{COVERAGE} }
            keys %{$hits->{$i}}) {
    print join ("\t", $out_name,
                      $dbname,
                      $cmap->{$j}->{$i},
                      $amap->{$i}->{$j} . "_" . $j, 
                      sprintf ("%.2f", $hits->{$i}->{$j}->{COVERAGE} * 100),
                      1000,
                      $hits->{$i}->{$j}->{DIFF},
                      "",
                      sprintf("%.3f", 1 - $hits->{$i}->{$j}->{IDENTITY}),
                      $hits->{$i}->{$j}->{QLENGTH},
                      0,
                      $i,
                      $j,
                      $hits->{$i}->{$j}->{ANNOTATION}), "\n";
    last;
  }
}
my $type;
if ($report_alleles) {
  foreach $i (keys %$hits) {
    foreach $j (sort
            { $hits->{$i}->{$b}->{IDENTITY} <=> $hits->{$i}->{$a}->{IDENTITY} ||
              $hits->{$i}->{$b}->{COVERAGE} <=> $hits->{$i}->{$a}->{COVERAGE} }
            keys %{$hits->{$i}}) {
      if ($hits->{$i}->{$j}->{IDENTITY} == 1) {
        $type = "consensus";
      } else {
        $type = "variant";
      }
      print ">$full_allele->{$i}->{$j}_$j.$type $out_name.$dbname\n";
      print $hits->{$i}->{$j}->{QSEQ}, "\n";
    }
  }
}
