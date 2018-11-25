#!/usr/bin/perl -w
#
# use gcovs.gz or gcov.gz files to get coverage
# take sample names from vcf file
# replace no-call positions with the right info
# default PASS
# change chromosome to numbers
#
use Getopt::Long;
Getopt::Long::Configure("pass_through");
my $verbose = 0;
my $data_dir = ".";	# where to look for gcov files
my $min_cov = 10;	# based on lofreq default
my $target = "";
my $chrom = ();
my $chrom_counter = 1;
my $max_name_length = 10;
my @map = ();
my $ref_vcf = "";
my $snponly = 1;

GetOptions (
  "datadir=s" => \$data_dir,
  "verbose" => \$verbose,
  "snponly!" => \$snponly,
  "ref=s" => \$ref_vcf
);

if (!-f $ref_vcf) {
  print "Usage: $0 -ref <ref_vcf> [ -datadir <gcov dir> ] <lofreq.gz file>\n";
  print "Can take vcf file and fill in with gcov file for missing genotypes\n";
  print "Will compare input vcf samples with all samples in ref_vcf\n";
  exit;
}

if (defined $ARGV[0] && -f $ARGV[0]) {
  if ($ARGV[0] =~ /\.gz$/) {
    open IN, "zcat $ARGV[0] |";
  } else {
    open IN, $ARGV[0];
  }
  while (<IN>) {
    if (/^#CHROM/) {
      chomp;
      @map = split /\t/, $_;
      last;
    }
  }
  @in = <IN>;
  close IN;
} else {
  while (<>) {
    if (/^#CHROM/) {
      chomp;
      @map = split /\t/, $_;
      last;
    }
  }
  @in = <>;
}
my $input;
my @f;
my $alt;
my @alt;
my $skip;
foreach $i (0..$#in) {
  chomp $in[$i];
  @f = split /\t/, $in[$i];
  @alt = split /,/, $f[4];
  $skip = 0;
  unshift @alt, $f[3];
  if ($snponly) {
    foreach $alt (@alt) {
      if (length($alt) > 1) {
        $skip = 1;
      }
    }
  }
  next if $skip;
  if (!defined $chrom->{$f[0]}) {
    $chrom->{$f[0]} = $chrom_counter;
    $chrom_counter++;
  }
  foreach $j (9..$#map) {
    $code = substr($f[$j], 0, 1);
    if ($code =~ /^\d+$/) {
      # need some logic here
      # sometimes there are two lines for the same position, one for snps and one for an indel
      # usually if it's one then the other will be "."
      # in this case if one has a variant called let the variant take precedence and ignore the "."
      $input->[$chrom->{$f[0]}]->{$f[1]}->{$j}->{$f[3]} = $alt[$code];
    } else {
      if (!defined $input->[$chrom->{$f[0]}]->{$f[1]}->{$j}->{$f[3]}) {
        $input->[$chrom->{$f[0]}]->{$f[1]}->{$j}->{$f[3]} = ".";
      }
    }
  }
}

my @gfiles;
foreach $i (9..$#map) {
  if (-f "$data_dir/$map[$i].gcovs.gz") {
    open ($gfiles[$i], "zcat $data_dir/$map[$i].gcovs.gz |");
  } elsif (-f "$data_dir/$map[$i].gcov.gz") {
    open ($gfiles[$i], "zcat $data_dir/$map[$i].gcov.gz |");
  } elsif (-f "$data_dir/$map[$i].gcovs") {
    open ($gfiles[$i], "$data_dir/$map[$i].gcovs");
  } elsif (-f "$data_dir/$map[$i].gcov") {
    open ($gfiles[$i], "$data_dir/$map[$i].gcov");
  } else {
    print STDERR "Can't find $map[$i].gcovs file in $data_dir...\n";
  }
}

if ($ref_vcf =~ /\.gz$/) {
  open R, "zcat $ref_vcf |";
} else {
  open R, $ref_vcf;
}
my @compare = ();
while (<R>) {
  next if /^$/;
  next if /^##/;
  if (/^#CHROM/) {
    chomp;
    @compare = split /\t/, $_;
    last;
  }
}
my $result;
$chrom_counter = 1;
my $last_cov;
my $cov;
my $fh;
if (scalar @compare) {
  while (<R>) {
    chomp;
    @f = split /\t/, $_;
    @alt = split /,/, $f[4];
    unshift @alt, $f[3];
    $skip = 0;
    if ($snponly) {
      foreach $alt (@alt) {
        if (length($alt) > 1) {
          $skip = 1;
        }
      }
    }
    next if $skip;
    if (defined $chrom->{$f[0]}) {
      $chromosome = $chrom->{$f[0]};
    } elsif ($f[0] =~ /^\d+$/) {
      $chromosome = $f[0];
    } else {
      $chrom->{$f[0]} = $chrom_counter;
      $chromosome = $chrom->{$f[0]};
      $chrom_counter++;
    }
    foreach $i (9..$#map) {
      foreach $j (9..$#f) {
        # $map[$i] is the sample name from input vcf
        # $f[$j] is the sample name from the reference vcf
        # @f are fields from the reference vcf, $f[1] is position
        # @g are fields from the gcov file
        $in = "";
        if (!defined $input->[$chromosome]->{$f[1]}->{$i}->{$f[3]}) {
          if (defined $last_cov->{$map[$i]}->{$f[1]}) {
            if ($last_cov->{$map[$i]}->{$f[1]} >= $min_cov) {
              $in = $f[3];
            } else {
              $in = ".";
            }
          } else {
            $fh = $gfiles[$i];
            while ($cov = <$fh>) {
              next if $cov =~ /^$/;
              next if $cov =~ /^#/;
              chomp $cov;
              @g = split /\t/, $cov;
              if ((defined $chrom->{$g[0]} && $chrom->{$g[0]} eq $chromosome || $g[0] eq $chromosome) && $g[1] == $f[1]) {
                $last_cov->{$map[$i]}->{$f[1]} = $g[2];
                if ($g[2] >= $min_cov) {
                  $in = $f[3];
                } else {
                  $in = ".";
                }
                last;
              }
            }
            if ($in eq "") {
              print "Couldn't find position $f[1] for gcov for $map[$i]\n";
              $in = ".";
            }
          }
        } else {
          $in = $input->[$chromosome]->{$f[1]}->{$i}->{$f[3]};
        }
        $code = substr($f[$j], 0, 1);
        if ($code =~ /^\d+$/) {
          if (defined $alt[$code]) {
            $comp = $alt[$code];
          } else {
            print STDERR "Error assigning allele on line $_\n";
          }
        } else {
          $comp = ".";
        }
        if ($in eq $comp) {
          if ($in eq ".") {
            $result->{$map[$i]}->{$compare[$j]}->{NOCOV}++;
          } else {
            $result->{$map[$i]}->{$compare[$j]}->{IDENT}++;
          }
        } elsif ($in eq "." || $comp eq ".") {
          $result->{$map[$i]}->{$compare[$j]}->{INDEL}++;
#print "$in $comp $_\n";
        } else {
          if (length $f[3] == length $f[4]) {
            $result->{$map[$i]}->{$compare[$j]}->{SNP}++;
          } else {
            $result->{$map[$i]}->{$compare[$j]}->{INDEL}++;
          }
        }
      }
    }
  }
  print "Compare";
  foreach $i (9..$#map) {
    print "\t", join ("\t", "$map[$i]-ID", "$map[$i]-SNP", "$map[$i]-INDEL", "$map[$i]-NOCOV");
    print "\t$map[$i]";
  }
  print "\n";
  foreach $j (9..$#compare) {
    print "$compare[$j]";
    foreach $i (9..$#map) {
      foreach $type (qw(IDENT SNP INDEL NOCOV)) {
        if (!defined $result->{$map[$i]}->{$compare[$j]}->{$type}) {
          $result->{$map[$i]}->{$compare[$j]}->{$type} = 0;
        }
      }
      my $snp = $result->{$map[$i]}->{$compare[$j]}->{SNP};
      my $id = $result->{$map[$i]}->{$compare[$j]}->{IDENT};
      my $indel = $result->{$map[$i]}->{$compare[$j]}->{INDEL};
      print "\t", join ("\t", $result->{$map[$i]}->{$compare[$j]}->{IDENT}, $result->{$map[$i]}->{$compare[$j]}->{SNP}, $result->{$map[$i]}->{$compare[$j]}->{INDEL}, $result->{$map[$i]}->{$compare[$j]}->{NOCOV});
      print "\t", ($snp + $indel)/($snp + $id + $indel);
    }
    print "\n";
  }
}
