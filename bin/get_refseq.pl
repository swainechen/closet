#!/usr/bin/perl -w
#
# take some unique string
# get the fna file from genbank's new ftp structure
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
use Net::FTP;
use File::Basename;
my $ext = "fna";
my $verbose = 0;
my $mapping = {
  "fna" => "genomic.fna.gz",
  "gbff" => "genomic.gbff.gz",
  "gff" => "genomic.gff.gz",
  "faa" => "protein.faa.gz",
  "gpff" => "protein.gpff.gz",
  "cds_from_genomic.fna" => "cds_from_genomic.fna.gz",
  "feature_table.txt" => "feature_table.txt.gz",
  "rna_from_genomic.fna" => "rna_from_genomic.fna.gz"
};

if (!defined $ARGV[0] || !length($ARGV[0])) {
  print "Usage: $0 <accession> [ -extension <ext> ]\n";
  print "  Typical accession would be GCA_000195995.1 or NC_003198.1\n";
  print "  Reasonable extensions are:\n";
  print "    fna (default)\n";
  print "    gbff\n";
  print "    gff\n";
  print "    faa\n";
  print "    gpff\n";
  print "    cds_from_genomic.fna\n";
  print "    feature_table.txt\n";
  print "    rna_from_genomic.fna\n";
  exit;
}

GetOptions (
  'extension=s' => \$ext,
  'verbose!' => \$verbose
);

if (!defined $ENV{GERMS_DATA} || !-d $ENV{GERMS_DATA}) {
  print STDERR "Need to define \$GERMS_DATA environment variable and make sure it exists\n";
  print STDERR "genomes_proks.txt needs to be there as well - download from http://www.ncbi.nlm.nih.gov/genome/browse/#\n";
  exit;
}

$lproks = "$ENV{GERMS_DATA}/genomes_proks.txt";

if (!-f $lproks) {
  print STDERR "Need to download genomes_proks.txt from http://www.ncbi.nlm.nih.gov/genome/browse/#\n";
  exit;
}

my $search = $ARGV[0];
my @hit = ();
my @f;
my $url;
my $suffix;
my $refseq;
my $i;
open F, $lproks;
while (<F>) {
  next if /^#/;
  next if /^$/;
  chomp;
  if (/$search/) {
    push @hit, $_;
  }
}
close F;
if (scalar @hit) {
  if ($verbose) {
    print STDERR "Found ", scalar @hit, " matching lines in $lproks\n";
  }
  foreach $i (@hit) {
    @f = split /\t/, $i;
    $base_url = $f[18];
    $suffix = $ext;
    if ($base_url =~ /^ftp:\/\/(.+?)(\/.+)$/) {
      $server = $1;
      $path = $2;
      @g = split /\//, $path;
      $refseq = $g[$#g];
      if (defined $mapping->{$ext}) {
        $suffix = $mapping->{$ext};
      }
      $url = "$path/$refseq" . "_" . $suffix;
      $out = File::Basename::basename($url);
      if ($verbose) {
        print STDERR "Retrieving ftp://$server/$url\n";
      }
      $ftp = Net::FTP->new($server) || die "Can't connect to $server: $@";
      $ftp->login();
      $ftp->binary();
      $ftp->get($url) || die "Can't get file $url", $ftp->message;
      if ($out =~ s/\.gz$//) {
        system("gunzip $out.gz");
      }
      print "$out\n";
    } else {
      next;
    }
    exit;
  }
} else {
  print STDERR "Can't find $search in genomes_proks.txt file\n";
  exit;
}
