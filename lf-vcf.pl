#!/usr/bin/perl -w
#
# take single sample vcf from lofreq
# add in header lines and genotyping info
#
# we'll add GT, DP, FT
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
my $sample = "";
my $mindp = 10;	# lofreq default
my $minaf = 0.5;	# i.e. consensus variants only - we are assuming we're using single whole genomes
my $help = 0;
&GetOptions (
  'sample=s' => \$sample,
  'mindp=i' => \$mindp,
  'minaf=f' => \$minaf,
  'help' => \$help
);

if ($help) {
  print "Usage: $0 [ -sample <sample name> ] [ -mindp <min depth> ] [ -minaf <min AF> ] <lofreq .vcf file>\n";
  print "       $0 [ -h ]\n";
  print "Can specify filename and sample name will cut off .lofreq or .vcf\n";
  print "Or can pipe from stdin and then should specify sample name\n";
  print "Will add GT, DP, FT for genotype information to the .vcf file\n";
  print "Default minimum depth is 10, minimum AF is 0.5\n";
  print "Often will do something like:\n  $0 foo.lofreq | bgzip -c > foo.gz && tabix -p vcf foo.gz\n";
  exit;
}

if (defined $ARGV[0] && -f $ARGV[0] && !length ($sample)) {
  if ($ARGV[0] =~ /(.+)\.lofreq$/) {
    $sample = $1;
  } elsif ($ARGV[0] =~ /(.+)\.vcf$/) {
    $sample = $1;
  } else {
    $sample = $ARGV[0];
  }
}
while (<>) {
  if (/^##/) {
    print; next;
  }
  if (/^#CHROM/) {
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";
    print "##FORMAT=<ID=FT,Number=2,Type=String,Description=\"Filter\">\n";
    chomp;
    print;
    print "\tFORMAT\t$sample\n";
    next;
  }
  chomp;
  @f = split /\t/, $_;
  $filter = $f[6];
  $info = $f[7];
  $depth = 0;
  if ($info =~ /DP=(\d+)/) {
    $depth = $1;
  }
  if ($info =~ /AF=([0-9.]+)/) {
    $af = $1;
  }
  next if $depth < $mindp;
  next if $af < $minaf;
  print join ("\t", @f, "GT:DP:FT", "1/1:$depth:$filter"), "\n";
}
