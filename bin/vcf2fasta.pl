#!/usr/bin/perl -w
#
# take snp or vcf file
# take orgcode or fasta file
# make new reference sequence
# should be ok for snp_by_ref.pl script from GERMS
#
use Orgmap;
use slchen;
use strict;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
my $reference_file = "";
my $field = 0;		# in case vcf file has multiple entries - this is just
                        # a field number, doesn't actually count entries
                        # the first entry is usually field 9
                        # if field is 0 then take all the lines as they come in
                        # and not just the ones positive for a given entry
my $all = 0;		# if $all, take everything; otherwise default to
                        # filtering on PASS entries
my $snponly = 1;	# whether to use only SNPs (ignore indels)
my $min_AF = 0.5;	# minimum allele frequency
my $rename = 0;		# if rename is true, then use the header line of
			# the vcf file, only do this if one chromosome
my $force_name = "";	# override with this for output fasta name if set
my $show_help = 0;

GetOptions (
  'ref=s' => \$reference_file,
  'field=i' => \$field,
  'minaf=f' => \$min_AF,
  'all!' => \$all,
  'snponly!' => \$snponly,
  'name=s' => \$force_name,
  'rename!' => \$rename,
  'help' => \$show_help
);

my ($chr, $pos, $ref, $alt);
my (@a, @f, @alleles, $allelenum);
my $AF;
my $indel;
my @names;

if ($show_help) {
  &print_usage;
}

if (!defined $ARGV[0] && (!defined $reference_file || !-f $reference_file)) {
  &print_usage;
}

if (defined $ARGV[0] && (!defined $reference_file || !-f $reference_file)) {
  &read_orgmap($ARGV[0]);
  $reference_file = $Orgmap::fnafile;
}
if (-f $reference_file) {
  open F, $reference_file;
  @a = ();
  while (<F>) {
    if ($_ =~ /^>\S/) {
      $_ =~ s/\s.*$//;
      $_ .= "\n";
    }
    push @a, $_;
  }
  $a = slchen::fasta2hash(@a);
  close F;
} else {
  &print_usage("Can't find reference fasta file $reference_file");
}

$indel = ();
@names = ();

while (<>) {
  if (/^#CHROM/) {
    chomp;
    @names = split /\t/, $_;
  }
  next if /^#/;
  next if /^$/;
  undef $chr;
  undef $pos;
  undef $ref;
  undef $alt;
  undef $AF;
  @f = split /\t/, $_;
  if ($f[7] =~ /AF=(.*?);/) {
    $AF = $1;
    $AF = 1 if !slchen::isfloat($AF);
  }
  if (defined $f[1] && $f[1] =~ /^\d+$/) {	# probably vcf
    if ($all || ($f[6] eq 'PASS' && (!defined $AF || $AF >= $min_AF))) {
      $chr = $f[0];
      $pos = $f[1];
      $ref = $f[3];
      @alleles = split(/,/, $f[4]);
      if (!$field) {
        $allelenum = 1;
      } else {
        $allelenum = (split(/\//, $f[$field]))[0];  # just take the first one
      }
      next if $allelenum !~ /^\d/;
      next if $allelenum <= 0;
      $alt = $alleles[$allelenum-1];
    } else {
      next;
    }
  } else {
    @f = split / /, $_;
    if (defined $f[1] && $f[1] =~ /^\d+$/) {	# probably snp from lofreq
      $chr = $f[0];
      $pos = $f[1];
      $f[2] =~ /(\w+)>(\w+)/;
      $ref = $1;
      $alt = $2;
    }
  }
  if (defined $chr && defined $a->{$chr} &&
      defined $pos && $pos > 0 && $pos <= length($a->{$chr}) &&
      defined $ref && uc(substr($a->{$chr}, $pos-1, length($ref))) eq $ref &&
      defined $alt) {
    if (length($ref) != length($alt)) {
      $indel->{$pos}->{CHR} = $chr;
      $indel->{$pos}->{REF} = $ref;
      $indel->{$pos}->{ALT} = $alt;
      next;
    } else {
      substr($a->{$chr}, $pos-1, length($ref)) = $alt;
    }
  } else {
    print STDERR "Error on line:\n$_";
    print STDERR "chr $chr length ", length($a->{$chr}), "\n";
    print STDERR "pos $pos\n";
    print STDERR "ref $ref\n";
    print STDERR "alt $alt\n";
    print STDERR "line @f";
    print STDERR "but I see ", uc(substr($a->{$chr}, $pos-1, length($ref))), " at this position\n";
  }
}
if (!$snponly) {
  foreach $pos (sort { $b <=> $a } keys %$indel) {
    substr($a->{$indel->{$pos}->{CHR}}, $pos-1, length($indel->{$pos}->{REF})) = $indel->{$pos}->{ALT};
  }
}
@a = slchen::hash2fasta($a);
if ($force_name ne "" && $#a == 1 && $a[0] =~ /^>/) {
  $a[0] = ">$force_name";
} elsif ($rename && $#a == 1 && $a[0] =~ /^>/ &&
    defined $names[$field] && length($names[$field])) {
  $a[0] = ">$names[$field]";
}
print join("\n", @a), "\n";

sub print_usage {
  my ($error) = @_;
  if (defined $error) {
    print STDERR $error, "\n";
  }

print <<__HELP__;
Usage: $0 [ <orgcode> | -ref <reference fasta> ] [options] <snp or vcf file>

Take a VCF or snp file (from lofreq, for example) and a reference sequence.
Output a sequence with the variants/SNPs incorporated.

VCF or snp file comes in on standard input.
Output goes to STDOUT, be sure to redirect to a file if needed.

Options:
  -field <int>       : For VCF files with genotype data, which can accommodate
                       multiple named samples. This option is 0-based, generally
                       -field 9 is the first sample.
                       Default 0 (use all fields)
  -rename|norename   : Whether to use VCF genotype column names for fasta names
                       This only makes sense if the reference has one contig
                       Default -norename (output fasta has the same fasta names
                       as the input reference sequences)
  -name <string>     : Force output fasta to have this name
                       This only makes sense if the reference has one contig
		       Default "" (no forced name)
  -minaf <float>     : Minimum allele frequency (AF) cutoff
                       Default 0.5
  -all|noall         : Whether to use the VCF FILTER column
                       Default -noall (only use PASS lines)
  -snponly|nosnponly : Whether to use only SNPs (and ignore INDELs)
                       Default -snponly (only use SNPs)
  -help              : Show this help and exit

__HELP__

  exit;
}
