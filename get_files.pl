#!/usr/bin/perl -w
#
# try to figure things out
# GIS is [DWR][A-Z]B###(#) usually
# URL should go to GenBank SRA
#
# use GERMS_DATA environment variable
#
use File::Basename;
use Cwd;
use Getopt::Long;
use DBI;
&Getopt::Long::Configure("pass_through");

# some globals
my $SRADB = "/mnt/genomeDB/ncbi/sra/SRAmetadb.sqlite";
my $DBH;
my %ATTR;

my $dir = "";
my $user = $ENV{USER};
my $fqd = "/mnt/software/stow/sra-toolkit-2.3.5-2/bin/fastq-dump";
my $combine = 1;	# whether to combine multiple files
my $delimiter = "\n";	# for output
my $delete = 1;		# whether to delete files after combining
my $debug = 0;
my @out;
my $out;
my $tries;
my $study;
my $sample;
my $url;
my $current;
my @f;

GetOptions (
  'debug!' => \$debug,
  'combine!' => \$combine,
  'delete!' => \$delete,
  'delimiter=s' => \$delimiter
);

if (defined $ENV{GERMS_DATA}) {
  $dir = $ENV{GERMS_DATA};
} else {
  die "Please set GERMS_DATA environment variable\n";
}
if (!defined $ARGV[0] || !length($ARGV[0])) {
  print "Usage: $0 <ID> [ -delimiter <char> ] [ -debug ] [ -nocombine ] [ -nodelete ]\n";
  print "ID is either some GIS sequencing library ID or a Genbank SRA URL.\n";
  print "This script relies on the GERMS_DATA environment variable being set\n";
  print "Default is to combine - for example HiSeq4K comes in 4 fastq files otherwise\n";
  print "Default is to delete files that are combined, leaving only one set of ID-combined_R[12].fastq.gz files\n";
  print "Default delimiter is newline, this is the separator for multiple filenames in output\n";
  exit;
}

$current = Cwd::getcwd();
chdir $dir;
@out = q12($ARGV[0], $dir);
if (scalar @out) {
  if ($combine && scalar @out > 2) {
    # don't delete since we didn't generate source files
    @out = combine($ARGV[0], $dir, $delete);
  }
  print join ($delimiter, @out), "\n";
  exit;
}

if ($ARGV[0] =~ /[SED]R[APSXR]\d+/) {
  $url = "";
  $db_query = sra_query($ARGV[0]);
  if (!defined $db_query || scalar keys %$db_query == 0) {
    die "Can't find data for $ARGV[0].\n";
  }
  if (scalar keys %$db_query > 1) {
    print STDERR "Ambiguous identifier $ARGV[0] - potential matches:\n";
    foreach $i (keys %$db_query) {
      print "$i\t$db_query->{$i}\t", make_url($i, $db_query->{$i});
    }
    die;
  }
  foreach $i (keys %$db_query) {
    $run = $i;
    $url = make_url($i, $db_query->{$i});
    last;
  }
  if (!length $url) {
    die "Couldn't figure out what $ARGV[0] is.\n";
  }
  print STDERR "$url\n" if $debug;
  @out = sra_files($run, $dir);
  if (scalar @out) {
    print join ($delimiter, @out), "\n";
    exit;
  }
  if ($debug) {
    print STDERR "Trying to get files from GenBank SRA. Running cluster command:\nssh -l $user aquila qsub -sync y -pe OpenMP 1 -l mem_free=1G,h_rt=1:59:00 -cwd -V -b y -q ionode.q -o /dev/null -e /dev/null -N $ARGV[0]-sra '\"wget $url --tries=5 --quiet -O $dir/$run.sra && $fqd --split-files --gzip -O $dir $dir/$run.sra && rm -f $dir/$run.sra\"' > /dev/null 2>&1\n";
  }
  system ("ssh -l $user aquila qsub -sync y -pe OpenMP 1 -l mem_free=1G,h_rt=1:59:00 -cwd -V -b y -q ionode.q -o /dev/null -e /dev/null -N $ARGV[0]-sra '\"wget $url --tries=5 --quiet -O $dir/$run.sra && $fqd --split-files --gzip -O $dir $dir/$run.sra && rm -f $dir/$run.sra\"' > /dev/null 2>&1");
  @out = sra_files($run, $dir);
  if (scalar @out) {
    print join ($delimiter, @out), "\n";
  }
} elsif ($ARGV[0] =~ /[DWR][A-Z]B\d\d\d+/) {
  if ($debug) {
    print STDERR "Trying to get files using SRAQuery. Running cluster command:\nssh -l $user aquila qsub -sync y -pe OpenMP 2 -l mem_free=4G,h_rt=1:59:00 -cwd -V -b y -q ionode.q -o /dev/null -e /dev/null -N $ARGV[0]-fastq \"java -Xmx3G -jar /mnt/software/bin/SRAQuery.jar -q \\\\\\\":library-id = \\\'$ARGV[0]\\\'\\\\\\\" -extract $dir\" > /dev/null 2>&1\n";
  }
  system ("ssh -l $user aquila qsub -sync y -pe OpenMP 2 -l mem_free=4G,h_rt=1:59:00 -cwd -V -b y -q ionode.q -o /dev/null -e /dev/null -N $ARGV[0]-fastq \"java -Xmx3G -jar /mnt/software/bin/SRAQuery.jar -q \\\\\\\":library-id = \\\'$ARGV[0]\\\'\\\\\\\" -extract $dir\" > /dev/null 2>&1");
  if ($? != 0) {
    die "Some error with SRAQuery job...\n";
  }
  @out = q12($ARGV[0], $dir);
  if ($combine && scalar @out > 2) {
    # delete source files because we should have generated them
    @out = combine($ARGV[0], $dir, $delete);
  }
  print join ($delimiter, @out), "\n";
} else {
  chdir $current;
  die "Couldn't figure out what $ARGV[0] is.\n";
}
chdir $current;

sub sra_files {
  my ($run, $dir) = @_;
  my $f;
  my @r = ();
  opendir(D, $dir);
  while ($f = readdir D) {
    push @r, "$dir/$f" if $f =~ /($run)_\d\.fastq\.gz/;
  }
  @r = sort { $a cmp $b } (@r);
  return @r;
}

sub sra_query {
  my ($search) = @_;
  my $q;
  my $sql;
  if (!defined $DBH) {
    $DBH = DBI->connect("DBI:SQLite:dbname=$SRADB", "", "", \%ATTR);
  }
  if ($search =~ /^.RR\d/) {
    $sql = "SELECT run_accession, experiment_accession FROM sra WHERE run_accession = ?";
  } elsif ($search =~ /^.RX\d/) {
    $sql = "SELECT run_accession, experiment_accession FROM sra WHERE experiment_accession = ?";
  } elsif ($search =~ /^.RS\d/) {
    $sql = "SELECT run_accession, experiment_accession FROM sra WHERE sample_accession = ?";
  }
  my $sth = $DBH->prepare($sql);
  $sth->execute($search);
  my $r;
  my @r;
  while (@r = $sth->fetchrow_array) {
    $r->{$r[0]} = $r[1];
  }
  return $r;
}

sub q12 {
  my ($lib, $dir) = @_;
  my $i = 0;
  my $j = 0;
  my @return = ();
  my @filelist = ();
  my $file;
  my @combined = ();
  return @return if !defined $lib || !length $lib || !-d $dir;
  if (-d "$dir/$lib.h5") {
    opendir D, "$dir/$lib.h5";
    while ($file = readdir(D)) {
      push @filelist, $file;
    }
    @filelist = sort { $a cmp $b } @filelist;
    foreach $file (@filelist) {
      # standard HiSeq
      if ($file =~ /_R[12]\.fastq$/     || $file =~ /_R[12]\.fastq.gz$/ ||
          $file =~ /_R[12]_\d+\.fastq$/ || $file =~ /_R[12]_\d+\.fastq.gz$/ ||
          $file =~ /_[12]\.fq$/      || $file =~ /_[12]\.fq.gz$/
         ) {
        push @return, "$dir/$lib.h5/$file";
        push @combined, "$dir/$lib.h5/$file" if $file =~ /-combined_R[12].fastq.gz$/;
      }
    }
  }
  if (scalar @combined == 2) {
    return @combined;
  } else {
    return (@return);
  }
}

sub combine {
  # need to know whether to delete files
  my ($lib, $dir, $delete) = @_;
  my @files = q12($lib, $dir);
  my $r1in = "";
  my $r2in = "";
  my $r1out = "$dir/$ARGV[0].h5/$ARGV[0]-combined_R1.fastq.gz";
  my $r2out = "$dir/$ARGV[0].h5/$ARGV[0]-combined_R2.fastq.gz";
  my $i;
  my @return = ($r1out, $r2out);
  open R1O, "| gzip > $r1out";
  open R2O, "| gzip > $r2out";
  foreach $i (0..$#files) {
    if ($files[$i] =~ /_R1\.fastq$/     ||
        $files[$i] =~ /_R1_\d+\.fastq$/ ||
        $files[$i] =~ /_1\.fq$/
       ) {
      open R1I, $files[$i];
      while (<R1I>) {
        print R1O;
      }
      close R1I;
      unlink $files[$i] if $delete;
    } elsif ($files[$i] =~ /_R1\.fastq.gz$/     ||
             $files[$i] =~ /_R1_\d+\.fastq.gz$/ ||
             $files[$i] =~ /_1\.fq\.gz$/
       ) {
      open R1I, "zcat $files[$i] |";
      while (<R1I>) {
        print R1O;
      }
      close R1I;
      unlink $files[$i] if $delete;
    } elsif ($files[$i] =~ /_R2\.fastq$/ ||
        $files[$i] =~ /_R2_\d+\.fastq$/  ||
        $files[$i] =~ /_2\.fq$/
       ) {
      open R2I, $files[$i];
      while (<R2I>) {
        print R2O;
      }
      close R2I;
      unlink $files[$i] if $delete;
    } elsif ($files[$i] =~ /_R2\.fastq.gz$/     ||
             $files[$i] =~ /_R2_\d+\.fastq.gz$/ ||
             $files[$i] =~ /_2\.fq\.gz$/
       ) {
      open R2I, "zcat $files[$i] |";
      while (<R2I>) {
        print R2O;
      }
      close R2I;
      unlink $files[$i] if $delete;
    }
  }
  close R1O;
  close R2O;
  if (-f $r1out && -f $r2out) {
    return @return;
  } else {
    die "Some problem, no $r1out or $r2out in combine...\n";
  }
}

sub make_url {
  my ($run, $exp) = @_;
  return undef if !defined $run || !length $run || !defined $exp || !length $exp;
  my $url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/";
  $url .= substr($exp, 0, 3) . "/" . substr($exp, 0, 6) . "/" . $exp . "/";
  $url .= $run . "/" . $run . ".sra";
  return($url);
}
