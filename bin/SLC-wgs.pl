#!/usr/bin/perl
#
# wrapper to run a whole WGS pipeline
# - species classification (kraken2)
# - assembly (velvetoptimizer or SPAdes)
# - scaffolding (OPERA, FinIS)
# - gene prediction (prokka)
#
use warnings;
use strict;
use Getopt::Long;
use LWP::Simple;
use File::Temp;
use File::Basename;
use File::Spec;
use File::Copy;
use File::Path;
use File::Which("which");
use Sys::MemInfo;
use Bio::SeqIO;
use Cwd;
use JSON;

# programs we need
my $bin_dir = "/usr/local/bin/";
my $programs;
my $ex;
my $original_commandline = join (" ", $0, @ARGV);
my @arguments = @ARGV;
# required programs first
$programs->{prokka} = which 'prokka';
$programs->{seqtk} = which 'seqtk';
foreach $ex (keys %$programs) {
  if (!defined $programs->{$ex} && !-f $programs->{$ex}) {
    print "Can't find required program/library $ex, expected at $programs->{$ex}.  Exiting...\n";
    die;
  }
}
# we need at least one assembler
my @assemblers = qw(skesa spades VelvetOptimiser.pl);
my $have_assembler = 0;
$programs->{skesa} = which 'skesa';
$programs->{spades} = which 'spades.py';
opendir VELVET, "/usr/local/src";
while (my $v = readdir VELVET) {
  if (-d "/usr/local/src/$v" && $v =~ /velvet/ &&
      -f "/usr/local/src/$v/contrib/shuffleSequences_fasta/shuffleSequences_fastq.pl" &&
      -f "/usr/local/src/$v/contrib/VelvetOptimiser-2.2.4/VelvetOptimiser.pl") {
    $programs->{VelvetShuffle} = "/usr/local/src/$v/contrib/shuffleSequences_fasta/shuffleSequences_fastq.pl";
    $programs->{VelvetOptimiser} = "/usr/local/src/$v/contrib/VelvetOptimiser-2.2.4/VelvetOptimiser.pl";
    last;
  }
}
closedir VELVET;
$programs->{velvetg} = which 'velvetg';
foreach $ex (@assemblers) {
  if (defined $programs->{$ex} && -f $programs->{$ex}) {
    $have_assembler = 1;
    last;
  }
}
if (!$have_assembler) {
  print "Can't find any assembler! Exiting...\n";
  die;
}
# optional programs
$programs->{kraken2} = which 'kraken2';
$programs->{fastqc} = which 'fastqc';
$programs->{OPERA} = which 'OPERA-LG';
$programs->{OPERA_preprocess} = which 'preprocess_reads.pl';
$programs->{FinIS} = which 'FinIS';

# parameters
my $tempdir;
my $final_sample_name = "final_output";
my $libname;
my $q1;
my $q2;
my $cleanup = 1;
my $temp_base = "/tmp";
my $output_dir = "/tmp";
my @output_whitelist = (
  "assembly",
  "err",
  "faa",
  "ffn",
  "fna",
  "gbk",
  "gff"
);

# These option hashes are going to be global...
my $default_options = {};
$default_options->{VELVETOPTIMIZER}->{s} = 0;	# 0 is auto
$default_options->{VELVETOPTIMIZER}->{e} = 99;	# fallback if auto fails
$default_options->{FASTQC}->{dofastqc} = 0;
$default_options->{threads} = `nproc`;
chomp $default_options->{threads};
if (!$default_options->{threads}) { $default_options->{threads} = 1; }
$default_options->{ASSEMBLER} = "velvetoptimizer";
$default_options->{SCAFFOLDER} = "opera";
$default_options->{OPERA}->{mapper} = "bwa";
$default_options->{contig_min} = 500;
$default_options->{ANNOTATOR} = "prokka";
$default_options->{MAXREADS} = 5000000;
$default_options->{MEMORYINGB} = int(Sys::MemInfo::totalmem()/1024/1024/1024);
chomp($default_options->{MEMORYINGB});
if (!$default_options->{MEMORYINGB}) {
  # this is just a guess for what a "reasonable" system should have
  $default_options->{MEMORYINGB} = 2;
} elsif ($default_options->{MEMORYINGB} >= 6) {
  $default_options->{MEMORYINGB} -= 2;
} elsif ($default_options->{MEMORYINGB} < 6) {
  $default_options->{MEMORYINGB} -= 1;
} if ($default_options->{MEMORYINGB} < 2) {
  $default_options->{MEMORYINGB} = 1;
}
$default_options->{SEED} = 11;
$default_options->{KRAKEN2}->{reads} = 100000;
$default_options->{KRAKEN2}->{DB} = "/usr/local/lib/Kraken2/minikraken2_v2_8GB_201904_UPDATE";
$default_options->{KRAKEN2}->{classification} = "NONE";

my $user_options = {};

# variables
my $coords;
my $rank;
my $novel;
my @reference;
my $reference;
my $delta;
my $snps;
my $genes;
my $contigs;
my $ordering;
my $assembly;
my @results;
my $index_dir;
my $current_dir;
my $lane;
my $pttfile;
my $ptt_ref;
my $best_reference_sequence;
my $mux_data;
my $mux;
my $plex;
my $global_threads = 0;
my $kraken_classification = "";
my $preferred_reference = "";
my $tempparent;

GetOptions (
  'tempdir=s' => \$temp_base,
  'q1=s' => \$q1,
  'q2=s' => \$q2,
  'name=s' => \$final_sample_name,
  't=i' => \$global_threads,
  'velvetoptimizer=s%' => \$user_options->{VELVETOPTIMIZER},
  'fastqc=s%' => \$user_options->{FASTQC},
  'qc!' => \$user_options->{FASTQC}->{dofastqc},
  'contig_min=i' => \$user_options->{contig_min},
  'cleanup!' => \$cleanup,
  'output_dir=s' => \$output_dir,
  'reference=s' => \$user_options->{REFERENCE},
  'assembler=s' => \$user_options->{ASSEMBLER},
  'scaffolder=s' => \$user_options->{SCAFFOLDER},
  'annotator=s' => \$user_options->{ANNOTATOR},
  'maxreads=i' => \$user_options->{MAXREADS},
  'memory=i' => \$user_options->{MEMORYINGB},
  'seed=i' => \$user_options->{SEED},
  'krakenreads=i' => \$user_options->{KRAKEN2}->{reads},
  'species=s' => \$user_options->{KRAKEN2}->{classification}
);

# postprocess some options

if ($global_threads) {
  $user_options->{threads} = $global_threads;
}

# make sure we have absolute paths
if (!File::Spec->file_name_is_absolute($temp_base)) {
  $temp_base = File::Spec->rel2abs($temp_base);
}
if (!File::Spec->file_name_is_absolute($output_dir)) {
  $output_dir = File::Spec->rel2abs($output_dir);
}
if (!-d $output_dir) {
  File::Path::make_path($output_dir) || die "Can't find or make output directory $output_dir\n";
}

# figure out how we're getting the files and info
# if we have MUX that's all we need and we go get the info and run it all
# to run things manually, we can also use this info: -q1 -q2 -name -insert -output_dir
# if running manually another common option would be -reference

if (!defined $q1 || !-f $q1 || (defined $q2 && !-f $q2) || !length($final_sample_name)) {
    print "Can't find fastq file(s) or have no sample name (-name option), exiting...\n";
    &print_usage;
    exit;
}


####################
# Do all the stuff #
####################

# easier with absolute filenames and paths, we'll check $q1 and $q2 later
$current_dir = getcwd;
$output_dir = File::Spec->rel2abs($output_dir);

$final_sample_name = clean_final_name($final_sample_name);

# set up files and directories we need
$user_options->{tempparent} = File::Temp::tempdir('WGS-XXXXXX', DIR => $temp_base, CLEANUP => $cleanup );
$user_options->{tempdir} = "$user_options->{tempparent}/$final_sample_name";
$tempdir = set_option("tempdir");
mkdir($tempdir);
open LOG, ">$tempdir/$final_sample_name-log.txt";
open SHORTLOG, ">$tempdir/$final_sample_name-summary.txt";
&shortlog("Your original command line was $original_commandline\n");

&shortlog(do_single_library($final_sample_name, $q1, $q2));
&post_process($final_sample_name);

# return to our original place so we can clean up temp directories if needed
chdir($current_dir);


#########################
# Subroutines down here #
#########################

sub clean_final_name {
  my ($name) = @_;
  if (defined $name && length $name) {
    $name =~ s/\s+/_/g;
    $name =~ s/\//_/g;
    $name =~ s/\\/_/g;
    $name =~ s/\?/_/g;
    $name =~ s/\!/_/g;
    $name =~ s/\(/_/g;
    $name =~ s/\)/_/g;
    $name =~ s/'/_/g;
    $name =~ s/"/_/g;
  } else {
    $name = "";
  }
  return ($name);
}

sub assemble_skesa {
  # we need fastq sequences
  my ($q1, $q2) = @_;
  my $tempdir = set_option("tempdir");
  my $expected_file = "$tempdir/$final_sample_name.skesa.fasta";
  my $command = "$programs->{skesa}";
  my @out;
  my $in;
  my $head;
  my $s;
  if ($q1 eq $q2) {
    # single-end
    $command .= " --reads $q1";
  } else {
    $command .= " --reads $q1,$q2";
  }
  $command .= " --cores " . set_option("threads");
  $command .= " --memory " . set_option("MEMORYINGB");
  $command .= " --min_contig ". set_option("contig_min");
  $command .= " > $expected_file";
  &shortlog($command);
  my $output = `$command 2>&1`;
  &log($output);
  chdir($tempdir);
  if (-f $expected_file) {
    return($expected_file);
  } else {
    &log("Couldn't find $expected_file after ".(caller(0))[3].", SPAdes run");
    return (0);
  }
}

sub assemble_SPAdes {
  # we need fastq sequences
  my ($q1, $q2) = @_;
  my $tempdir = set_option("tempdir");
  my $expected_file = "$tempdir/$final_sample_name.spades.fasta";
  my $spades_scaffolds = "$tempdir/spades_working/scaffolds.fasta";
  my $command = "$programs->{spades}";
  my @out;
  my $in;
  my $head;
  my $s;
  if ($q1 eq $q2) {
    # single-end
    $command .= " --s1 $q1";
  } else {
    $command .= " --pe1-1 $q1 --pe1-2 $q2";
  }
  $command .= " -t " . set_option("threads");
  $command .= " -m " . set_option("MEMORYINGB");
  $command .= " --cov-cutoff auto --isolate -o $tempdir/spades_working";
  &shortlog($command);
  my $output = `$command 2>&1`;
  &log($output);
  chdir($tempdir);
  if (-f $spades_scaffolds) {
    # we have to do our own length filtering for SPAdes assembly
    open SPADES_SCAF, $spades_scaffolds;
    open SPADES_FILT, ">$expected_file";
    @out = ();
    $head = "";
    $s = "";
    while ($in = <SPADES_SCAF>) {
      chomp $in;
      if ($in =~ /^>/) {
        if ($head ne "" && length($s) > set_option("contig_min")) {
          print SPADES_FILT "$head\n";
          print SPADES_FILT "$s\n";
        }
        $head = $in;
        $s = "";
      } else {
        $s .= $in;
      }
    }
    if ($head ne "" && length($s) > set_option("contig_min")) {
      print SPADES_FILT "$head\n";
      print SPADES_FILT "$s\n";
    }
    close SPADES_FILT;
    close SPADES_SCAF;
    return($expected_file);
  } else {
    &log("Couldn't find $spades_scaffolds after ".(caller(0))[3].", SPAdes run");
    return (0);
  }
}

sub assemble_velvet {
  # use velvet optimiser
  # options should only be start and end kmer to use
  my ($q1, $q2) = @_;
  my $tempdir = set_option("tempdir");
  my $kstart = set_option("VELVETOPTIMIZER", "s");
  my $kend = set_option("VELVETOPTIMIZER", "e");
  my $readlength = `head -n 4 $q1 | grep -A1 '^\@' | tail -n 1 | wc | awk '{printf "\%.0f", \$3 - 1}'`;
  my $best_middle = 0;
  my $velvetoptimizer_output_dir = "$tempdir/velvetoptimiser";
  my $velvetlog = "";
  my $vlog;
  my $velvetg_options;
  my $velveth_k;
  my $single_fastq;
  my $file_spec;
  my $shuffle;
  my $expected_file;
  my $command;
  my $velvet_output_contigs;
  my $output;
  # we need to make a single file with all the reads for velvet
  if ($q1 eq $q2) {
    $single_fastq = 1;
    $file_spec = "'-short -fastq $q1'";
  } else {
    $single_fastq = 0;
    $shuffle = "$tempdir/$final_sample_name.shuffle-fastq.txt";
    $expected_file = $shuffle;
    $command = "$programs->{VelvetShuffle} $q1 $q2 $shuffle";
    &shortlog($command);
    $output = `$command 2>&1`;
    &log($output);
    if (!-f $expected_file) {
      &log("Couldn't find $expected_file after ".(caller(0))[3].", velvet shuffle run");
      return (0);
    }
    $file_spec = "'-shortPaired -fastq $shuffle'";
  }
  # VelvetOptimizer with multiple threads tries to use OpenMP which isn't available on the SMPs.  Force threads to be 1 (-t option) for now
  if ($kstart == 0) {
    # try do to an optimization
    # just guess the kmer should be around 70% of the read length - and odd
    # try to adjust up or down appropriately
    # do 10 kmer values at a time
    # max kmer should be 255, but let's say reasonable starting point no matter
    # what the read length is 99, and a reasonable ending point is 31
    # when we find the right one we just set kstart and kend to the optimal
    # and run the original version
    $kstart = int(($readlength * 0.7) / 2) * 2 + 1 - 10;
    $kstart = 99 if $kstart > 99;
    &shortlog("Beginning auto optimization of velvet kmer");
    while (!$best_middle) {
      $kend = $kstart + 20;
      $kend = 255 if $kend > 255;
      $command = "$programs->{VelvetOptimiser} -s $kstart -e $kend -d $tempdir/velvetoptimiser -f $file_spec -t 1 --o '-min_contig_lgth " . set_option("contig_min") . "'";
      &shortlog("Auto optimization of velvet kmer: $command");
      $output = `$command 2>&1`;
      &log($output);
      opendir DIR, $velvetoptimizer_output_dir;
      while ($velvetlog = readdir(DIR)) {
        last if $velvetlog =~ /Logfile.txt$/;
      }
      closedir DIR;
      if (-f "$velvetoptimizer_output_dir/$velvetlog") {
        open VLOG, "$velvetoptimizer_output_dir/$velvetlog";
        while ($vlog = <VLOG>) {
          # the final options should be the last one listed
          if ($vlog =~ /^Velvetg parameter string: (.*)$/) {
            $velvetg_options = $1;
          }
        }
        close VLOG;
        if ($velvetg_options =~ /auto_data_(\d+)/) {
          $velveth_k = $1;
        }
        # stop if we hit our limits or are done
        if ($velveth_k == 31 || $velveth_k == 255 ||
            ($kstart < $velveth_k && $velveth_k < $kend) ) {
          $best_middle = 1;
          $kstart = $velveth_k;
          $kend = $velveth_k;
          &shortlog("Best kmer: $velveth_k");
        } else {
          if ($velveth_k == $kstart) {
            $kstart -= 18;	# don't move the full 20 so we overlap by 2,
				# in case the optimal is the edge of one of our
				# intervals
            $kstart = 31 if $kstart < 31;
          } elsif ($velveth_k == $kend) {
            $kstart += 18;
            $kstart = 255 if $kstart > 255;	# though this shouldn't happen
          }
        }
      } else {
        # we failed at this velvetoptimiser run for some reason...
        $best_middle = 1;
        $kend = set_option("VELVETOPTIMIZER", "e");
        $kstart = $kend - 20;
        &log("Couldn't auto-optimize velvet assembly, reverting to range $kstart to $kend for kmer");
      }
      # clean things up - space is one of the reasons to do this optimization
      # VelvetOptimiser should already clean up the auto_data dirs
      if (-d $velvetoptimizer_output_dir) {
        File::Path::remove_tree($velvetoptimizer_output_dir);
      }
    }
  }

  # the final VelvetOptimiser run
  $command = "$programs->{VelvetOptimiser} -s $kstart -e $kend -d $tempdir/velvetoptimiser -f $file_spec -t 1 --o '-min_contig_lgth " . set_option("contig_min") . "'";
  $velvet_output_contigs = "$velvetoptimizer_output_dir/contigs.fa";
  $expected_file = "$tempdir/$final_sample_name.contigs.txt";
  &shortlog($command);
  $output = `$command 2>&1`;
  &log($output);
  # VelvetOptimizer runs velvetg with the -clean yes switch.  This deletes the
  # LastGraph file that's needed for FinIS finishing...
  # the last velvetg command is in the log file in that directory
  # velvetg on the final directory
  $velvetg_options = "-clean no";
  $velveth_k = 0;
  my $velvetg_exp_cov = 0;
  my $velvetg_cov_cutoff = 0;
  opendir DIR, $velvetoptimizer_output_dir;
  while ($velvetlog = readdir(DIR)) {
    last if $velvetlog =~ /Logfile.txt$/;
  }
  closedir DIR;
  if (-f "$velvetoptimizer_output_dir/$velvetlog") {
    open VLOG, "$velvetoptimizer_output_dir/$velvetlog";
    while ($vlog = <VLOG>) {
      # the final options should be the last one listed
      if ($vlog =~ /^Velvetg parameter string: (.*)$/) {
        $velvetg_options = $1;
      }
    }
    close VLOG;
    if ($velvetg_options =~ /auto_data_(\d+)/) {
      $velveth_k = $1;
    }
    if ($velvetg_options =~ /-exp_cov\s+(\d+)/) {
      $velvetg_exp_cov = $1;
    }
    if ($velvetg_options =~ /-cov_cutoff\s+(\d+\.?\d*)/) {
      $velvetg_cov_cutoff = $1;
    }
    $velvetg_options =~ s/-clean yes/-clean no/;
    $velvetg_options =~ s/auto_data_\d+//;
  }

  $command = "$programs->{velvetg} $velvetoptimizer_output_dir $velvetg_options";
  &shortlog("velveth kmer: $velveth_k");
  &shortlog("velvetg exp_cov: $velvetg_exp_cov");
  &shortlog("velvetg cov_cutoff: $velvetg_cov_cutoff");
  &shortlog($command);
  $output = `$command 2>&1`;
  &log($output);
  
  # this gets put into $tempdir/velvetoptimiser/
  if (-f $velvet_output_contigs) {
    File::Copy::copy("$velvet_output_contigs", $expected_file);
    if (!$single_fastq && defined $shuffle && -f $shuffle) {
      unlink($shuffle);
    }
    return ($expected_file);
  } else {
    &log("Couldn't find $velvet_output_contigs after ".(caller(0))[3].", velvet optimiser run");
    if (!$single_fastq && defined $shuffle && -f $shuffle) {
      unlink($shuffle);
    }
    return (0);
  }
}

sub scaffold_OPERA {
  my ($q1, $q2, $contigs) = @_;
  my $tempdir = set_option("tempdir");
  my $opera_output_dir = "$tempdir/OPERA";
  my $mapper = set_option("OPERA", "mapper");
  my $expected_file = "$tempdir/$final_sample_name.gapFilled";
  my $map_file = "$tempdir/$final_sample_name.read-on-contig.bam";

  # OPERA requires samtools 0.1.19 or below (see wiki)
  # look for these here
  my $samtools_dir = "";
  if (-f "/usr/local/src/samtools-0.1.19/samtools") {
    $samtools_dir = "/usr/local/src/samtools-0.1.19";
  } elsif (-f "/usr/local/src/samtools-0.1.18/samtools") {
    $samtools_dir = "/usr/local/src/samtools-0.1.18";
  } elsif (-f "/usr/local/src/samtools-0.1.17/samtools") {
    $samtools_dir = "/usr/local/src/samtools-0.1.17";
  }

  if ($q1 eq $q2) {
    File::Copy::copy($contigs, $expected_file);
    &shortlog("q1 and q2 files the same.  Skipping scaffolding.");
  } else {
    # OPERA makes an assumption about where perl is that's not true for us
    my $command = "perl $programs->{OPERA_preprocess} --contig $contigs --illumina-read1 $q1 --illumina-read2 $q2 --out $map_file --map-tool $mapper";
    if ($samtools_dir ne "") {
      $command .= " --samtools-dir $samtools_dir";
    }
    &shortlog($command);
    my $output = `$command 2>&1`;
    # OPERA's preprocess leaves this read.sai file
    unlink("$tempdir/read.sai");
    &log($output);
    if (!-f $map_file) {
      &log("Couldn't find mapping file $map_file after ".(caller(0))[3]." run");
      return(0);
    }
    $command = "$programs->{OPERA} $contigs $map_file $opera_output_dir $samtools_dir";
    &shortlog($command);
    $output = `$command 2>&1`;
    &log($output);
    if (-f "$opera_output_dir/scaffoldSeq.fasta") {
      File::Copy::copy("$opera_output_dir/scaffoldSeq.fasta", $expected_file);
    } else {
      &log("Couldn't find $expected_file after ".(caller(0))[3]." run");
      return(0);
    }
  }

  return($expected_file);
}

sub finish_FinIS {
  my ($velvet_folder, $scaffolds, $output_dir) = @_;
  my $expected_file;
  my $command = "$programs->{FinIS} $velvet_folder $scaffolds $output_dir";
  my $output;
  my $finis_file = "$output_dir/scaffolds.filled.fasta";
  # check for what we need
  if ((-f "$velvet_folder/contigs.fa" && -f "$velvet_folder/LastGraph") ||
      (-f "$velvet_folder/$final_sample_name.scafSeq")
     ) {
    mkdir($output_dir);
    $expected_file = "$tempdir/$final_sample_name.gapFilled.FinIS";
    &shortlog($command);
    $output = `$command 2>&1`;
    &log($output);
    if (-f $finis_file) {
      File::Copy::copy($finis_file, $expected_file);
    } else {
      &log("Couldn't find $finis_file after ".(caller(0))[3]." run");
      return(0);
    }
    return($expected_file);
  } else {
    &shortlog("Couldn't find contigs.fa or LastGraph file in $velvet_folder");
    return 0;
  }
}

sub species_kraken2 {
  my ($q1, $q2) = @_;
  my $expected_file = "kraken2-out.report";
  my $output_file = "kraken2-out.tmp";
  my $command = "$programs->{kraken2}";
  my $kraken_ds1_filename;
  my $kraken_ds2_filename;
  my $kraken_dscommand;
  my @f;
  my @g;
  my $i;
  my $max;
  my $backup_max;
  my $output;
  my $return;
  my $return_index;
  my $backup_return;
  my $backup_return_index;
  my $n = set_option("KRAKEN2", "reads");

  # just return if this was set on the command line
  $return = set_option("KRAKEN2", "classification");
  if (defined $return && uc($return) ne "NONE" && $return ne "") {
    return $return;
  }

  $command .= " --threads " . set_option("threads");
  $command .= " --db " . set_option("KRAKEN2", "DB");
  $command .= " --report $expected_file";
  $kraken_ds1_filename = $q1 . "-" . $n . "_" . set_option("SEED") . ".fastq";
  $kraken_ds2_filename = $q2 . "-" . $n . "_" . set_option("SEED") . ".fastq";
  $kraken_dscommand = $programs->{seqtk} . " sample -s " . set_option("SEED") . " $q1 " . $n . " > $kraken_ds1_filename";
  &shortlog($kraken_dscommand);
  $output = `$kraken_dscommand 2>&1`;
  if ($q1 eq $q2) {
    $command .= " $kraken_ds1_filename 2>&1 > $output_file";
  } else {
    $kraken_dscommand = $programs->{seqtk} . " sample -s " . set_option("SEED") . " $q2 " . $n . " > $kraken_ds2_filename";
    &shortlog($kraken_dscommand);
    $output = `$kraken_dscommand 2>&1`;
    $command .= " --paired $kraken_ds1_filename $kraken_ds2_filename 2>&1 > $output_file";
  }
  &shortlog($command);
  # for the kraken run, we're using the trick of 'command 2>&1 > output'
  # this puts stdout into the output file, and stderr shows up on stdout
  # so we can capture stderr with backticks
  # this is because we want the classifications into the output file
  # and the summary comes on stderr which we want to parse
  $output = `$command`;
  &log($output);
  $output =~ s/\r/\n/g;
  @f = split /\n/, $output;
  foreach $i (@f) {
    if ($i =~ /(\d+) sequences .* processed/) {
      &shortlog("Kraken2: $1 total sequences processed");
    } elsif ($i =~ /^\s*\d+ sequences classified/) {
      $i =~ s/^\s+//;
      &shortlog("Kraken2: $i");
    } elsif ($i =~ /^\s*\d+ sequences unclassified/) {
      $i =~ s/^\s+//;
      &shortlog("Kraken2: $i");
    }
  }
  if (-f $expected_file) {
    open REPORT, $expected_file;
    @f = <REPORT>;
    close REPORT;
    $max = 0;
    $backup_max = 0;
    $return = "";
    $return_index = -1;
    $backup_return = "";
    $backup_return_index = -1;
    foreach $i (0..$#f) {
      $f[$i] =~ s/^\s+//;
      @g = split /\t/, $f[$i];
      if ($g[3] eq "G") {
        if ($g[0] > $backup_max) {
          $backup_max = $g[0];
          $backup_return = $g[5];
          $backup_return =~ s/^\s+//;
          $backup_return_index = $i;
        }
      }
      if ($g[3] eq "S") {
        if ($g[0] > $max) {
          $max = $g[0];
          $return = $g[5];
          $return =~ s/^\s+//;
          $return_index = $i;
        }
      }
    }
    &log("Cleaning up kraken intermediate files $output_file, $expected_file, and downsampled inputs...");
    unlink($output_file);
    unlink($expected_file);
    unlink($kraken_ds1_filename);
    unlink($kraken_ds2_filename);
    if ($return_index >= 0) {
      &shortlog("Kraken classification line: $f[$return_index]");
      return($return);
    } elsif ($backup_return_index >= 0) {
      &shortlog("Kraken classification line: $f[$backup_return_index]");
      return($backup_return);
    }
  } else {
    &shortlog("Couldn't find $expected_file file after initial Kraken2 run, skipping...");
    return("");
  }
  return("");
}

sub select_reference {
  my ($assembly, $num_to_test, $preferred_reference, $initial) = @_;
  if (!defined $num_to_test || $num_to_test <= 0) {
    $num_to_test = 200;
  }
  my $command = "$programs->{closest_species} -num_test $num_to_test $assembly";
  if ($initial ne "") {
    $command .= " -initial $initial";
  }
  my $output = `$command`;
  my @reference = split /\n/, $output;
  # if there was a reference specified, we'll use that as the "best"
  # move it around if it's already in the list, if not then just stick it
  # at the front
  if (defined $preferred_reference && length $preferred_reference && -f $preferred_reference) {
    unshift @reference, $preferred_reference;
    foreach my $i (1..$#reference) {
      if ($preferred_reference eq $reference[$i]) {
        splice @reference, $i, 1;
        last;
      }
    }
  }
  return (@reference);
}

sub annotate_prokka {
  my ($assembly, $prefix, $refspecies) = @_;
  my $genus;
  my $species;
  my @wanted = qw(err faa ffn fna gbk gff);
  my $outdir = "$tempdir/prokka";
  my $species_flag;
  my $command;
  my $output;
  my $i;
  my $ext;
  my $error;
  undef $genus;
  undef $species;
#  ($genus, $species) = split /_/, $refspecies;
  ($genus, $species) = split /\s+/, $refspecies;
  if (defined $genus && length $genus > 0) {
    if (defined $species && length $species > 0) {
      $species_flag = "--species $species";
    }
    $command = "$programs->{prokka} --outdir $outdir --prefix $prefix --genus $genus $species_flag --quiet --cpus " . set_option("threads") . " --strain $prefix --kingdom Bacteria --gcode 0 --addgenes --locustag $prefix --mincontiglen " . set_option("contig_min") . " --usegenus $assembly";
  } else {
    # if no genus, use neither genus nor species
    $command = "$programs->{prokka} --outdir $outdir --prefix $prefix --quiet --cpus " . set_option("threads") . " --strain $prefix --kingdom Bacteria --gcode 0 --addgenes --locustag $prefix --mincontiglen " . set_option("contig_min") . " $assembly";
  }
  &shortlog($command);

  # overall output log
  $output = `$command 2>&1`;
  # should be no output here actually - get log from the log file
  if (-f "$outdir/$prefix.log") {
    $output = "";
    open PROKKAF, "$outdir/$prefix.log";
    while ($i = <PROKKAF>) {
      $output .= $i;
    }
    close PROKKAF;
  }
  &log($output);

  # summary of the run
  if (-f "$outdir/$prefix.txt") {
    $output = "\n==== Prokka Summary ====\n";
    open PROKKAF, "$outdir/$prefix.txt";
    while ($i = <PROKKAF>) {
      $output .= $i;
    }
    close PROKKAF;
  }
  &shortlog($output);

  $error = 0;
  foreach $ext (@wanted) {
    if (-f "$outdir/$prefix.$ext") {
      File::Copy::copy("$outdir/$prefix.$ext", "$tempdir/$prefix.$ext");
    } else {
      &log("Couldn't find $prefix.$ext after ".(caller(0))[3]." run");
      $error++;
    }
  }
  if ($cleanup) {
    File::Path::remove_tree($outdir);
  }
  return($error);
}

sub join_assembly {
  my ($assembly, $prefix) = @_;
  my @pf;
  my $pfseqs;
  my @pfout;
  my %used;
  my $key;
  my $expected_file = "$assembly-joined";
  &log("Joining assembly for gene calling\n");
  open PF, $assembly;
  open PS, ">$expected_file";
  print PS ">", substr($prefix,0,20), "\n";
  @pf = <PF>;
  # get rid of spaces for other info like length, just keep contig name
  foreach my $i (0..$#pf) {
    if ($pf[$i] =~ /^>/) {
      my @f = split /\s+/, $pf[$i];
      $pf[$i] = $f[0];
    }
  }
  $pfseqs = fasta2hash(@pf);
  @pfout = ();
  # go from longest to shortest
  foreach $key (sort {length($pfseqs->{$b}) <=> length($pfseqs->{$a})} keys %$pfseqs) {
    push @pfout, $pfseqs->{$key};
  }
  print PS join ("N" x 100, @pfout);
  close PS;
  close PF;
  return ($expected_file) if -f $expected_file;
}

sub qc_FASTQC {
  my (@files) = @_;
  my $tempdir = set_option("tempdir");
  my $f;
  my $command = "$programs->{fastqc} -o $tempdir --noextract";
  $command .= " -t " . set_option("threads");
  $command .= " " . join (" ", @files);
  &shortlog($command);
  my $output = `$command 2>&1`;
  &log($output);
}

sub shortlog {
  my ($message) = @_;
  # this will also do a regular log of everything as well
  chomp $message;
  &log($message);
  print SHORTLOG "$message\n";
}

sub log {
  my ($message) = @_;
  print LOG "# ", basename($0), " ## ", scalar(localtime), "\n";
  chomp $message;
  print LOG "# ", basename($0), " ## $message\n";
}

sub do_single_library {
  my ($final_sample_name, $q1, $q2) = @_;

  my $tempparent = set_option("tempparent");
  my $tempdir = set_option("tempdir");
  my $contigs;
  my $assembly;
  my @reference;
  my ($delta, $coords, $snps);
  my $novel;
  my $genes;
  my $contig_stats;
  my $scaffold_stats;
  my $current_dir;
  my $final_information = "\n==== Final Summary ====\n";
  my $tempinfo;
  my $single_fastq = 0;

  if (defined $q2) {
    print STDERR "STARTED $final_sample_name using files $q1 and $q2 in $tempdir by $ENV{USER} on ", `hostname`;
    &shortlog("STARTED $final_sample_name using files $q1 and $q2 in $tempdir by $ENV{USER} on " . `hostname`);
    &shortlog("Files:\n" . `md5sum $0 $q1 $q2`);
  } else {
    print STDERR "STARTED $final_sample_name using file $q1 in $tempdir by $ENV{USER} on ", `hostname`;
    $q2 = $q1;	# this is how we'll deal with single end internally
    $single_fastq = 1;
    &shortlog("STARTED $final_sample_name using file $q1 in $tempdir by $ENV{USER} on ", `hostname`);
    &shortlog("Files:\n" . `md5sum $0 $q1`);
  }
  print STDERR "You can follow the progress by running 'tail -f $tempdir/$final_sample_name-log.txt'; you may need to be logged in to ", `hostname`;
  &shortlog("You can follow the progress using \'tail -f $tempdir/$final_sample_name-log.txt\; you may need to be logged in to " . `hostname`);

  if (!defined $q1 || !-f $q1 || !defined $q2 || !-f $q2) {
    &log("Can't find file $q1 and/or $q2 for $final_sample_name, skipping...\n");
    close LOG;
    close SHORTLOG;
    return;
  }

  # partly this gives us absolute paths for $q1 and $q2 so we can move around
  $q1 = File::Spec->rel2abs($q1);
  $q2 = File::Spec->rel2abs($q2);
  $current_dir = getcwd;
  chdir($tempdir);
  my $cleanup_fastq1 = 0;
  my $cleanup_fastq2 = 0;
  if ($q1 =~ /\.gz$/) {
    my $origq = $q1;
    $q1 =~ s/\.gz$//;
    $q1 = $tempdir . "/" . basename($q1);
    my $unzip = `zcat $origq > $q1`;
    &shortlog("Unzipping gzipped $origq into temp directory ($tempdir)");
    $cleanup_fastq1 = 1;
  }
  if ($single_fastq) {
    $q2 = $q1;
  } elsif ($q2 =~ /\.gz$/) {
    my $origq = $q2;
    $q2 =~ s/\.gz$//;
    $q2 = $tempdir . "/" . basename($q2);
    my $unzip = `zcat $origq > $q2`;
    &shortlog("Unzipping gzipped $origq into temp directory ($tempdir)");
    $cleanup_fastq2 = 1;
  }

  # collect info
  $tempinfo = `wc -l $q1 | awk '{printf "\%.0f", \$1/4}'`;
  chomp $tempinfo;
  if ($single_fastq) {
    $final_information .= "Number of single-end reads: $tempinfo\n";
  } else {
    $final_information .= "Number of paired-end reads: $tempinfo\n";
  }
  $kraken_classification = species_kraken2($q1, $q2);
  $final_information .= "Kraken classification: $kraken_classification\n";
  if (set_option("MAXREADS") > 0 && set_option("MAXREADS") < $tempinfo) {
    &shortlog("Downsampling to " . set_option("MAXREADS") . " reads with seed " . set_option("SEED") . "\n");
    $final_information .= "Reads used for assembly: " . set_option("MAXREADS") . "\n";
    my $ds_filename = $q1 . "-" . set_option("MAXREADS") . "_" . set_option("SEED");
    my $command = $programs->{seqtk} . " sample -s " . set_option("SEED") . " $q1 " . set_option("MAXREADS") . " > $ds_filename";
    &shortlog($command);
    my $output = `$command`;
    if ($cleanup_fastq1) {
      unlink("$tempdir/".basename($q1));
    }
    $q1 = $ds_filename;
    if ($single_fastq) {
      $q2 = $ds_filename;
    } else {
      $ds_filename = $q2 . "-" . set_option("MAXREADS") . "_" . set_option("SEED");
      $command = $programs->{seqtk} . " sample -s " . set_option("SEED") . " $q2 " . set_option("MAXREADS") . " > $ds_filename";
      &shortlog($command);
      $output = `$command`;
      if ($cleanup_fastq2) {
        unlink("$tempdir/".basename($q2));
      }
      $q2 = $ds_filename;
    }
  }
  # this should be read length +1 because of newline
  $tempinfo = `head -n 4 $q1 | grep -A1 '^\@' | tail -n 1 | wc | awk '{printf "\%.0f", \$3 - 1}'`;
  chomp $tempinfo;
  $final_information .= "Read length (based on first read): $tempinfo\n";

  #
  # FastQC
  #
  if (set_option("FASTQC", "dofastqc")) {
    if ($single_fastq) {
      &qc_FASTQC($q1);
    } else {
      &qc_FASTQC($q1, $q2);
    }
  }

  #
  # assembly
  #
  if (uc(set_option("ASSEMBLER")) eq "VELVETOPTIMIZER" ||
      uc(set_option("ASSEMBLER")) eq "VELVETOPTIMISER") {
     $contigs = assemble_velvet($q1, $q2);
    if (!-f $contigs) {
      $contigs = assemble_SPAdes($q1, $q2);
      $final_information .= "Assembly: SPAdes (Velvet failed)\n";
    } else {
      $final_information .= "Assembly: VelvetOptimizer\n";
    }
  } elsif (uc(set_option("ASSEMBLER")) eq "SPADES") {
    $contigs = assemble_SPAdes($q1, $q2);
    if (!-f $contigs) {
      $contigs = assemble_velvet($q1, $q2);
      $final_information .= "Assembly: VelvetOptimizer (SPAdes failed)\n";
    } else {
      $final_information .= "Assembly: SPAdes\n";
    }
  } elsif (uc(set_option("ASSEMBLER")) eq "SKESA") {
    $contigs = assemble_skesa($q1, $q2);
  }
  $contig_stats = assembly_stats($contigs, set_option("contig_min"));
  &shortlog("\n----- Contig stats for initial assembly -----\n$contig_stats->{text}\n-----");
  if ($contig_stats->{total_length} == 0) {
    &log("Couldn't do initial assembly for $final_sample_name, aborting...");
    return;
  }

  #
  # scaffolding and finishing
  #
  if ($q1 ne $q2) {
    &log("Running scaffolding\n");
    if (uc(set_option("SCAFFOLDER")) eq "OPERA") {
      $assembly = scaffold_OPERA($q1, $q2, $contigs);
      $final_information .= "Scaffolding: OPERA\n";
    }
    #
    # finishing if we used velvet and OPERA
    #
    if (!-f $assembly) {
      $assembly = $contigs;
      &shortlog("Some problem with scaffolder; continuing with contigs.");
    } elsif (uc(set_option("SCAFFOLDER")) eq "OPERA" && uc(set_option("ASSEMBLER")) eq "VELVETOPTIMIZER" || uc(set_option("ASSEMBLER")) eq "VELVETOPTIMISER") {
      &log("Running Finishing\n");
      my $finished = "";
      $finished = finish_FinIS("$tempdir/velvetoptimiser", "$tempdir/OPERA/scaffolds.scaf", "$tempdir/FinIS");
      if (-f $finished) {
        $assembly = $finished;
        $final_information .= "Finishing: FinIS\n";
      } else {
        $final_information .= "Finishing: Failed for some reason\n";
      }
    } else {
      $final_information .= "Finishing: None\n";
    }
  } else {
    &log("Single end.  No scaffolding or finishing.\n");
    $assembly = $contigs;
  }

  $scaffold_stats = assembly_stats($assembly, set_option("contig_min"));
  &shortlog("\n---- Final assembly stats -----\n$scaffold_stats->{text}\n-----");
  if ($scaffold_stats->{total_length} == 0) {
    &log("Got zero total length during scaffolding for $final_sample_name, aborting...");
    return;
  }

  &shortlog("Renaming final assembly to $final_sample_name.assembly\n");
  rename($assembly, "$tempdir/$final_sample_name.assembly");
  $assembly = "$tempdir/$final_sample_name.assembly";
  if ($cleanup) {
    &shortlog("Cleaning up intermediate assembly files - contigs.txt, gapFilled, etc.\n");
    unlink("$tempdir/$final_sample_name.contigs.txt");
    unlink("$tempdir/$final_sample_name.gapFilled");
    unlink("$tempdir/$final_sample_name.spades.fasta");
  }
  if ($cleanup_fastq1) {
    unlink("$tempdir/".basename($q1));
  }
  if ($cleanup_fastq2) {
    unlink("$tempdir/".basename($q2));
  }
  $final_information .= $scaffold_stats->{text};

  #
  # Annotation
  #
  # old way - run prodigal and do a blast-based annotation to best reference
  # this unfortunately doesn't have all the files we might want, though they
  # can be generated - ex. .gff
  # new way - prokka - this by default gives all the ncbi refseq files
  # including .fna, .faa, .ffn, .gff, .gbk
  # pick genus / species from closest_species
  #
  my $joined_assembly = join_assembly($assembly, $final_sample_name);
  if (uc(set_option("ANNOTATOR")) eq "PROKKA") {
    # we need the directory which has the genus/species name as the refspecies
    &annotate_prokka($joined_assembly, $final_sample_name, $kraken_classification);
  }
  unlink($joined_assembly);

  chdir($current_dir);
  return($final_information);
}

sub post_process {
  my ($final_sample_name) = @_;
  my $tempdir = set_option("tempdir");
  my $tempparent = set_option("tempparent");
  my $f;
  my $test;
  my $keep;
  if (!$cleanup) {
    print STDERR "NOTE -- Intermediate files have been left in:\n$tempdir\n";
    &shortlog("NOTE -- Intermediate files have been left in:\n$tempdir\n");
  } else {
    opendir CLEANUP, $tempdir;
    while ($f = readdir CLEANUP) {
      if (-f "$tempdir/$f") {
        $keep = 0;
        foreach $test (@output_whitelist) {
          if ($f =~ /\.$test$/) {
            $keep = 1;
            last;
          }
        }
        if ($f =~ /-summary.txt$/ || $f =~ /-log.txt$/) {
          $keep = 1;
        }
        unlink("$tempdir/$f") if !$keep;
      } elsif (-d "$tempdir/$f" && $f ne "." && $f ne "..") {
        File::Path::remove_tree("$tempdir/$f");
      }
    }
  }
  &log("tar cvzf $output_dir/$final_sample_name.tgz $final_sample_name");
  close LOG;
  close SHORTLOG;
  chdir($tempparent);
  `tar cvzf $output_dir/$final_sample_name.tgz $final_sample_name`;
  if (-f "$output_dir/$final_sample_name.tgz" && -s "$output_dir/$final_sample_name.tgz") {
    print STDERR "FINISHED Final output is in $output_dir/$final_sample_name.tgz\n";
    # we don't shortlog here because the file is closed
  } else {
    if (-f "$output_dir/$final_sample_name.tgz") {
      print STDERR "Error - somehow didn't create final output file $output_dir/$final_sample_name.tgz\n";
      # we don't shortlog here because the file is closed
    } else {
      print STDERR "Error - final output file $output_dir/$final_sample_name.tgz has length 0!\n";
      # we don't shortlog here because the file is closed
    }
  }
  return ("$output_dir/$final_sample_name.tgz");
}

sub assembly_stats {
  # we will delete N's in the assembly and not count them - scaffold style
  my $fasta = shift;
  my $contig_cutoff;
  if (scalar @_) {
    $contig_cutoff = shift;
  } else {
    $contig_cutoff = 0;
  }

  my $key;
  my $seq;
  my @a;
  my $subtotal;
  my $total;
  my $return = ();
  my @length = ();

  $return->{number_contigs} = 0;
  $return->{total_length} = 0;
  $return->{avg_length} = 0;
  $return->{max_length} = 0;
  $return->{min_length} = 0;
  $return->{N50_length} = 0;
  $return->{N50_number} = 0;
  $return->{N90_length} = 0;
  $return->{N90_number} = 0;
  $return->{text} = "";

  if (-f $fasta) {
    open FASTA, $fasta;

    # length filter first
    @a = <FASTA>;
    close FASTA;
    $seq = fasta2hash(@a);
    $total = 0;
    foreach $key (keys %$seq) {
      $seq->{$key} =~ tr/nN//d;
      if (length($seq->{$key}) >= $contig_cutoff) {
        push @length, length($seq->{$key});
        $total += length($seq->{$key});
      }
    }
    @length = sort { $a <=> $b } @length;
    $subtotal = 0;
    foreach my $i (reverse 0..$#length) {
      $subtotal += $length[$i];
      if ($subtotal >= $total * 0.9 && $return->{N90_number} == 0) {
        $return->{N90_number} = $#length - $i + 1;
        $return->{N90_length} = $length[$i];
      }
      if ($subtotal >= $total/2 && $return->{N50_number} == 0) {
        $return->{N50_number} = $#length - $i + 1;
        $return->{N50_length} = $length[$i];
      }
    }
    $return->{number_contigs} = scalar(@length);
    $return->{total_length} = $total;
    $return->{avg_length} = int($total/scalar(@length)) if scalar(@length) > 0;
    $return->{max_length} = $length[$#length];
    $return->{min_length} = $length[0];
    foreach $key (qw(number_contigs total_length avg_length max_length min_length N50_length N50_number N90_length N90_number)) {
      $return->{text} .= "$key: $return->{$key}\n";
    }
  }

  return $return;
}

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
#      undef $name;
    }
  }
  return $seq;
}

sub hash2fasta {
  # returns an array
  # doesn't have newlines
  my ($seq) = @_;
  my $name;
  my @return;
  foreach $name (keys %$seq) {
    push @return, ">$name";
    push @return, $seq->{$name};
  }
  return @return;
}

sub revcomp {
  my ($inline) = $_[0];
  my $outline = reverse ($inline);
  $outline =~ tr/ABCDGHKMNRSTVWXYabcdghkmnrstvwxy/TVGHCDMKNYSABWXRtvghcdmknysabwxr/;
  return $outline;
}

sub set_option {
  my (@key) = @_;
  my $return = "";
  my $done = 0;
  my $ref = $user_options;
  foreach my $i (0..$#key) {
    if (defined $ref->{$key[$i]}) {
      $ref = $ref->{$key[$i]};
      if ($i == $#key) {
        $done = 1;
        $return = $ref;
      }
    }
  }
  if (!$done) {
    $ref = $default_options;
    foreach my $i (0..$#key) {
      if (defined $ref->{$key[$i]}) {
        $ref = $ref->{$key[$i]};
        if ($i == $#key) {
          $done = 1;
          $return = $ref;
        }
      }
    }
  }
  return $return;
}

sub print_usage {
  my $name = basename($0);
  print <<__USAGE__;
Usage example:
---------------
   $name -q1 <fastq file> -q2 <fastq file> -name <output file name> -output_dir <output directory> -tempdir <temp directory>

This will give you a file named <output file name>.tgz in the directory <output directory>.  If you choose velvetoptimizer and opera, it will run FinIS for you also.

Common options:
---------------
  -output_dir <dir>     --  where output files go - you probably always want
                              to specify this
  -tempdir              --  temporary directory to use -
                              default $temp_base
  -assembler <velvetoptimizer|spades|skesa>
                        --  which assembler to use - default $default_options->{ASSEMBLER}
  -scaffolder <opera>   --  which scaffolder - default $default_options->{SCAFFOLDER}
  -annotator <prokka>   --  how to annotate - default $default_options->{ANNOTATOR}
  -reference <fna file> --  force a certain reference sequence - this assumes
                              NCBI refseq files are available - ptt, faa, ffn
  -t <threads>          --  number of CPUs to use - default $default_options->{threads}
  -memory <int>         --  memory to specify (in GB) - default $default_options->{MEMORYINGB}
  -contig_min <int>     --  contig minimum cutoff - default $default_options->{contig_min}
  -cleanup|nocleanup    --  whether to keep intermediate files - default cleanup
  -qc|noqc              --  whether to run FastQC (default -noqc)
  -maxreads             --  max number of reads to use for assembly -
                              default $default_options->{MAXREADS}
                              if set to 0, then use everything
  -seed                 --  seed used for seqtk to downsample if we need to
                              use fewer reads
                              default $default_options->{SEED}
  -species <species>    --  If set, won't rerun Kraken 2 to get a species
                              identification
                              default $default_options->{KRAKEN2}->{classification} - meaning run Kraken 2

Advanced options:
-----------------------------------
  -velvetoptimizer option=value

    For each of these, specify the exact option for that program. This script does no checking on these options/values. If you need to specify more than one option, use multiple switches. For example, the defaults for Velvetoptimizer are currently set as follows:

    -velvetoptimizer s=51 -velvetoptimizer e=99

    N.B. For velvetoptimizer there are special options for start and end. The default start is 0, which indicates to do auto-selection and optimization of the k-mer based on the read length. This also will try to use less temp drive space by running velvetoptimizer multiple times if needed. If this fails, the fallback is to just take whatever comes out of velvetoptimizer based on the value of the end parameter (default 99; the start value will be 20 less than this). (i.e. default behavior is --velvetoptimizer s=0 -velvetoptimizer e=99)
__USAGE__
}
