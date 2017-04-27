#!/usr/bin/perl -w
#
# run reticulate program
# import into reticulate table
#
use warnings;
use strict;
use File::Temp qw(tempdir);
use File::Copy;
use File::Basename;
use slchen;

my $currentdir;
my $tempdir;
my $tempfile;
my $output_basename;

my $reticulate_p_value;
my @output;

my $run_silent = 1;	# by default run silent
my $silent_switch;
my @sum;
my @frags;
my $out;
my $inner_p;
my $outer_p;
my $global_inner;
my $global_outer;
my $sequence;
my $sequence_length;
my @break;
my $alignment;
my $fragment_outputfile;

my $name;
my @s;
my $i;
my @f;

use Getopt::Long;
&Getopt::Long::Configure("pass_through");
GetOptions (
  "silent!" => \$run_silent	# whether to use /r switch for GENECONV
);

if (!@ARGV || !-f $ARGV[0]) {
  print "Usage: $0 <fasta file>\n";
  print "Output will be reticulate and sawyer output files and p-values\n";
  exit;
}

$tempdir = tempdir( CLEANUP => 1);
$tempfile = "rs_seqinput";
File::Copy::copy($ARGV[0], "$tempdir/$tempfile");
$currentdir = `pwd`;
chomp $currentdir;
$output_basename = File::Basename::basename($ARGV[0]);

##################
# run reticulate #
##################

chdir $tempdir;
open T, ">input";
print T join ("\n", $tempfile, 'n', $tempfile, '3', 'n', 'p', '', $tempfile, $tempfile, 'r', '1', 'y', '10000', '0', "$tempfile.out", 'q');
close T;
system "reticulate-noX < input > /dev/null";
if (!-f "$tempfile.out") {
  print "Error running reticulate, skipping to GENECONV...\n";
} else {
  open T, "$tempfile.out";
  @output = <T>;
  close T;

  # parse output file for p value
  foreach $i (0..$#output) {
    if ($output[$i] =~ /^The P value is (.*)/) {
      $reticulate_p_value = $1;
      chomp $reticulate_p_value;
    }
  }
  print "Reticulate P-value: $reticulate_p_value\n";
  File::Copy::move("$tempdir/$tempfile.ps", "$currentdir/$output_basename.ps");
  File::Copy::move("$tempdir/$tempfile.out", "$currentdir/$output_basename.out");
}

##############
# run sawyer #
##############

if ($run_silent) {
  $silent_switch = "/r";
} else {
  $silent_switch = "";
}
system "geneconv $tempfile $silent_switch > /dev/null";
open T, "$tempfile.sum";
@sum = <T>;
close T;
open T, "$tempfile.frags";
@frags = <T>;
close T;

# parse sum file for p values
foreach $i (0..$#sum) {
  if ($sum[$i] =~ /^\s*(.*) I and (.*) O significant global fragments/) {
    $global_inner = $1;
    $global_outer = $2;
    if ($global_inner =~ /no/i) { $global_inner = 0; }
    if ($global_outer =~ /no/i) { $global_outer = 0; }
  }
  if ($sum[$i] =~ /^\s*frags\s+Score\s+P-value/) {
    if ($i+1 <= $#sum && $sum[$i+1] =~ /^\s*SCORE\s+/) {
      @f = split /\s+/, $sum[$i+1];
      $inner_p = $f[2];
      if ($inner_p =~ /found/i) {
        undef $inner_p;
      }
    }
  }
  if ($sum[$i] =~ /^\s*OuterSeq/) {
    if ($i+1 <= $#sum && $sum[$i+1] =~ /^\s*frags\s+/) {
      @f = split /\s+/, $sum[$i+1];
      $outer_p = $f[3];	# because of leading space
      if ($outer_p =~ /found/i) {
        undef $outer_p;
      }
    }
  }
}

print "GENECONV Inner P-value: $inner_p\n";
print "GENECONV Outer P-value: $outer_p\n";
File::Copy::move("$tempdir/$tempfile.sum", "$currentdir/$output_basename.sum");
File::Copy::move("$tempdir/$tempfile.frags", "$currentdir/$output_basename.frags");

chdir $currentdir;

exit if ($inner_p > 0.05 && $outer_p > 0.05);

# parse silent output frags to get breaks
# we need @frags output and we need to read in the sequences
open F, $ARGV[0];
@s = <F>;
close F;
$sequence = slchen::fasta2hash(@s);
foreach $i (keys %$sequence) {
  $sequence_length = length($sequence->{$i});
  last;
}
@break = ();
foreach $out (@frags) {
  next if $out =~ /^#/;
  next if $out !~ /^G/;
  chomp $out;
  $out =~ s/>\s+//g;	# eliminate "> 1.0" which messes up our fields
  @f = split /\s+/, $out;
  push @break, $f[4];
  push @break, ($f[5] + 1);	# gives begin and end
				# end is end of fragment
				# but we want to use beginning of next fragment
}
push @break, 1;
push @break, $sequence_length+1;
@break = slchen::sortu @break;
foreach $i (1..$#break) {
  $alignment = "";
  $fragment_outputfile = "$output_basename.frag.$i\n";
  
  foreach $name (keys %$sequence) {
    $alignment .= ">$name\n";
    $alignment .= substr ($sequence->{$name}, $break[$i-1] - 1, $break[$i] - $break[$i-1]) . "\n";
    print STDERR "error, length ", $break[$i] - $break[$i-1], " from $break[$i-1] to $break[$i] not divisible by 3\n", join (" ", @break), "\n" if ($break[$i] - $break[$i-1]) % 3 != 0;
  }

  open F, ">$fragment_outputfile";
  print F $alignment;
  close F;
}
print join ("\t", @break), "\n";
