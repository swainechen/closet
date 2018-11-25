#!/usr/bin/perl -w
#
# fix-coverage.pl changes names if too long
# fix these in the tree
# need original vcf and the tree that comes from it
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $vcf = "";
GetOptions (
  'vcf=s' => \$vcf
);

if (!-f $vcf) {
  print "Usage: $0 -vcf <vcf file> <tree file>\n";
  print "Looks for '##Name_change=...' in vcf file and replaces the remapped\n";
  print "names in the tree. Output goes to stdout.\n";
  exit;
}

my $line;
my %undo;
my $g;
my @f;
my $old;
my $new;
my $field;
my $count;

$line = 0;
open V, $vcf;
while (<V>) {
  last if !/^#/;
  $line++;
  if (/^##Name_change=Field (\d+) has name longer .*?; (\S+) changed to (Remap_\d+)$/) {
    $field = $1;
    $old = $2;
    $new = $3;
    # sanity check
    $g = $new;
    $g =~ s/^Remap_//;
    if ($field != $g) {
      print STDERR "Failed to parse line $line, mismatch in field and Remap name\n";
      exit;
    }
    $undo{$new} = $old;
  }
}
close V;

$line = 0;
while (<>) {
  $line++;
  foreach $g (keys %undo) {
    $count = $_ =~ s/$g:/$undo{$g}:/g;
    if ($count > 1) {
      print STDERR "Made $line substitutions on line for $g, exiting\n";
      exit;
    }
  }
  print;
}
