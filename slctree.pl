#!/usr/bin/perl -w
#
# take standard tree file
# optionally take information from another file to change how labels are drawn
# output will be SVG
#
use warnings;
use strict;
use Slctree;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my ($treefile, $optionfile, $outfile);
my (%treedata, %optiondata);
my @tree;
my $scalebar;
my %probability;
%optiondata = (
  'xsize' => 0,
  'ysize' => 0,
  'xoffset' => 100,
  'yoffset' => 100,
  'xscale' => 0,
  'yscale' => 15,
  'linewidth' => 2,
  'treespacing' => 5,	# space between trees, will get multiplied by yscale
  'numtrees' => 0,
  'numleaves' => 0,
  'branchlengthmax' => 0,
  'branchlengthmean' => 0,
  'bootstrapoffset' => 00,
  'bootstrapcutoff' => 50
);
$treefile = "";
$optionfile = "";
$outfile = 'out.svg';
$scalebar = 0;
GetOptions (
  'treefile=s' => \$treefile,
  'outfile=s' => \$outfile,
  'option=s' => \$optionfile,
  'xsize=i' => \$optiondata{xsize},
  'ysize=i' => \$optiondata{ysize},
  'xscale=f' => \$optiondata{xscale},
  'yscale=f' => \$optiondata{yscale},
  'scalebar!' => \$scalebar,
  'bootstrapoffset=f' => \$optiondata{bootstrapoffset},
  'bootstrapcutoff=f' => \$optiondata{bootstrapcutoff}
);
if (!length $treefile && $ARGV[0]) {
  $treefile = shift;
}

@tree = ();
if (!$treefile || $treefile eq '-') {
  @tree = <>;
} else {
  open T, $treefile;
  @tree = <T>;
  close T;
}
if ($outfile !~ /\.svg$/) { $outfile .= '.svg'; }

# parse options, get filenames, etc

# parse optional file to change labels
if (-f $optionfile) {
  &Slctree::parse_option($optionfile, \%optiondata);
}

# parse tree
#&Slctree::parse_tree($treefile, \%treedata, \%probability, \%optiondata);
&Slctree::parse_tree(\@tree, \%treedata, \%probability, \%optiondata);

# draw tree
&Slctree::draw_tree(\%treedata, \%optiondata, \%probability, $scalebar, $outfile);
