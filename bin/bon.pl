#!/usr/bin/perl -w
#
#  bon.pl - make a heatmap with labels from a text file
#  bon stands for box of numbers
#  will search for a largest box of numbers in lower right hand corner of file
#  will give options to take row/column immediately adjacent to data as
#  labels or join all rows/columns not part of data
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#   If this program helps you in your research and you publish images which
#   you've used slcview.pl to create, please cite:
#     S.L. Chen, unpublished.  http://slcview.sourceforge.net
#
############################################################################

use Image::Magick;
use Slcview;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

# initialize some variables

my %options;
my $outfile = "";
my $force = 0;
my $delimiter = "\t";
my ($printhelp, $printcolors, $printfonts, $printgnu) = (0,0,0,0);
my $debug = 0;
my $zero = 0;	# to center around some other value

# make xsize and ysize 10 by default so we can draw labels
$options{xsize} = 10;
$options{ysize} = 10;

GetOptions (
  'outfile=s' => \$outfile,
  'f' => \$force,
  'poscolor=s' => \$options{poscolor},
  'negcolor=s' => \$options{negcolor},
  'absentcolor=s' => \$options{absentcolor},
  'xsize=f' => \$options{xsize},
  'ysize=f' => \$options{ysize},
  'width=i' => \$options{width},
  'height=i' => \$options{height},
  'legend=s' => \$options{legend},
  'legsize=i' => \$options{legsize},
  'legnumber=i' => \$options{legnumber},
  'rowlabels|genelabels=i' => \$options{genelabels},
  'collabels|arraylabels=i' => \$options{arraylabels},
  'font=s' => \$options{font},
  'spacing=i' => \$options{spacing},
  'bgcolor|backgroundcolor=s' => \$options{bgcolor},
  'delimiter=s' => \$delimiter,
  'help' => \$printhelp,
  'listcolors|printcolors' => \$printcolors,
  'listfonts|printfonts' => \$printfonts,
  'listgnu|printgnu' => \$printgnu,
  'scaling=s' => \$options{scaling_function},
  'debug' => \$debug,
  'zero=f' => \$zero
);

# do help commands first, which will short circuit everything else
if ($printhelp) { &printhelp; }
if ($printcolors) { &Slcview::printcolors; }
if ($printfonts) { &Slcview::printfonts; }
if ($printgnu) { &Slcview::printgnu; }

%options = set_options (%options);

if (!$force && -e $outfile) {
  print STDOUT "$outfile exists, overwrite? (y/N) ";
  $response = <STDIN>;
  if ($response !~ m/^y/i) { exit (1); }
}

# parse the data file

my @input = <>;
my @data = ();	# this will be an array of arrays
		# syntax close to Programming Perl, 3rd ed. p.270

# first put @input into @data
foreach my $i (0..$#input) {
  chomp $input[$i];
  $data[$i] = [ split /$delimiter/, $input[$i] ];
}

# do some sanity checking on the data - look for same number of columns,
# same number of rows
# don't count columns in header lines which start with # - these get -1 cols
my @numcols = ();
my @colstart = ();	# where good data starts in that row
foreach my $row (0..$#data) {
  if ($data[$row][0] =~ /^#/) {
    $numcols[$row] = -1;
    $colstart[$row] = -1;
    next;
  }
  $numcols[$row] = $#{$data[$row]};

  $colstart[$row] = 0;		# if all data are good, will fall through below
  foreach my $j (reverse 0 .. $numcols[$row]) {
    # match C float based on Programming Perl, 3rd ed., p.190
    # also accept blank, nan, #div/0 - change these all to blank (absent)
    if (!defined $data[$row][$j]	||
        $data[$row][$j] =~ /^\s*$/	||
        $data[$row][$j] =~ /^nan$/i	||
        $data[$row][$j] =~ /^.?div.?0?/i) {
      $data[$row][$j] = "";
    } elsif (defined $data[$row][$j] &&
        $data[$row][$j] !~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
      $colstart[$row] = $j+1;
      last;	# it's defined and not a number and not a bad value as above
    }
  }
}
if ($debug) {
  print "numcols: ", join (" ", @numcols), "\n";
  print "colstart: ", join (" ", @colstart), "\n";
}

my $maxcols = arraymax(@numcols);
foreach my $row (0..$#numcols) {
  next if $numcols[$row] < 0;
  if ($numcols[$row] < $maxcols) {
    print STDERR "Row ", $row+1, " only has ", $numcols[$row]+1, " columns,";
    print STDERR " but ", $maxcols+1, " were expected.  Padding at the end";
    print STDERR " with missing values.\n";
    foreach my $j ($numcols[$row]+1 .. $maxcols) {
      $data[$row][$j] = "";
    }
  }
}

# doing same thing for columns, but a little trickier
# also allow for column comments - basically if element in first row has a #
# at the beginning of the field, consider the whole column a comment
# this should be ok for most people who just put a # at the beginning of the
# first line, and want the first column to be labels
my @numrows = ();
my @rowstart = ();	# where good data starts in that col
foreach $col (0..$maxcols) {
  if ($data[0][$col] =~ /^#/) {
    $numrows[$col] = -1;
    $rowstart[$col] = -1;
    next;
  }
  my $i = 0;
  # The loop below has the difficulty that if there is an undefined value
  # in the middle of the array it won't look past it - so don't leave an
  # undefined there, use a blank field
  while (defined $data[$i][$col]) { $i++; }
  $numrows[$col] = $i-1;
  
  $rowstart[$col] = 0;
  foreach my $i (reverse 0 .. $numrows[$col]) {
    if (!defined $data[$i][$col]	||
        $data[$i][$col] =~ /^\s*$/	||
        $data[$i][$col] =~ /^nan$/i	||
        $data[$i][$col] =~ /^.?div.?0?/i) {
      $data[$i][$col] = "";	# this should have been done above, but...
    } elsif (defined $data[$i][$col] &&
        $data[$i][$col] !~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
      $rowstart[$col] = $i+1;
      last;
    }
  }
}
if ($debug) {
  print "numrows: ", join (" ", @numrows), "\n";
  print "rowstart: ", join (" ", @rowstart), "\n";
}

my $maxrows = arraymax(@numrows);
foreach my $col (0..$#numrows) {
  next if $numrows[$col] < 0;
  if ($numrows[$col] < $maxrows) {
    print STDERR "Column ", $col+1, " only has ", $numrows[$col]+1, " rows,";
    print STDERR " but ", $maxrows+1, " were expected.  Padding at the end";
    print STDERR " with missing values.\n";
    foreach my $i ($numcols[$col]+1 .. $maxrows) {
      $data[$i][$col] = "";
    }
  }
}

# now we have data between colstart and numcols, rowstart and numrows
# look for max rowstart among columns that are valid in the last row
# and max colstart among rows that are valid in the last column
# and that should give us our box of numbers
my $toprow = arraymax (@rowstart[$colstart[$maxrows] .. $maxcols]);
my $leftcol = arraymax (@colstart[$rowstart[$maxcols] .. $maxrows]);
if ($debug) { print "toprow: $toprow\nleftcol: $leftcol\n"; }
my @commentrows = ();		# these are only comments interspersed with
my @commentcols = ();		# data - there may be other comments also
my $countcols = 0;
my @box = ();
foreach my $i ($toprow .. $maxrows) {
  if ($data[$i][0] =~ /^#/) {
    push @commentrows, $i;
    next;
  }
  foreach my $j ($leftcol .. $maxcols) {
    if ($data[0][$j] =~ /^#/) {
      if (!$countcols) {	# only count these once
        push @commentcols, $j;
        $countcols = 1;
      }
      next;
    }
    if ($data[$i][$j] =~ /^\s*$/) {
      push @box, "";
    } else {
      push @box, $data[$i][$j] - $zero;
    }
  }
}
my $boxrows = $maxrows - $toprow + 1 - scalar @commentrows;
my $boxcols = $maxcols - $leftcol + 1 - scalar @commentcols;

# have to mangle xsize/ysize and width/height
if ($options{width} == -1) {
  $options{width} = sprintf "%i", $boxcols * $options{xsize};
} else {
  $options{xsize} = $options{width}/$boxcols;
}
if ($options{height} == -1) {
  $options{height} = sprintf "%i", $boxrows * $options{ysize};
} else {
  $options{ysize} = $options{height}/$boxrows;
}

my $heatmap = Image::Magick->new;
if ($debug) { print "row $boxrows\ncol $boxcols\ndata:\n", join (" ", @box), "\n"; }
$heatmap = draw_heatmap (\@box, $boxrows, $boxcols, %options);

# get gene/row labels together
@rowlabels = ();
foreach my $i ($toprow .. $maxrows) {
  my $comment = 0;
  foreach my $j (@commentrows) {
    $comment = 1 if $i == $j;
  }
  next if $comment;
  my @label = ();
  foreach my $j (0..$maxcols) {		# catch all the comment columns
    push @label, $data[$i][$j] if $data[0][$j] =~ /^#/;
  }
  push @rowlabels, join (" ", @label);
}
if ($debug) { print "rowlabels: ", join (" ", @rowlabels), "\n"; }
# autoscale if not specified
if ($options{genelabels} < 0) {
  $options{genelabels} = label_scale(\@rowlabels, $options{ysize});
}
my $rlabimage = Image::Magick->new;
$rlabimage = draw_labels (\@rowlabels, $options{genelabels}, $options{height}, %options);
$debug && print "rlab ", $rlabimage->Get('width', 'height'), "\n";

# get array/column labels together
@collabels = ();
foreach my $j ($leftcol .. $maxcols) {
  my $comment = 0;
  foreach my $i (@commentcols) {
    $comment = 1 if $j == $i;
  }
  next if $comment;
  my @label = ();
  foreach my $i (0..$maxrows) {
    push @label, $data[$i][$j] if $data[$i][0] =~ /^#/;
  }
  push @collabels, join (" ", @label);
}
if ($debug) { print "collabels: ", join (" ", @collabels), "\n"; }
my $clabimage = Image::Magick->new;
# autoscale if not specified
if ($options{arraylabels} < 0) {
  $options{arraylabels} = label_scale(\@collabels, $options{xsize});
}
$clabimage = draw_labels (\@collabels, $options{arraylabels}, $options{width}, %options);
$clabimage->Rotate(degrees=>270);
$debug && print "clab ", $clabimage->Get('width', 'height'), "\n";

$debug && print "alab ", $options{arraylabels}, " glab ", $options{genelabels}, "\n";

# put these together
# should look like this
# ---------------
# |    1    | 2 |	1 = clabimage
# |---------|---|	2 = urspacer
# |         |   |	3 = heatmap
# |    3    | 4 |	4 = rlabimage
# |         |   |
# ---------------
# put together 1 and 3 first, then 2 and 4, then 1-3 and 2-4

my $urheight = $options{arraylabels};
my $urwidth = $options{genelabels};
my $urspacer = Image::Magick->New;
$urspacer->Read("xc:$options{bgcolor}");
$urspacer->Scale($urwidth.'x'.$urheight.'!');
$debug && print "urspacer width: ", $urspacer->Get('width'), "\n";

my $imageleft = Image::Magick->New;
$imageleft = stack_TB ($options{spacing}, $options{bgcolor}, $clabimage, $heatmap);
$debug && print "imageleft width: ", $imageleft->Get('width'), "\n";
my $imageright = Image::Magick->New;
$imageright = stack_TB ($options{spacing}, $options{bgcolor}, $urspacer, $rlabimage);
$debug && print "imageright width: ", $imageright->Get('width'), "\n";

my $outputimage = Image::Magick->New;
$outputimage = stack_LR ($options{spacing}, $options{bgcolor}, $imageleft, $imageright);

$outputimage->Write(filename=>$outfile, dither=>'False');

sub arraymax {
  my @a = @_;
  my $max = $a[0];
  foreach my $i (0..$#a) {
    if ($a[$i] > $max) {
      $max = $a[$i];
    }
  }
  return ($max);
}

sub label_scale {
  my ($labels, $boxsize) = @_;
  my @len = ();
  my $mean = 0;
  my $stdev = 0;

  # do auto-detection of $xres for draw_labels if needed - try to make it
  # long enough for mean length of label + 1 standard deviation
  foreach my $lab (@$labels) {
    push @len, length $lab;
  }
  # calculate mean and stdev
  foreach my $len (@len) {
    $mean += $len;
  }
  $mean /= scalar @len;
  foreach my $len (@len) {
    $stdev += ($len - $mean)*($len - $mean);
  }
  $stdev /= scalar @len;
  $stdev = sqrt ($stdev);
  return (int ($boxsize * $mean));
}

sub stack_LR {
  my ($spacing, $color, @pieces) = @_;
  my $index = 0;
  my $height = 0;
  my $piece;
  my $image = new Image::Magick;
  my $spacer;
  $image->Set(Adjoin=>'True');

  return (undef) if !scalar @pieces;	# no images given
  $piece = shift @pieces;
  while (!defined $piece || !defined ($piece->Get('width'))) {
    last if !scalar @pieces;
    $piece = shift @pieces;
  }

  if (defined $piece) {
    $height = $piece->Get('height');
    return ($height) if ($height <= 0);
    $image->[$index] = $piece;

    while ($piece = shift @pieces) {
      if (defined $piece && $piece->Get('width') > 0) {
        if ($piece->Get('height') == $height) {
          $spacer = new Image::Magick;
          $spacer->Read("xc:$color");
          $spacer->Scale($spacing.'x'.$height.'!');
          ++$index;
          $image->[$index] = $spacer;
          ++$index;
          $image->[$index] = $piece;
        } else {
          return (-1);		# not all images were the same height;
        }
      }
    }
    if ($index > 0) {
      return ($image->Append(stack=>'False'));
    } else {
      return ($image->[0]);
    }
  } else {
    return (undef);			# no good images found
  }
}

sub stack_TB {
  my ($spacing, $color, @pieces) = @_;
  my $index = 0;
  my $width = 0;
  my $piece;
  my $image = new Image::Magick;
  my $spacer;
  $image->Set(Adjoin=>'True');

  return (undef) if !scalar @pieces;	# no images given
  $piece = shift @pieces;
  while (!defined $piece || !defined ($piece->Get('height'))) {
    last if !scalar @pieces;
    $piece = shift @pieces;
  }

  if (defined $piece) {
    $width = $piece->Get('width');
    return ($width) if ($width <= 0);
    $image->[$index] = $piece;

    while ($piece = shift @pieces) {
      if (defined $piece && $piece->Get('height') > 0) {
        if ($piece->Get('width') == $width) {
          $spacer = new Image::Magick;
          $spacer->Read("xc:$color");
          $spacer->Scale($width.'x'.$spacing.'!');
          ++$index;
          $image->[$index] = $spacer;
          ++$index;
          $image->[$index] = $piece;
        } else {
$debug && print "error $width not equal ", $piece->Get('width'), "\n";
          return (-1);		# not all images were the same width;
        }
      }
    }
    if ($index > 0) {
      return ($image->Append(stack=>'True'));
    } else {
      return ($image->[0]);
    }
  } else {
    return (undef);			# no good images found
  }
}

sub printhelp {
  print "Usage: $0 <options>\n";
  exit;
}
