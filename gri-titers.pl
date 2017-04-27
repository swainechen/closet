#!/usr/bin/perl -w
#
# use gri to graph titers
# use symbols for each point, log scale on y axis
# pick x axis labels from top
#
use slchen;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
my $outfile = "";
my $xlog = 0;
my $ylog = 1;
my $mediancolor = "0 0 0";
my $symbolsize = 0.25;
my $logclip = 1;	# if log scale, and value is <= 0, set it to this number
my $drawmedian = 1;
my $symbol = "bullet";
my $yaxislabel = "CFU";
my $xaxislabel = "Strain";
my $logtransform = 0;
my $scale = 1;
my $fontsize = 10;
my $boxwhisker = 0;	# whether to draw dots or box and whisker plot
GetOptions (
  'outfile=s' => \$outfile,
  'median!' => \$drawmedian,
  'logclip=f' => \$logclip,
  'mediancolor=s' => \$mediancolor,
  'symbolsize=f' => \$symbolsize,
  'symbol=s' => \$symbol,
  'xlog!' => \$xlog,
  'ylog!' => \$ylog,
  'logtransform!' => \$logtransform,
  'scale=f' => \$scale,
  'font=i' => \$fontsize,
  'xname=s' => \$xaxislabel,
  'yname=s' => \$yaxislabel,
  'xline=f' => \$xline,
  'yline=f' => \$yline,
  'ymin=f' => \$ymin,
  'box' => \$boxwhisker
);

if (!@ARGV) {
  print "Usage: $0 [ <options> ] <data file> -o <output file>\n";
  print "Options:\n";
  print <<__END__;
  -median|nomedian	whether to draw line at median (default true)
  -mediancolor <r g b>	color of median line (default 0 0 0 (black))
  -logclip <float>	what value to clip values at if log scale (default 1)
  -symbolsize <float>	size of symbols (default 0.1 cm)
  -symbol <symbol>	symbol to use (default bullet) (try plus, box, circ)
  -xlog|noxlog		whether x axis is log scale (default no)
  -ylog|noylog		whether y axis is log scale (default yes)
  -logtransform|nologtransform	whether to log transform data (default no)
  -scale <float>	multiply all y-values by this number (default 1)
  -font <fontsize>	font size for axis labels (default 10)
  -xname <text>		label for x-axis (default "Strain")
  -yname <text>		label for y-axis (default "CFU")
  -xline <float>	draw a vertical line at this x value (default undefined)
  -yline <float>	draw a horizontal line at this y value (default undefined)
  -ymin <float>		make sure y-axis includes this value (default undefined)
  -box			whether to do box-and-whisker plot (default no)
__END__

  exit;
}

if ($outfile !~ /\.ps$/) { $outfile .= ".ps"; }

my $xcoord;
my $ycoord;
my $gritext;
my $ncols = 0;
my @x = ();
$columnoffset = 0;
my @label = ();
my $addcols = 1;
my $min_positive;	# for automatic minimum of clipped log values
while (<>) {
  next if /^##/;
  chomp;
  if (s/^#\s*//) {
    $columnoffset = scalar @label;
    $addcols = 1;
    push @label, split (/\t/, $_);
    next;
  }
  @f = split /\t/, $_;
  if ($addcols) {
    $ncols += scalar @f;
    $addcols = 0;
  }
  push @x, $f[0];
  foreach $i (0..$#f) {
    next if !isfloat($f[$i]);
    $f[$i] *= $scale;
    if ($logtransform) {
      if ($f[$i] > 0) {
        $f[$i] = log($f[$i])/log(10);
        $min_positive = $f[$i] if !defined $min_positive || $f[$i] < $min_positive;
      } else {
        if ($logclip == 0) {	# do auto-figuring of this later
          $f[$i] = -1;
        } else {
          $f[$i] = log($logclip * $scale)/log(10);
        }
      }
    }
    push @{"y".($i+$columnoffset)}, $f[$i];
  }
}
open GRI, "| gri -nowarn_offpage -batch -output $outfile";
if ($xlog) {
  print GRI "set x type log\n";
}
if ($ylog) {
  print GRI "set y type log\n";
} else {
  $yaxislabel = 'Log$_{10}$ CFU' if $yaxislabel eq 'CFU';
}
print GRI "set axes style 1\n";
#print GRI "set x axis 0 ", $ncols+1, " 1\n";
print GRI "set x axis 0 ", scalar @label + 1, " 1\n";
if (@label) {
  print GRI "set x axis labels ";
#  foreach $i (0..$ncols - 1) {
  foreach $i (0..scalar @label - 1) {
    $gritext = $label[$i];
    $gritext =~ s/delta\s*/\$\\Delta\$/g;
    print GRI $i+1, " \"$gritext\" ";
  }
  print GRI "\n";
}
# get autoscaling of y axis
print GRI "set y axis unknown\n";
print GRI "read columns x y\n";
#foreach $i (0..$ncols-1) {
foreach $i (0..scalar @label - 1) {
  foreach $j (0..$#{"y$i"}) {
    ${"y$i"}[$j] = $min_positive if $logclip == 0 && ${"y$i"}[$j] == -1;
    $ycoord = ${"y$i"}[$j];
    if ($ycoord <= 0 && $ylog) { $ycoord = $logclip * $scale; }
    print GRI $i+1, " ", $ycoord, "\n";
  }
  if (defined $ymin) {
    print GRI $i+1, " ", $ymin, "\n";
  }
}
print GRI "\n";
print GRI "set font size $fontsize\n";
print GRI "set x name \"$xaxislabel\"\n";
print GRI "set y name \"$yaxislabel\"\n";
print GRI "draw axes\n";

# draw extra lines if we need them
if (defined $xline) {
  print GRI "draw line from $xline ..ybottom.. to $xline ..ytop..\n";
}
if (defined $yline) {
  if ($logtransform) { $yline = log($yline)/log(10); }
  print GRI "draw line from ..xleft.. $yline to ..xright.. $yline\n";
}

# draw the data
if ($boxwhisker) {
  foreach $i (0..scalar @label - 1) {
    $xcoord = $i + 1;
    print GRI gri_boxwhisker(\@{"y$i"}, $xcoord, 0.25, 0, 0) if scalar @{"y$i"};
  }
} else {
  print GRI "set symbol size $symbolsize\n";
  #foreach $i (0..$ncols-1) {
  foreach $i (0..scalar @label - 1) {
    $xcoord = $i+1;
    foreach $j (0..$#{"y$i"}) {
      $ycoord = ${"y$i"}[$j];
      if ($ylog && $ycoord <= 0) {
        $ycoord = $logclip;
      }
      print GRI "draw symbol $symbol at $xcoord $ycoord\n";
    }
  }

  # draw the medians
  print GRI "set line width 1.5\n";
  if ($drawmedian) {
  #  foreach $i (0..$ncols-1) {
    foreach $i (0..scalar @label - 1) {
      $xcoord = $i+1;
      my ($xl, $xr) = ($xcoord - 0.15, $xcoord + 0.15);
      $ycoord = array_median(@{"y$i"});
      if (defined $ycoord) {
        if (length $mediancolor) {
          print GRI "set color hsb $mediancolor\n";
        }
        print GRI "draw line from $xl $ycoord to $xr $ycoord\n";
      }
    }
  }
}

close GRI;
