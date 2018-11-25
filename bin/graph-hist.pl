#!/usr/bin/perl -w
#
# use gri to draw a graph of a histogram, using output from hist-graph.pl or
# dual-hist.pl
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
my $outfile = "";
my $nolegend = 0;
my $xlog = 0;	# log scale x axis
my $ylog = 0;	# log scale y axis
my $type = 'bar';	# line or bar
my $mark = "";
my $xname = "Bin";
my $yname = "Frequency";
my $textbin = 0;	# whether to treat bins as text instead of numbers
my $symbolsize = 0.1;
my $fontsize = 12;	# this is the default GRI font size
my $bar_width = 0.8;	# fraction of space - 1 makes bars touch, less than
			# 1 makes bars have space between them
GetOptions (
  'outfile=s' => \$outfile,
  'nolegend' => \$nolegend,
  'mark=s' => \$mark,
  'symbolsize=f' => \$symbolsize,
  'type=s' => \$type,
  'barwidth=f' => \$bar_width,
  'xlog' => \$xlog,
  'ylog' => \$ylog,
  'xname=s' => \$xname,
  'yname=s' => \$yname,
  'textbin' => \$textbin,
  'fontsize=f' => \$fontsize
);

if ($outfile !~ /\.ps$/) { $outfile .= ".ps"; }

my ($xleft, $xright, $ybottom, $ytop);
my $width;
my $ncols = 0;
my @x = ();
while (<>) {
  next if /^##/;
  chomp;
  if (s/^#\s+//) {
    @label = split /\t/, $_;
    next;
  }
  @f = split /\t/, $_;
  if (!$ncols) { $ncols = scalar @f; }
  if (scalar @f != $ncols) {
    print "Line $#x has ", scalar (@f), " columns but first line has $ncols\n";
    exit;
  }
  push @x, $f[0];
  foreach $i (1..$#f) {
    push @{"y$i"}, $f[$i];
  }
}
open GRI, "| gri -nowarn_offpage -batch -output $outfile";
if ($xlog) {
  print GRI "set x type log\n";
}
if ($ylog) {
  print GRI "set y type log\n";
}

# get GRI to do autoscaling
# if we are doing a bar graph, we need an extra space on each side
# assume x values are sorted and evenly spaced
print GRI "read columns x y\n";
foreach $i (0..$#x) {
  foreach $j (1..$ncols - 1) {
    if ($textbin) {
      print GRI $i+1, "\t", ${"y$j"}[$i], "\n";
    } else {
      print GRI "$x[$i]\t", ${"y$j"}[$i], "\n";
    }
  }
}
if (defined $x[0] && defined $x[1]) {
  if ($type eq "bar") {		# extend x-axis a little bit
    if ($textbin) {
      print GRI "0\t", ${"y1"}[0], "\n";
      print GRI $#x+2, "\t", ${"y1"}[0], "\n";
    } else {
      print GRI $x[0] - ($x[1] - $x[0]), "\t", ${"y1"}[0], "\n";
      print GRI $x[$#x] + ($x[1] - $x[0]), "\t", ${"y1"}[0], "\n";
    }
  }
}
print GRI "\n";
print GRI "set x name \"$xname\"\n";
print GRI "set y name \"$yname\"\n";
if ($textbin) {
  # set x labels if they are not numbers
  print GRI "set x axis labels ";
  foreach $i (0..$#x) {
    print GRI $i+1, " \"$x[$i]\" ";
  }
  print GRI "\n";
}
print GRI "set font size $fontsize\n";
print GRI "draw axes\n";
print GRI "delete columns\n";

# the real graphing
foreach $j (1..$ncols - 1) {
  print GRI "set color hsb ", ($j-1)/($ncols-1)*2/3, " 1 1\n";
  if ($type eq 'line') {
    print GRI "read columns x y\n";
    foreach $i (0..$#x) {
      print GRI "$x[$i]\t", ${"y$j"}[$i], "\n";
    }
    print GRI "\n";
    if ($mark ne "") {
      print GRI "set symbol size $symbolsize\n";
      print GRI "draw symbol $mark\n";
    } else {
      print GRI "draw curve\n";
    }
  } elsif ($type eq 'bar') {
    if (defined $x[0] && defined $x[1]) {
      if ($textbin) {
        $width = 1;
      } else {
        $width = $x[1] - $x[0];
      }
      $width *= $bar_width;
    } else {
      die "Less than 2 x-values\n";
    }
    foreach $i (0..$#x) {
      # first calculate the coordinates
      if ($textbin) {
        $xleft = $i+1 - $width/2 + ($j-1) * ($width/($ncols-1));
        $xright = $xleft + $width/($ncols-1);
      } else {
        $xleft = $x[$i] - $width/2 + ($j-1) * ($width/($ncols-1));
        $xright = $xleft + $width/($ncols-1);
      }
#      $ybottom = 0;	# sometimes axis may not be at 0
      $ytop = ${"y$j"}[$i];
#      print GRI "draw box filled $xleft $ybottom $xright $ytop\n";
      print GRI "draw box filled $xleft ..ybottom.. $xright $ytop\n";
    }
  }
  if (!$nolegend) {
    print GRI ".labx. = { rpn ..xleft.. 0.98 ..xright.. ..xleft.. - * + }\n";
    if ($ylog) {
      print GRI ".laby. = { rpn ..ytop.. ..ytop.. ..ybottom.. / 0.05 power $j power / }\n";
    } else {
      print GRI ".laby. = { rpn ..ytop.. 0.05 $j * ..ytop.. ..ybottom.. - * - }\n";
    }
    if (defined $label[$j]) {
      print GRI "draw label \"$label[$j]\" rightjustified at .labx. .laby.\n";
    } else {
      print GRI "draw label \"Column $j\" rightjustified at .labx. .laby.\n";
    }
  }
  print GRI "set color hsb 0 0 0\n";
  print GRI "delete columns\n";
}
close GRI;
