#!/usr/bin/perl -w #
# use gri to draw a graph of a histogram, using output from hist-graph.pl or
# dual-hist.pl
#
use slchen;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$outfile = "";
$legend = 1;
$symbolsize = 0.1;
$mark = "bullet";
$color_col = -1;
$scale_color = 0;	# if true, then take numbers in color_col and autoscale
			# to range 0 (red) to 2/3 (blue)
$huescale = 1;		# colors will be distributed from 0 to $huescale
$xlog = 0;
$ylog = 0;
$fit = 1;		# whether to draw a regression line
$background = 0;	# whether to draw a background color
$line = 0;		# whether to draw line between points, they need to
			# be ordered if so
$xaxislabel = "";
$yaxislabel = "";

GetOptions ('outfile=s' => \$outfile,
            'xlog' => \$xlog,
            'ylog' => \$ylog,
            'xname=s' => \$xaxislabel,
            'yname=s' => \$yaxislabel,
            'symbolsize=f' => \$symbolsize,
            'mark=s' => \$mark,
            'colorcol=i' => \$color_col,
            'scalecolor!' => \$scale_color,
            'huescale=f' => \$huescale,
            'fit!' => \$fit,
            'background' => \$background,
            'line' => \$line,
            'legend!' => \$legend);
if ($outfile !~ /\.ps$/) { $outfile .= ".ps"; }

$ncols = 0;
@x = ();
@label = ();
while (<>) {
#  next if /^#/;
  chomp;
  if (s/^#\s+?//) {
    @label = split /\t/, $_;
    next;
  }
  @f = split /\t/, $_;
  if (!$ncols) { $ncols = scalar @f; }
  if (scalar @f != $ncols) {
    print "Line $#x has ", scalar (@f), " columns but first row has $ncols\n";
    exit;
  }
  push @x, $f[0];
  foreach $i (1..$#f) {
    push @{"y$i"}, $f[$i];
  }
}

if ($color_col >= 0) {
  %color = ();
  &set_color_hash (\%color, \@{"y$color_col"});
}

open GRI, "| gri -nowarn_offpage -batch -output $outfile";
if ($background) {
  print GRI "set bounding box 0 0 15 15 cm\n";
  print GRI "set color hsb 1 0 1\n";
  print GRI "draw box filled 0 0 15 15 cm\n";
}
if ($xlog) {
  print GRI "set x type log\n";
}
if ($ylog) {
  print GRI "set y type log\n";
}
print GRI "read columns x y\n";
foreach $i (0..$#x) {
  foreach $j (1..$ncols - 1) {
    next if $j == $color_col;
    print GRI "$x[$i]\t", ${"y$j"}[$i], "\n";
  }
}
print GRI "\n";
#if (defined $xaxislabel && length($xaxislabel)) {
#  print GRI "set x name \"$xaxislabel\"\n";
#} elsif (defined $label[0]) {
#  print GRI "set x name \"$label[0]\"\n";
#}
#if (defined $yaxislabel && length($yaxislabel)) {
#  print GRI "set y name \"$yaxislabel\"\n";
#} elsif (defined $label[1]) {
#  print GRI "set y name \"$label[1]\"\n";
#}
print GRI "draw axes\n";
if ($line) {
  print GRI "draw curve\n";
} else {
  print GRI "delete columns\n";
  print GRI ".labx. = { rpn ..xleft.. 0.98 ..xright.. ..xleft.. - * + }\n";
  print GRI ".laby. = ..ytop..\n";
  print GRI "set symbol size $symbolsize\n";
  foreach $j (1..$ncols - 1) {
    next if $j == $color_col;
    if ($color_col >= 0) {
      foreach $i (0..$#x) {
        print GRI "set color hsb ", $huescale * $color{${"y$color_col"}[$i]}, " 1 1\n";
        print GRI "draw symbol $mark at $x[$i] ", ${"y$j"}[$i], "\n";
      }
      print GRI "read columns x y\n";
      foreach $i (0..$#x) {
        print GRI "$x[$i]\t", ${"y$j"}[$i], "\n";
      }
      print GRI "\n";
    } else {
      print GRI "read columns x y\n";
      foreach $i (0..$#x) {
        print GRI "$x[$i]\t", ${"y$j"}[$i], "\n";
      }
      print GRI "\n";
      print GRI "set color hsb ", ($j-1)/($ncols-1), " 1 1\n";
      print GRI "draw symbol $mark\n";
    }

    if ($legend) {
      if ($ylog) {
        print GRI ".laby. = { rpn .laby. ..ytop.. ..ybottom.. / 0.05 power / }\n";
      } else {
        print GRI ".laby. = { rpn .laby. 0.05 ..ytop.. ..ybottom.. - * - }\n";
      }
#      print GRI "set color hsb 0 0 0\n";
      if (defined ($label[$j]) && defined ($label[0])) {
        print GRI "draw label \"$label[$j] (y) vs $label[0] (x)\" rightjustified at .labx. .laby.\n";
      } else {
        print GRI "draw label \"col $j (y) vs col 0 (x)\" rightjustified at .labx. .laby.\n";
      }

      # labels for colors
      if ($color_col >= 0) {
        foreach $colorkey (sort keys %color) {
          if ($ylog) {
            print GRI ".laby. = { rpn .laby. ..ytop.. ..ybottom.. / 0.05 power / }\n";
          } else {
            print GRI ".laby. = { rpn .laby. 0.05 ..ytop.. ..ybottom.. - * - }\n";
          }
          print GRI "set color hsb $color{$colorkey} 1 1\n";
          print GRI "draw label \"$colorkey\" rightjustified at .labx. .laby.\n";
          print GRI "set color hsb 0 0 0\n";
        }
      }
    }

    # do the regression line
    if ($fit) {
      print GRI "regress y vs x\n";
      print GRI ".ly. = ..ybottom..\n";
      print GRI ".lx. = { rpn ..ybottom.. ..coeff0.. - ..coeff1.. / }\n";
      print GRI "if { rpn .lx. ..xleft.. > }\n";
      print GRI ".lx. = ..xleft..\n";
      print GRI ".ly. = { rpn ..coeff0.. ..coeff1.. ..xleft.. * + }\n";
      print GRI "else if { rpn .lx. ..xright.. < }\n";
      print GRI ".lx. = ..xright..\n";
      print GRI ".ly. = { rpn ..coeff0.. ..coeff1.. ..xright.. * + }\n";
      print GRI "end if\n";
      print GRI ".ry. = ..ytop..\n";
      print GRI ".rx. = { rpn ..ytop.. ..coeff0.. - ..coeff1.. / }\n";
      print GRI "if { rpn .rx. ..xright.. < }\n";
      print GRI ".rx. = ..xright..\n";
      print GRI ".ry. = { rpn ..coeff0.. ..coeff1.. ..xright.. * + }\n";
      print GRI "else if { rpn .rx. ..xleft.. > }\n";
      print GRI ".rx. = ..xleft..\n";
      print GRI ".ry. = { rpn ..coeff0.. ..coeff1.. ..xleft.. * + }\n";
      print GRI "end if\n";
      print GRI "\\r_value = \"..R2..\"\n";
      print GRI "set color hsb 0 0 0\n";
      print GRI "draw line from .lx. .ly. to .rx. .ry.\n";
    }

    if ($legend) {
      if ($ylog) {
        print GRI ".laby. = { rpn .laby. ..ytop.. ..ybottom.. / 0.05 power / }\n";
      } else {
        print GRI ".laby. = { rpn .laby. 0.05 ..ytop.. ..ybottom.. - * - }\n";
      }
      if ($fit) {
        print GRI "draw label \"R2 = \\\@r_value\" rightjustified at .labx. .laby.\n";
      }
    }
    print GRI "delete columns\n";
  }
}
close GRI;

sub set_color_hash {
  my ($c_ref, $d_ref) = (@_);
  # if they are all numbers between 0 and 1 then just use them
  # if not then count how many different labels and assign different colors
  # data should be in @$d_ref array
  # assign colors to %$c_ref hash
  my $already_ok = 1;
  my $all_float = 1;
  my $max = $d_ref->[0];
  my $min = $d_ref->[0];
  my $i;
  foreach $i (0..$#{$d_ref}) {
    if (!isfloat($d_ref->[$i])) {
      $all_float = 0;
      $already_ok = 0;
    } else {
      $max = $d_ref->[$i] if $d_ref->[$i] > $max;
      $min = $d_ref->[$i] if $d_ref->[$i] < $min;
      if ($d_ref->[$i] < 0 || $d_ref->[$i] > 1) {
        $already_ok = 0;
      }
    }
    $c_ref->{$d_ref->[$i]} = 0;
  }
  if ($already_ok) {
    foreach $i (0..$#{$d_ref}) {
      $c_ref->{$d_ref->[$i]} = $d_ref->[$i];
    }
  } else {
    if ($scale_color && $all_float) {
      my $scale = $max - $min;
      $scale = 1 if $scale == 0;
      foreach $i (0..$#{$d_ref}) {
        $c_ref->{$d_ref->[$i]} = ($d_ref->[$i] - $min) / $scale * 2/3;
      }
      print "Min (red):\t$min\n";
      print "Max (blue):\t$max\n";
    } else {
      my $total = scalar keys %$c_ref;
      $i = 0;
      foreach my $key (sort keys %$c_ref) {
        $c_ref->{$key} = $i/$total;
        $i++;
      }
    }
  }
}
