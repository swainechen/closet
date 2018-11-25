#!/usr/bin/perl -w
#
# make a circular map
# mark off beginning and end
# label in the middle
# mark specified locations with tick marks, in color
# positive marks to outside, negative to inside
# can specify color with number after location
# default colors are red on outside, cyan on inside (negative)
# if a line starts with "label" then the next field will be taken as
#   a text label, and it will be written below the genome name in the
#   center with the hue specified by the third field of that line
#
use Orgmap qw(:DEFAULT read_sequence $sequence $topology $orgname $subname $pttfile $genefield);
&read_orgmap;
&read_sequence;
$genome_length = length($sequence);
open PTT, $pttfile;
@ptt= <PTT>;
close PTT;
foreach $line (@ptt) {
  if ($line =~ m/^\s*(\d+)\.\.(\d+)\s/) {
    $center = ($1 + $2)/2;
    @f = split /\t/, $line;
    $f[1] =~ s/\s+//;	# strand
    $ptt{$f[$genefield]} = $f[1].$center;
  }
}
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$size = 10;
$outfile = '';
$tickscale = 1;
$huescale = 1;
$circlescale = 1;
$numcircles = 1;
$linewidth = 0.709;
GetOptions ('size=i' => \$size,
            'tickscale=f' => \$tickscale,
            'huescale=f' => \$huescale,
            'circlescale=f' => \$circlescale,
            'linewidth=f' => \$linewidth,
            'numcircles=i' => \$numcircles,
            'outfile=s' => \$outfile);
if ($outfile !~ m/\.ps$/) { $outfile .= '.ps'; }

$ticklength = 0.5;
$label_increment = 0.2;
$positive_hue = 0;
$negative_hue = 0.5;
$hue = -1;
$pi = 3.14159265359;
$radius = $size/2 - $ticklength - 2 * $ticklength * $numcircles;
$x_off = $radius + $ticklength + 2 * $ticklength * $numcircles;
$y_off = $radius + $ticklength + 2 * $ticklength * $numcircles;
$text_x = $x_off;
$text_y = $y_off;
$text_off = $radius/10;
%circle_drawn = ();

open GRI, "| gri -nowarn_offpage -batch -output $outfile";
if ($topology == 0) {
  print GRI "draw circle with radius $radius at $x_off $y_off\n";
  $circle_drawn{0} = 1;
  print GRI "draw label \"$orgname genome, $genome_length bp\" centered at $x_off $y_off cm\n";
  ($x1, $y1) = cartesian($radius + 1, 0);
  ($x2, $y2) = cartesian($radius - 1, 0);
  print GRI "draw line from $x1 $y1 to $x2 $y2 cm\n";
  $x1 += $label_increment;
  print GRI "draw label \"1\" at $x1 $y1 cm\n";
  $x1 -= 2*$label_increment;
  print GRI "draw label \"$genome_length\" rightjustified at $x1 $y1 cm\n";
  print GRI "set line width $linewidth\n";
  while (<>) {
    next if /^#/;
    chomp;
    s/^\s+//;
    s/\s+$//;
    @f = split;
    $cmd = shift @f;
    if ($cmd eq "label") {	# label <something> <hue>
      $text_hue = pop @f;
      if ($text_hue =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/) {	# decimal (perl book)
        while ($text_hue > 1) { $text_hue -= 1; }
      } else {
        push @f, $text_hue;
        $text_hue = -1;
      }
      $text_label = join " ", @f;
      $text_y -= $text_off;
      if ($text_hue < 0) {
        print GRI "set color hsb 0 0 0\n";
      } else {
        print GRI "set color hsb $text_hue 1 1\n";
      }
      print GRI "draw label \"$text_label\" centered at $text_x $text_y cm\n";
      next;
    } elsif ($cmd =~ /^(\d+)\.\.(\d+)$/) {	# 100..200 <hue|length|circle>
      $coordinate = ($1+$2)/2;
      if ($1 > $2) { $negative = 1; } else { $negative = 0; }
      ($hue, $length, $circle) = hlc($f[0]);
    } elsif ($cmd =~ /^([-+]?\d+)$/) {		# -500 <hue|length|circle>
      if ($1 < 0) { $negative = 1; } else { $negative = 0; }
      $coordinate = abs($1);
      ($hue, $length, $circle) = hlc($f[0]);
    } else {			# assume this is CC0400 <hue|length|circle>
      $cmd =~ s/^[-+]//;
      if (defined $ptt{$cmd}) {
        $coordinate = abs ($ptt{$cmd});
        if ($ptt{$cmd} < 0) { $negative = 1; } else { $negative = 0; }
        if (defined $f[0]) { ($hue, $length, $circle) = hlc($f[0]); }
        else { ($hue, $length, $circle) = (-1, 1, 0); }
      } else {
        next;
      }
    }

    if (!$circle_drawn{$circle}) {
      print GRI "set color hsb 0 0 0\n";
      print GRI "draw circle with radius ", $radius + 2 * $ticklength * $circle, " at $x_off $y_off\n";
      $circle_drawn{$circle} = 1;
    }
    ($x1, $y1) = cartesian($radius + $circle*2*$ticklength, $coordinate/$genome_length);
    if (!$negative) {
      ($x2, $y2) = cartesian($radius + $circle*2*$ticklength + $ticklength*$length, $coordinate/$genome_length);
      if ($hue == -1) { $hue = $positive_hue; }
    } else {
      ($x2, $y2) = cartesian($radius + $circle*2*$ticklength - $ticklength*$length, $coordinate/$genome_length);
      if ($hue == -1) { $hue = $negative_hue; }
    }
    print GRI "set color hsb $hue 1 1\n";
    print GRI "draw line from $x1 $y1 to $x2 $y2 cm\n";
    $hue = -1;
  }
}
close GRI;

sub cartesian {		# a little different because start 0 rad at (0,1)
			# go clockwise as radians increase
			# actually don't use radians, but relative distance
			# around the circle, i.e. 1 is all the way around again
  my ($r, $theta) = @_;
  my $x = $r * sin(2*$pi*$theta) + $x_off;
  my $y = $r * cos(2*$pi*$theta) + $y_off;
  return ($x, $y);
}

sub hlc {
  my ($s) = @_;
  $s = "" if !defined $s;
  my @s = split /\|/, $s;
  if (defined $s[0] && $s[0] =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/) {
    $s[0] *= $huescale;
    while ($s[0] > 1) { $s[0] -= 1; }
  } else {
    $s[0] = -1;
  }
  if (defined $s[1] && $s[1] =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/) {
    $s[1] *= $tickscale;
  } else {
    $s[1] = 1;
  }
  if (defined $s[2] && $s[2] =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/) {
    $s[2] *= $circlescale;
  } else {
    $s[2] = 0;
  }
  return ($s[0], $s[1], $s[2]);
}
