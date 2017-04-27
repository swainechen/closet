#!/usr/bin/perl -w
#
# use gri to draw a multiple line graphs
# use some column to tell which data series it goes to
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
use slchen;
$outfile = "";
$nolegend = 0;
$log = 0;	# log scale y axis
$xcol = 0;
$ycol = 1;
$seriescol = 2;
$colorcol = 3;
$xsize = 10;
$ysize = 10;
GetOptions ('outfile=s' => \$outfile,
            'nolegend' => \$nolegend,
            'xcol=i' => \$xcol,
            'ycol=i' => \$ycol,
            'seriescol=s' => \$seriescol,
            'colorcol=i' => \$colorcol,
            'xsize=f' => \$xsize,
            'ysize=f' => \$ysize,
            'log' => \$log);
if ($outfile !~ /\.ps$/) { $outfile .= ".ps"; }

@seriescol = parse_list($seriescol);
%series = ();
my $seriesname;

while (<>) {
  next if /^#/;
  chomp;
  @f = split /\s+/, $_;
  $seriesname = make_key(@f[@seriescol]);
  $series{$seriesname} = 1;
  push @{"x$seriesname"}, $f[$xcol];
  push @{"y$seriesname"}, $f[$ycol];
  if (defined $f[$colorcol] && isfloat($f[$colorcol])) {
    $color{$seriesname} = $f[$colorcol];
  }
}

open GRI, "| gri -nowarn_offpage -batch -output $outfile";
print GRI "set x size $xsize\n";
print GRI "set y size $ysize\n";
if ($log) {
  print GRI "set y type log\n";
}
print GRI "read columns x y\n";
foreach $seriesname (keys %series) {
  foreach $i (0..$#{"x$seriesname"}) {
    print GRI ${"x$seriesname"}[$i], "\t", ${"y$seriesname"}[$i], "\n";
  }
}
print GRI "\n";
print GRI "draw axes\n";
print GRI "delete columns\n";

@label = ();
$label = "";
$j = 0;
foreach $seriesname (keys %series) {
  print GRI "read columns x y\n";
  foreach $i (0..$#{"x$seriesname"}) {
    print GRI ${"x$seriesname"}[$i], "\t", ${"y$seriesname"}[$i], "\n";
  }
  print GRI "\n";
  print GRI "set color hsb ", $j/(scalar(keys %series)), " 1 1\n";
  if (defined $color{$seriesname}) {
    print GRI "set color hsb $color{$seriesname} 1 1\n";
  }
  print GRI "draw curve\n";
  $j++;
  if (!$nolegend) {
    print GRI ".labx. = { rpn ..xleft.. 0.98 ..xright.. ..xleft.. - * + }\n";
    if ($log) {
      print GRI ".laby. = { rpn ..ytop.. ..ytop.. ..ybottom.. / 0.05 power $j power / }\n";
    } else {
      print GRI ".laby. = { rpn ..ytop.. 0.05 $j * ..ytop.. ..ybottom.. - * - }\n";
    }
    @label = break_key($seriesname);
    $label = join ",", @label;
    print GRI "draw label \"Series $label\" rightjustified at .labx. .laby.\n";
  }
  print GRI "set color hsb 0 0 0\n";
  print GRI "delete columns\n";
}
close GRI;
