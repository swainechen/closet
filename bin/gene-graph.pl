#!/usr/bin/perl -w
#
# make a simple one-line graphic to show direction and location of genes
#

if ((!defined ($ARGV[0])) || ($ARGV[0] eq '-h') || ($ARGV[0] eq '--help')) { &PrintUsage; }

use Orgmap qw(:DEFAULT read_sequence $sequence $pttfile $genefield get_desc $orgname $subname);
&read_orgmap;

use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$outfile = '';
$noannotation = 0;
$xsize = 20;
$ysize = 20;
GetOptions ('outfile=s' => \$outfile,
            'noannotation' => \$noannotation,
            'x=i' => \$xsize,
            'y=i' => \$ysize);
if ($outfile !~ /\.ps$/) { $outfile .= ".ps"; }

# GRI location parameters
$left_border = 1;
$right_border = $xsize - $left_border;
$vertical_increment = 3;
$label_increment = 1.5;
$annotation_increment = 2;
$textline_increment = 0.5;
$gene_height = 0.2;
$negative_correction = 0.25;
$annotationy = $ysize + $vertical_increment;

$position = -1;
$gene = '';
$length = -1;

$pos_hue = -1;
$neg_hue = -1;
sub pos_hue {	# use colors close to red
  ++$pos_hue;
  if ($pos_hue >= 3) { $pos_hue -= 3; }
  my $return = $pos_hue/10 + 9/10;
  if ($return > 1) { --$return; }
  return ($return);
}
sub neg_hue {	# use colors close to blue
  ++$neg_hue;
  if ($neg_hue >= 3) { $neg_hue -= 3; }
  return (2/3 - 1/10 + $neg_hue/10);
}

open PTT, $pttfile;
@ptt = <PTT>;
close PTT;
&read_sequence;
$genome_length = length ($sequence);
if ($length > $genome_length) { print "$length greater than genome length ($genome_length).  Exiting.\n"; exit (1); }

open GRI, "| gri -nowarn_offpage -batch -output $outfile";
print GRI "set x size $xsize\n";
print GRI "set y size $ysize\n";
print GRI "set x axis 0 $xsize\n";
print GRI "set y axis 0 $ysize\n";
print GRI "set axes style none\n";
while (defined ($in = <>)) {
  next if $in =~ m/^#/;
  $position = -1; $gene = ''; $length = -1;
  chomp $in;
  ($input, $length) = split /\s+/, $in;
  if ((!defined $input) || (!defined $length)) { next; }
  if ($input =~ m/(\d+)\.\.(\d+)/) {
    @highlight = ($1, $2);
    $position = int (($1 + $2)/2);
  } elsif ($input =~ m/^[+-]?\d+$/) {
    $position = abs $input;
    @highlight = ($input);
  } else {
    $gene = $input;
    @highlight = ();
  }
  if ($position == -1) {
    foreach $line (@ptt) {
      if ($line =~ m/^\s*?(\d+)\.\.(\d+).*$gene/) { 
        ($s, $e) = ($1, $2);
        # a last sanity check
        @line = split /\t/, $line;
        if ($line[$genefield] eq $gene) {
          $position = int (($s + $e)/2);
        }
      }
    }
  }
  if ($position > -1) {
    $x = $left_border;
    $y = $annotationy - $vertical_increment;
    $labely = $y + $label_increment;
    $annotationy = $y - $annotation_increment;
    if ($noannotation) { $annotationy = $y - $textline_increment + $label_increment; }
    print GRI "set color hsb 0 0 0\n";
    print GRI "draw line from $left_border $y to $right_border $y\n";
    if (!$noannotation) {
      if ($gene ne '') {
        print GRI "draw label \"$orgname genome, $length bp centered on $gene\" at $x $labely\n";
      } else {
        print GRI "draw label \"$orgname genome, $length bp centered at $position\" at $x $labely\n";
      }
    }
    $draw_start = int($position - $length/2);
    $draw_end = int($position + $length/2);
    if ($draw_start < 0) { $wrap = 1; }			# left side neg
    elsif ($draw_start > $genome_length) { $wrap = 2; }	# right side big
    else { $wrap = 0; }
    foreach $line (@ptt) {
      if ($line =~ m/^\s*?(\d+)\.\.(\d+)/) { 
        $s = $1; $e = $2;
        if (!$wrap) {
          next if ($e < $draw_start);
          next if ($s > $draw_end);
        } elsif ($wrap == 1) {
          next if ($s > $draw_end && $e - $genome_length < $draw_start);
          if ($s > $draw_end) { $s -= $genome_length; $e -= $genome_length; }
        } elsif ($wrap == 2) {
          next if ($s + $genome_length > $draw_end && $e < $draw_start);
          if ($e > $draw_start) { $e += $genome_length; $s += $genome_length; }
        }
        chomp $line;
        @line = split /\t/, $line;
        if ($s < $draw_start) { $s = $draw_start; }
        if ($e > $draw_end) { $e = $draw_end; }
        $leftx = ($s - $draw_start)/$length * ($right_border - $left_border) + $left_border;
        $rightx = ($e - $draw_start)/$length * ($right_border - $left_border) + $left_border;
        $centerx = ($leftx + $rightx) / 2;
        @line = split /\t/, $line;
        $annotation = get_desc ($line);
        if ($line[1] =~ m/\+/) {
          $geney = $y + $gene_height;
          $labely = $geney + $gene_height;
          $hue = pos_hue();
          print GRI "set color hsb $hue 1 1\n";
          if ($gene eq $line[$genefield]) {
            print GRI "set color hsb 0 0 0\n";
          }
          print GRI "draw box filled $leftx $y $rightx $geney\n";
          print GRI "draw label \"$line[$genefield]\" centered at $centerx $labely\n";
          if (!$noannotation) {
            print GRI "draw label \"$line[$genefield] $annotation\" at $x $annotationy\n";
            $annotationy -= $textline_increment;
          }
        } else {
          $geney = $y - $gene_height;
          $labely = $geney - $gene_height - $negative_correction;
          $hue = neg_hue();
          print GRI "set color hsb $hue 1 1\n";
          if ($gene eq $line[$genefield]) {
            print GRI "set color hsb 0 0 0\n";
          }
          print GRI "draw box filled $leftx $geney $rightx $y\n";
          print GRI "draw label \"$line[$genefield]\" centered at $centerx $labely\n";
          if (!$noannotation) {
            print GRI "draw label \"$line[$genefield] $annotation\" at $x $annotationy\n";
            $annotationy -= $textline_increment;
          }
        }
      }
      if ($gene eq '') {
        print GRI "set color hsb 0 0 0\n";
        $hy1 = $y - 2*$gene_height;
        $hy2 = $y + 2*$gene_height;
        ($lastx, $lasty1, $lasty2) = (-1) x 3;
        foreach $highlight (@highlight) {
          $hx = ($highlight - $draw_start)/$length * ($right_border - $left_border) + $left_border;
          print GRI "draw line from $hx $hy1 to $hx $hy2\n";
          if ($lastx != -1 && $lasty1 != -1 && $lasty2 != -1) {
            print GRI "draw line from $hx $hy1 to $lastx $lasty1\n";
            print GRI "draw line from $hx $hy2 to $lastx $lasty2\n";
          }
          $lastx = $hx;
          $lasty1 = $hy1;
          $lasty2 = $hy2;
          $hy1 = $y;
          $hy2 = $y;
#          $hy1 += $gene_height;
#          $hy2 -= $gene_height;
        }
#        $centerx = ($right_border + $left_border)/2;
#        $bottomy = $y - 2*$gene_height;
#        $topy = $y + 2*$gene_height;
#        print GRI "set color hsb 0 0 0\n";
#        print GRI "draw line from $centerx $bottomy to $centerx $topy\n";
      }
    }
  }
}
close GRI;


sub PrintUsage {
  print "Usage: gene-organization.pl ORGCODE -o <outputfile.ps> [ -x <width in cm> ] [ -y <height in cm> ] [-noannotation ]\n";
  print "Then it takes two kinds of input on STDIN:\n";
  print "1:  [ position | gene ] length\n";
  print "2:  <start>..<end> length\n";
  print "For the first kind of input, you should only specify one of position or gene.  You must specify length.\n";
  print "For the second kind, either start or end can be smaller, and it will show the direction.  You must specify length.\n";
  print "The -noannotation flag causes only the line of genes to be drawn with their gene labels, no gene annotations below or location description above.\n";
  print "Output will be a postscript file.  Genes transcribed to the right will be above the line, and in orange, red, or magenta.  Genes transcried to the left will be below the line and in cyan, blue, or purple.\n";
  print "If you have specified a gene, it will be drawn in black above or below the line, accordingly.\n";
  exit (-1);
}
