#!/usr/bin/perl -w
#
# take standard tree file
# optionally take information from another file to change how labels are drawn
# output will be SVG
#
package Slctree;
require Exporter;
use warnings;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $nodeprefix %special_key);

@ISA = qw(Exporter);
@EXPORT = qw();
#@EXPORT = qw(parse_tree fill_tree next_level draw_tree draw_subtree svg_text svg_line set_leaves svg_stuff numerically sortfn sortme parse_option);
@EXPORT_OK = qw(parse_tree fill_tree next_level draw_tree draw_subtree svg_text svg_line set_leaves svg_stuff numerically sortfn sortme parse_option);
$nodeprefix = "Nodenumber";
%special_key = (
  'branchlength' => 1,
  'sortvalue' => 1,
  'bootstrap' => 1
);

sub parse_option {
  my ($f, $oref) = @_;
  # expect tab delimited file
  # first field is name/label
  # second is color
  my @f;
  $oref->{'color'} = {};
  open F, $f;
  while (<F>) {
    chomp;
    next if /^$/;
    next if /^#/;
    @f = split /\t/, $_;
    $oref->{'color'}->{$f[0]} = $f[1];
  }
  close F;
}

sub full_tree_text {
  # this calls tree_text
  # expects a top level tree data hash where first level keys are tree number
  # adds on the outside parentheses and probabilities and semicolon
  my ($dref, $pref) = @_;
  my $key;
  my $tree;
  my @tree;
  foreach $key (keys %$dref) {
    $tree = "(" . tree_text($dref->{$key}) . ")";
    if (defined $pref->{$key} && $pref->{$key} > 0) {
      $tree .= "[" . $pref->{$key} . "]";
    }
    $tree .= ";";
    push @tree, $tree;
  }
  return @tree;
}

sub tree_text {
  # print out tree text from tree data hash
  my ($dref) = @_;
  my $return;
  my $key;
  my $node;
  my @node;
  # each level should have at least special keys - branchlength, sortvalue, bootstrap
  foreach $key (keys %$dref) {
    next if $Slctree::special_key{$key};
    if ($key =~ /$nodeprefix/) {
      $node = "(" . tree_text($dref->{$key}) . ")";
      if (defined $dref->{$key}->{bootstrap} && $dref->{$key}->{bootstrap} > 0) {
        $node .= $dref->{$key}->{bootstrap};
      }
      if (defined $dref->{$key}->{branchlength} && $dref->{$key}->{branchlength} > 0) {
        $node .= ":" . $dref->{$key}->{branchlength};
      }
    } else {
      # we are at a terminal leaf
      $node = $key;
      # only look for branchlength, since we can't have a bootstrap for a terminal node
      if (defined $dref->{$key}->{branchlength} && $dref->{$key}->{branchlength} > 0) {
        $node .= ":" . $dref->{$key}->{branchlength};
      }
    }
    push @node, $node;
  }
  return join (",", @node);
}

sub all_children {
  # take tree data ref
  # return names of all children terminal leaves
  my ($dref) = @_;
  my $key;
  my @children = ();
  foreach $key (keys %$dref) {
    next if $Slctree::special_key{$key};
    if ($key =~ /$nodeprefix/) {
      push @children, all_children($dref->{$key});
    } else {
      push @children, $key;
    }
  }
  return @children;
}

sub parse_tree {
  my ($tref, $dref, $pref, $oref) = @_;
  # reference to tree array ($tref), usually one line per element, we will combine them here if needed
  #   tree data hash ($dref), first keys are tree number (for multiple trees)
  #   probability hash ($pref)
  #   option hash ($oref)
  # option hash to store number of trees and number of leaves for autoscaling
  my @tree = ();	# this is one tree string per element
  my $i = 0;	# lines in @$tref array
  my $j = 0;	# number of whole trees in the array
  foreach $i (0..$#{$tref}) {
    # combine trees if needed across lines (elements)
    chomp $tref->[$i];
    next if ($tref->[$i] =~ /^$/);
    next if ($tref->[$i] =~ /^#/);
    $tree[$j] .= $tref->[$i];
    if ($tref->[$i] =~ /(?:\[\d+\.?\d*\])?;$/) {
      $j++;
    }
  }
  $oref->{numtrees} = $j;
  foreach $j (0..$#tree) {
    $tree[$j] =~ s/\s//g;
    my $index = 0;	# we use this to number internal nodes
    $pref->{$j} = fill_tree($tree[$j], \%{$dref->{$j}}, \$index, $oref);
  }
  ($oref->{branchlengthmax}, $oref->{branchlengthmean}) = branchstats($dref);
}

sub branchstats {
  my ($dref) = @_;
  # get max and average branch length from root to leaves over all trees
  my ($max, $mean) = (0,0);
  my ($treemax, $treemean, $treetotal, $treeleaves);
  my $j;
  my $key;
  foreach $j (keys %$dref) {
    # the top level has no branchlengths, but we need some hash keys
    # so that the recursion works
    foreach $key (keys %Slctree::special_key) {
      $dref->{$j}->{$key} = 0;
#      $dref->{$j}->{branchlength} = 0;
#      $dref->{$j}->{sortvalue} = 0;
#      $dref->{$j}->{bootstrap} = 0;
    }
    ($treemax, $treetotal, $treeleaves) = treebranchstats(\%{$dref->{$j}});
    $mean += $treetotal/$treeleaves;
    $max = $treemax if $treemax > $max;
  }
  $mean /= scalar (keys %$dref);
  return ($max, $mean);
}

sub treebranchstats {
  my ($dref) = @_;
  # get max and average branch length from root to leaves for a single tree
  # always return absolute value of branch lengths because we might have
  # negative branch length if there is no branch length in the tree file
  my ($max, $total, $leaves) = 0;
  my ($submax, $subtotal, $subleaves);
  my $tempmax = 0;
  my $key;
  $total = 0;
  if (scalar (keys %$dref) > scalar keys %Slctree::special_key) {
  # we have a subtree
    foreach $key (keys %$dref) {
      next if $Slctree::special_key{$key};
#      next if $key eq 'sortvalue';
#      next if $key eq 'branchlength';
#      next if $key eq 'bootstrap';
      ($submax, $subtotal, $subleaves) = treebranchstats(\%{$dref->{$key}});
      $tempmax = $submax if $submax > $tempmax;
      if (defined $dref->{branchlength}) {
        $total += $subleaves * $dref->{branchlength} + $subtotal;
      } else {
        $total += $subtotal;
      }
      $leaves += $subleaves;
    }
    if (defined $dref->{branchlength}) {
      $max += $tempmax + $dref->{branchlength};
    } else {
      $max += $tempmax;
    }
  } else {	# we have a leaf, this is easy
    if (defined $dref->{branchlength}) {
      return (abs $dref->{branchlength}, abs $dref->{branchlength}, 1);
    } else {
      return (0, 0, 1);
    }
  }
  return ($max, $total, $leaves);
}

sub fill_tree {
  my ($tree, $dref, $iref, $oref) = @_;
  # tree in one line, reference to tree data hash, index for numbering internal
  # nodes, option data to keep track of number of terminal leaves
  # this returns the probability of the tree
  my @decompose = ();
  my @branch = ();
  my ($nodename, $trimmed, $branchlength, $probability);
  my @sort = ();
  my @f;
  undef $probability;
  if ($tree =~ /;$/) {	# this is the first iteration, we have the whole tree
    if ($tree =~ /\[(.*?)\];$/) {
      $probability = $1;
    } else {
      $probability = 1;
    }
    $tree =~ /^\((.*)\).*?;/;	# get out the most stuff inside parentheses
    $tree = $1;
  }
  # expect no surrounding parentheses
  # pull out subtrees
  @decompose = next_level($tree);

  # now we should be able to split on commas
  $trimmed = shift @decompose;
  @branch = split /,/, $trimmed;
  foreach my $i (0..$#branch) {
    if ($branch[$i] =~ /^\(\)/) {	# there is a subtree
      $nodename = "$nodeprefix".${$iref};
      ${$iref}++;
      # we expect after a close parentheses we may have two values
      # ()52:0.32 -- the 52 is a bootstrap value, 0.32 is branch length
      if ($branch[$i] =~ /\(\)(.*?:.*?)$/) {	# use branch length in node name
        @f = split /:/, $1;
        if (length $f[0]) {
          $dref->{$nodename}->{'bootstrap'} = $f[0];
        } else {
          $dref->{$nodename}->{'bootstrap'} = -1;
        }
        if (length $f[1]) {
          $dref->{$nodename}->{'branchlength'} = $f[1];
        } else {
          $dref->{$nodename}->{'branchlength'} = -1;
        }
      } else {
        $dref->{$nodename}->{'bootstrap'} = -1;
        $dref->{$nodename}->{'branchlength'} = -1;
      }
      fill_tree(shift (@decompose), \%{$dref->{$nodename}}, $iref, $oref);
      @sort = ();
      foreach my $key (keys %{$dref->{$nodename}}) {
        next if $key eq 'sortvalue';
        next if $key eq 'branchlength';
        next if $key eq 'bootstrap';
        push @sort, $dref->{$nodename}->{$key}->{'sortvalue'};
      }
      $dref->{$nodename}->{'sortvalue'} = sortfn($dref->{$nodename});
#      $dref->{$nodename}->{'sortvalue'} = sortme(@sort);
#      $dref->{$nodename}->{'sortvalue'} = $dref->{$nodename}->{'branchlength'};
    } else {			# no subtree - this shouldn't have bootstraps
      undef $branchlength;
      ($nodename, $branchlength) = split /:/, $branch[$i];
      $dref->{$nodename}->{'bootstrap'} = -1;
      if (defined $branchlength) {	# use branch length if available
        $dref->{$nodename}->{'branchlength'} = $branchlength;
      } else {
        $dref->{$nodename}->{'branchlength'} = -1;
      }
#      $dref->{$nodename}->{'sortvalue'} = sortme($nodename);
      $dref->{$nodename}->{'sortvalue'} = sortfn($dref->{$nodename});
      if ($branchlength != -1) {
        $dref->{$nodename}->{'sortvalue'} = $branchlength;
      } else {
        $dref->{$nodename}->{'sortvalue'} = 0;
      }
      $oref->{numleaves}++;	# increment number of terminal leaves
    }
  }
  if (defined $probability) { return ($probability); }
}

sub next_level {
  my ($tree) = @_;
  my $level = 0;
  my $i;
  my @tree = ();
  my @start = ();
  my @end = ();
  my @return = ();
  @tree = split //, $tree;
  foreach $i (0..$#tree) {
    if ($level < 0) {
      die "Got a negative level!  Position $i in $tree\n";
    }
    if ($tree[$i] eq '(') {
      $level++;
      if ($level == 1) {
        push @start, $i;
      }
    }
    if ($tree[$i] eq ')') {
      $level--;
      if ($level == 0) {
        push @end, $i;
      }
    }
  }
  if ($level != 0) {
    die "Got a level of $level after looking at whole tree line\n";
  }
  # return the tree with the next levels cut out, and the next levels
  foreach $i (reverse 0..$#start) {
    $return[$i] = substr ($tree, $start[$i]+1, $end[$i] - $start[$i] - 1);
    substr ($tree, $start[$i]+1, $end[$i] - $start[$i] - 1) = "";
  }
  return ($tree, @return);
}

sub draw_tree {
  my ($tref, $oref, $pref, $scalebar, $outfile) = @_;
  open O, ">$outfile";
  my ($svgheader, $svgfooter) = &svg_stuff($oref);
  my $pscale;	# scale size or color based on tree probability or something
  print O $svgheader;
  if ($scalebar) {
    $scalebar = sprintf("%.2f", ($oref->{branchlengthmax})/10);
    svg_line (0, 0, $scalebar, 0, "black", $oref);
    svg_text (0, 1, $scalebar, "black", $oref);
    $oref->{'yoffset'} += 4 * $oref->{'yscale'};
  }
  # the tref hash should start with integer keys for subtrees, at least has 0
  my $current = 0;
  my %leafpos;	# keep track of terminal leaf positions (vertical)
  foreach my $key (sort { $a <=> $b } keys %$tref) {
    %leafpos = ();
    $current = set_leaves(\%{$tref->{$key}}, \%leafpos, $current) + $oref->{treespacing};
    $pscale = $pref->{$key};
    &draw_subtree(\%{$tref->{$key}}, $oref, $pscale, \%leafpos, 0);
  }
  print O $svgfooter;
}

sub draw_subtree {
  my ($tref, $oref, $pscale, $lref, $x) = @_;
  # each node will have a sortvalue key
  # each node should have a branchlength key also
  # each node should also have a bootstrap key also
  # terminal leaves should have a position in $lref
  my @key = ();
  my $key;
  my ($max, $min, $color);
  undef $max;
  undef $min;
  foreach $key (keys %$tref) {
    next if $Slctree::special_key{$key};
#    next if $key eq 'branchlength';
#    next if $key eq 'sortvalue';
#    next if $key eq 'bootstrap';
    (!defined $max) && ($max = $lref->{$key});
    (!defined $min) && ($min = $lref->{$key});
    ($lref->{$key} > $max) && ($max = $lref->{$key});
    ($lref->{$key} < $min) && ($min = $lref->{$key});
    $color = "black";
    if (defined $oref->{'color'}->{$key}) {
      $color = $oref->{'color'}->{$key};
    }
    svg_line($x, $lref->{$key}, $x + $tref->{$key}->{'branchlength'}, $lref->{$key}, $color, $oref);
if ($key =~ /Urine/) { $color = "red"; }
if ($key =~ /Peri/) { $color = "magenta"; }
    if ($key !~ /$nodeprefix/) {
      svg_text($x + $tref->{$key}->{'branchlength'}, $lref->{$key}, $key, $color, $oref);
    }
    if (keys %{$tref->{$key}} > scalar keys %Slctree::special_key) {
    # all should by default have sortvalue and branchlength and bootstrap
      draw_subtree(\%{$tref->{$key}}, $oref, $pscale, $lref, $x + $tref->{$key}->{'branchlength'});
    }
  }
  if (defined $tref->{'bootstrap'} && $tref->{'bootstrap'} > 0 &&
      $tref->{'bootstrap'} >= $oref->{'bootstrapcutoff'}) {
    svg_text($x + $oref->{'bootstrapoffset'}/$oref->{'xscale'}, ($max + $min)/2, $tref->{'bootstrap'}, "black", $oref);
  }
  svg_line($x, $max, $x, $min, "black", $oref);
}

sub svg_line {
  # draw a line in svg using the output handle O
  my ($x1, $y1, $x2, $y2, $color, $oref) = @_;
  my $xscale = $oref->{'xscale'};
  my $yscale = $oref->{'yscale'};
  my $xoffset = $oref->{'xoffset'};
  my $yoffset = $oref->{'yoffset'};
  my $linewidth = $oref->{'linewidth'};
  $x1 = $x1 * $xscale + $xoffset;
  $x2 = $x2 * $xscale + $xoffset;
  $y1 = $y1 * $yscale + $yoffset;
  $y2 = $y2 * $yscale + $yoffset;
  $color = svg_color_spec($color);
#  if ($color ne "black") {
#    $color = 0 if $color < 0;
#    $color = 1 if $color > 1;
#    $red = int ($color * 255);
#    $green = int ((1-$color) * 255);
#    $color = "rgb($red,$green,255)";
#  }
  print O "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" style=\"stroke-width:$linewidth;stroke:$color;\" opacity=\"1\"/>\n";
}

sub svg_text {
  # draw text in svg
  my ($x, $y, $text, $color, $oref) = @_;
  my $xscale = $oref->{'xscale'};
  my $yscale = $oref->{'yscale'};
  my $xoffset = $oref->{'xoffset'};
  my $yoffset = $oref->{'yoffset'};
  $x = $x * $xscale + $xoffset + $yscale/2;
#  $y = $y * $yscale + $yoffset + $yscale/2;
  $y = $y * $yscale + $yoffset;
  $color = svg_color_spec($color);
  print O "<text x=\"$x\" y=\"$y\" fill=\"$color\" style=\"alignment-baseline:central\">\n";
  print O $text;
  print O "\n</text>\n";
}

sub svg_color_spec {
  # convert random stuff into a valid color spec
  my ($x) = @_;
  my %color_keynames = (
    aliceblue  => 1,
    antiquewhite => 1,
    aqua => 1,
    aquamarine => 1,
    azure => 1,
    beige => 1,
    bisque => 1,
    black => 1,
    blanchedalmond => 1,
    blue => 1,
    blueviolet => 1,
    brown => 1,
    burlywood => 1,
    cadetblue => 1,
    chartreuse => 1,
    chocolate => 1,
    coral => 1,
    cornflowerblue => 1,
    cornsilk => 1,
    crimson => 1,
    cyan => 1,
    darkblue => 1,
    darkcyan => 1,
    darkgoldenrod => 1,
    darkgray => 1,
    darkgreen => 1,
    darkgrey => 1,
    darkkhaki => 1,
    darkmagenta => 1,
    darkolivegreen => 1,
    darkorange => 1,
    darkorchid => 1,
    darkred => 1,
    darksalmon => 1,
    darkseagreen => 1,
    darkslateblue => 1,
    darkslategray => 1,
    darkslategrey => 1,
    darkturquoise => 1,
    darkviolet => 1,
    deeppink => 1,
    deepskyblue => 1,
    dimgray => 1,
    dimgrey => 1,
    dodgerblue => 1,
    firebrick => 1,
    floralwhite => 1,
    forestgreen => 1,
    fuchsia => 1,
    gainsboro => 1,
    ghostwhite => 1,
    gold => 1,
    goldenrod => 1,
    gray => 1,
    grey => 1,
    green => 1,
    greenyellow => 1,
    honeydew => 1,
    hotpink => 1,
    indianred => 1,
    indigo => 1,
    ivory => 1,
    khaki => 1,
    lavender => 1,
    lavenderblush => 1,
    lawngreen => 1,
    lemonchiffon => 1,
    lightblue => 1,
    lightcoral => 1,
    lightcyan => 1,
    lightgoldenrodyellow => 1,
    lightgray => 1,
    lightgreen => 1,
    lightgrey => 1,
    lightpink => 1,
    lightsalmon => 1,
    lightseagreen => 1,
    lightskyblue => 1,
    lightslategray => 1,
    lightslategrey => 1,
    lightsteelblue => 1,
    lightyellow => 1,
    lime => 1,
    limegreen => 1,
    linen => 1,
    magenta => 1,
    maroon => 1,
    mediumaquamarine => 1,
    mediumblue => 1,
    mediumorchid => 1,
    mediumpurple => 1,
    mediumseagreen => 1,
    mediumslateblue => 1,
    mediumspringgreen => 1,
    mediumturquoise => 1,
    mediumvioletred => 1,
    midnightblue => 1,
    mintcream => 1,
    mistyrose => 1,
    moccasin => 1,
    navajowhite => 1,
    navy => 1,
    oldlace => 1,
    olive => 1,
    olivedrab => 1,
    orange => 1,
    orangered => 1,
    orchid => 1,
    palegoldenrod => 1,
    palegreen => 1,
    paleturquoise => 1,
    palevioletred => 1,
    papayawhip => 1,
    peachpuff => 1,
    peru => 1,
    pink => 1,
    plum => 1,
    powderblue => 1,
    purple => 1,
    red => 1,
    rosybrown => 1,
    royalblue => 1,
    saddlebrown => 1,
    salmon => 1,
    sandybrown => 1,
    seagreen => 1,
    seashell => 1,
    sienna => 1,
    silver => 1,
    skyblue => 1,
    slateblue => 1,
    slategray => 1,
    slategrey => 1,
    snow => 1,
    springgreen => 1,
    steelblue => 1,
    tan => 1,
    teal => 1,
    thistle => 1,
    tomato => 1,
    turquoise => 1,
    violet => 1,
    wheat => 1,
    white => 1,
    whitesmoke => 1,
    yellow => 1,
    yellowgreen => 1
  );
  # if it's already valid, just return
  if ($x =~ /#[0-9A-Fa-f]{3}/ ||
      $x =~ /#[0-9A-Fa-f]{6}/ ||
      $x =~ /rgb\(.*\)/ ||
      defined $color_keynames{$x}
     ) {
    return $x;
  }
  # check if float
  if ($x =~ /^([+-])?(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
    $x = 0 if $x < 0;
    $x = 1 if $x > 1;
    my ($r, $g, $b);
    if ($x < 1/3) {
      $r = int(255 - 255*$x*3);
      $g = int(255*$x*3);
      $b = 0;
    } elsif ($x >= 1/3 && $x < 2/3) {
      $r = 0;
      $g = int(255 - 255*($x-1/3)*3);
      $b = int(255*($x-1/3)*3);
    } else {
      $r = int(255*($x-2/3)*3);
      $g = 0;
      $b = int(255 - 255*($x-2/3)*3);
    }
    return "rgb($r,$g,$b)";
  }
  return "black";
}

sub set_leaves {
  my ($tref, $lref, $current) = @_;
  my @key = ();
  my $key;
  foreach $key (keys %$tref) {
    next if $Slctree::special_key{$key};
#    next if $key eq 'sortvalue';
#    next if $key eq 'branchlength';
#    next if $key eq 'bootstrap';
    push @key, $key;
  }
  foreach $key (sort { $tref->{$a}->{'sortvalue'} <=> $tref->{$b}->{'sortvalue'} } @key) {
    if (keys %{$tref->{$key}} > scalar keys %Slctree::special_key) {
    # each will already have sortvalue and branchlength and bootstrap
      $current = set_leaves(\%{$tref->{$key}}, $lref, $current);
      my $nodepos = 0;
      foreach my $key2 (keys %{$tref->{$key}}) {
        next if $Slctree::special_key{$key2};
#        next if $key2 eq 'sortvalue';
#        next if $key2 eq 'branchlength';
#        next if $key2 eq 'bootstrap';
        $nodepos += $lref->{$key2};
      }
      # subtract number of special keys - sortvalue and branchlength and bootstrap
      $nodepos /= scalar (keys %{$tref->{$key}}) - (scalar keys %Slctree::special_key);
      $lref->{$key} = $nodepos;
    } else {
      $lref->{$key} = $current;
      $current++;
    }
  }
  return $current;
}

sub sortfn {
  my ($ref) = @_;
  my $return = 0;
  my $key;

  # this one sorts by number of terminal leaves
  foreach $key (keys %$ref) {
    next if $Slctree::special_key{$key};
    $return += $ref->{$key}->{'sortvalue'};
  }
  $return = 1 if !$return;

  return $return;
}

sub numerically {
  $a <=> $b;
}

sub sortme {
  my (@x) = @_;
  my $return = 0;
#my %r = qw(O157 1 ecoE 2 sfl2 3 sfl3 4 ecol 5 uti8 6 ecoC 7);
#if ($#x == 0) { return $r{$x[0]}; }
  foreach my $i (0..$#x) {
    $return += $x[$i];
  }
  if (scalar @x) {
    $return /= scalar @x;
  }
  return $return;
}

sub svg_stuff {
  my ($oref) = @_;

  my $svgheader =
'<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN"    "http://www.w3.org/TR/2001/REC
-SVG-20010904/DTD/svg10.dtd" [
  <!ENTITY ns_flows "http://ns.adobe.com/Flows/1.0/">
  <!ENTITY ns_extend "http://ns.adobe.com/Extensibility/1.0/">
  <!ENTITY ns_ai "http://ns.adobe.com/AdobeIllustrator/10.0/">
  <!ENTITY ns_svg "http://www.w3.org/2000/svg">
  <!ENTITY ns_xlink "http://www.w3.org/1999/xlink">
]>
  <!--width="629.285" height="799.257" viewBox="0 0 629.285 799.257" overflow="visible" enable-background="new 0 0 629.285 799.257"-->
<svg xmlns:x="&ns_extend;" xmlns:i="&ns_ai;" xmlns="&ns_svg;" xmlns:xlink="&ns_xlink;" xmlns:a="http://ns.adobe.com/AdobeSVGViewerExtensions/3.1/" width="__WIDTH__" height="__HEIGHT__" viewBox="0 0 __WIDTH__ __HEIGHT__" overflow="visible" enable-background="new 0 0 __WIDTH__ __HEIGHT__ " xml:space="preserve" id="applicationBase" >
  <g transform="translate(0,0),scale(1,1)">
';

  my $svgfooter = 
'</g>
</svg>
';

  my $width = $oref->{'xsize'};
  my $height = $oref->{'ysize'};
  if (!$height) {	# try to do some autoscaling
    $height = $oref->{numtrees} * $oref->{treespacing} + $oref->{numleaves} + 1;
    # do a +1 here to make room for a potential scalebar
    $height *= $oref->{yscale};
    $height += 2 * $oref->{yoffset};
    $oref->{ysize} = $height;
  }
  if (!$width) {	# try to do some autoscaling, based on height
    # see if we'll also need to set xscale
    if (!$oref->{xscale}) {
      # we can use either use branchlengthmax and branchlengthaverage
      $width = $height - 2 * $oref->{yoffset} + 2 * $oref->{xoffset};
      $oref->{xsize} = $width;
      $oref->{xscale} = ($width - 2 * $oref->{xoffset}) / $oref->{branchlengthmax};
    } else {
      $width = $oref->{branchlengthmax} * $oref->{xscale} + 2 * $oref->{xoffset};
    }
  }
  $svgheader =~ s/__WIDTH__/$width/g;
  $svgheader =~ s/__HEIGHT__/$height/g;
  return ($svgheader, $svgfooter);
}
