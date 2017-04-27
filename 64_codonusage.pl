#!/usr/bin/perl -w
#
# general script to draw boxes for 64 values for each codon/nucleotide triplet
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
use vars qw($A $C $D $E $F $G $H $I $K $L $M $N $P $Q $R $S $T $V $W $X $Y);

$outfile = 'test.ps';
GetOptions ('outfile=s' => \$outfile);
open GRI, "| gri -batch -nowarn_offpage -output $outfile";
my %data = ();
my %count = (); 
while (<>) {
  next if /^[#>]/;
  chomp;
  my @f = split /\t/, $_;
  my $key = $f[0];
  my $aa = $f[1];
  my $number = $f[3];
  $usage{$key} = $number;
  if (aa($key) ne $aa) { die aa($key)." ne $aa, $key translation error!"; }
  if (defined $count{$aa}) { ++$count{$aa}; }
  else { $count{$aa} = 1; }
  if (defined ${$aa}) { ${$aa} += $number; }
  else { ${$aa} = $number; }
}
my $total = 0;
foreach my $aa (keys %count) {
  $total += $count{$aa};
}
if ($total != 64) { die "wrong number of data points"; }
foreach $codon (keys %usage) {
  if (${aa($codon)}) {
    $data{$codon} = $usage{$codon}/${aa($codon)} - 1/$count{aa($codon)};
  } else {
    $data{$codon} = 0 - 1/$count{aa($codon)};
  }
}
&aa_color;
&draw_boxes ("test", \%data, 2, 0);
close GRI;

sub draw_boxes {
  # take a name as first argument to label graph
  # next take a ref to a hash where keys are trinucleotides
  # finally take a scale and offset - subtract offset, divide by scale
  my ($name, $ref, $scale, $offset) = @_;
  my (@f, $key, $uckey, $centerx, $centery, $llx, $urx, $lly, $ury);
  my ($side, $hue, $sat, $bright);
  print GRI "delete scale\n";
  print GRI "set axes style default\n";
  print GRI "set color rgb 0 0 0\n";
  print GRI "set x axis 0.5 4.5 1\n";
  print GRI "set x name \"First position\"\n";
  print GRI "set y axis 0.5 4.5 1\n";
  print GRI "set x name \"Second position\"\n";
  print GRI "set font size 0\n";
  print GRI "draw axes\n";
  print GRI "set font size 12\n";
  print GRI "draw label \"G\" centered at 1 0.15\n";
  print GRI "draw label \"A\" centered at 2 0.15\n";
  print GRI "draw label \"T\" centered at 3 0.15\n";
  print GRI "draw label \"C\" centered at 4 0.15\n";
  print GRI "draw label \"First position\" centered at 2.5 -0.15\n";
  print GRI "draw label \"G\" centered at 0.2 0.93\n";
  print GRI "draw label \"A\" centered at 0.2 1.93\n";
  print GRI "draw label \"T\" centered at 0.2 2.93\n";
  print GRI "draw label \"C\" centered at 0.2 3.93\n";
  print GRI "draw label \"Second position\" centered at 0.05 2.5 rotated 90\n";

  print GRI "draw label \"$name\" at 0 5\n";
  print GRI "draw label \"scale: $scale\" at 0 4.8\n";
  print GRI "draw label \"offset: $offset\" at 0 4.6\n";

  foreach $key (keys %{$ref}) {
    $uckey = uc $key;
    if ($uckey =~ m/^[GATC]{3}$/) {
      @f = split //, $uckey;
    } else {
      @f = ("X", "X", "X");
    }
    if ($f[0] eq 'X') { $centerx = 1/3; }
    elsif ($f[0] eq 'G') { $centerx = 1; }
    elsif ($f[0] eq 'A') { $centerx = 2; }
    elsif ($f[0] eq 'T') { $centerx = 3; }
    elsif ($f[0] eq 'C') { $centerx = 4; }
    if ($f[1] eq 'X') { $centery = 1/3; }
    elsif ($f[1] eq 'G') { $centery = 1; }
    elsif ($f[1] eq 'A') { $centery = 2; }
    elsif ($f[1] eq 'T') { $centery = 3; }
    elsif ($f[1] eq 'C') { $centery = 4; }
    if ($f[2] eq 'G') { $centerx -= 0.2; $centery += 0.2; }
    elsif ($f[2] eq 'A') { $centerx += 0.2; $centery += 0.2; }
    elsif ($f[2] eq 'T') { $centerx -= 0.2; $centery -= 0.2; }
    elsif ($f[2] eq 'C') { $centerx += 0.2; $centery -= 0.2; }

    # center values around offset - hopefully will provide more dynamic range
    $side = (${$ref}{$key} - $offset)/$scale;
#    $hue = 0;
#    $sat = 0;
#    $bright = 0;
    ($hue, $sat, $bright) = (hsb($uckey));

    $llx = $centerx - $side;
    $lly = $centery - $side;
    $urx = $centerx + $side;
    $ury = $centery + $side;

    print GRI "set color hsb $hue $sat $bright\n";
    if ($side > 0) {
      print GRI "draw box filled $llx $lly $urx $ury\n";
    } else {
      print GRI "draw box $llx $lly $urx $ury\n";
    }
  }
}

sub hsb {
  my ($l) = @_;
  if (($l =~ m/[gatcGATC]{3}$/) && (aa($l) ne 'X')) {
    return (${aa($l)}, 1, 1);
  }
  return (0, 0, 0);
}

sub aa {
  my ($codon) = @_;
  use vars qw($TTT $TTC $TTA $TTG $TCT $TCC $TCA $TCG $TAT $TAC $TAA $TAG $TGT $TGC $TGA $TGG $CTT $CTC $CTA $CTG $CCT $CCC $CCA $CCG $CAT $CAC $CAA $CAG $CGT $CGC $CGA $CGG $ATT $ATC $ATA $ATG $ACT $ACC $ACA $ACG $AAT $AAC $AAA $AAG $AGT $AGC $AGA $AGG $GTT $GTC $GTA $GTG $GCT $GCC $GCA $GCG $GAT $GAC $GAA $GAG $GGT $GGC $GGA $GGG);
  $codon =~ tr/gatc/GATC/;
  $TTT="F";
  $TTC="F";
  $TTA="L";
  $TTG="L";
  $TCT="S";
  $TCC="S";
  $TCA="S";
  $TCG="S";
  $TAT="Y";
  $TAC="Y";
  $TAA="X";
  $TAG="X";
  $TGT="C";
  $TGC="C";
  $TGA="X";
  $TGG="W";
  $CTT="L";
  $CTC="L";
  $CTA="L";
  $CTG="L";
  $CCT="P";
  $CCC="P";
  $CCA="P";
  $CCG="P";
  $CAT="H";
  $CAC="H";
  $CAA="Q";
  $CAG="Q";
  $CGT="R";
  $CGC="R";
  $CGA="R";
  $CGG="R";
  $ATT="I";
  $ATC="I";
  $ATA="I";
  $ATG="M";
  $ACT="T";
  $ACC="T";
  $ACA="T";
  $ACG="T";
  $AAT="N";
  $AAC="N";
  $AAA="K";
  $AAG="K";
  $AGT="S";
  $AGC="S";
  $AGA="R";
  $AGG="R";
  $GTT="V";
  $GTC="V";
  $GTA="V";
  $GTG="V";
  $GCT="A";
  $GCC="A";
  $GCA="A";
  $GCG="A";
  $GAT="D";
  $GAC="D";
  $GAA="E";
  $GAG="E";
  $GGT="G";
  $GGC="G";
  $GGA="G";
  $GGG="G";
  return ($$codon);
}

sub aa_color {
  # draw amino acid legend - colored
  $G = 0;
  $A = 0.05;
  $V = 0.10;
  $L = 0.15;
  $I = 0.20;
  $F = 0.25;
  $W = 0.30;
  $M = 0.35;
  $C = 0.40;
  $S = 0.45;
  $T = 0.50;
  $Y = 0.55;
  $N = 0.60;
  $Q = 0.65;
  $D = 0.70;
  $E = 0.75;
  $R = 0.80;
  $K = 0.85;
  $H = 0.90;
  $P = 0.95;
  $X = 0;
}

sub draw_aa_labels {
  print GRI "set axes style none\n";
  print GRI "set symbol size 0.3\n";
  print GRI ".leg_x. = {rpn ..xleft.. ..xright.. ..xleft.. - 0.5 * +}\n";
  print GRI ".leg_y. = ..ybottom..\n";
  foreach $aa('G','A','V','L','I','F','W','M','C','S','T','Y','N','Q','D','E','R','K','H','P') {
    my $hue = ${$aa};
    print GRI "set color hsb $hue 1 1\n";
    print GRI "draw symbol legend bullet \"$aa\" at .leg_x. .leg_y.\n";
    print GRI ".leg_y. += {rpn ..ytop.. ..ybottom.. - 20 /}\n";
  }
}
