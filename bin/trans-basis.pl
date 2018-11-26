#!/usr/bin/perl -w
#
# svd gives A = USV(t)
# so given another matrix A' will calculate U' with the same S and V(t)
# i.e. U = AV(t)(-1)S(-1)
# so U' = A'V(t)(-1)S(-1)
#
# note V(-1) = V(t) so this is actually U' = A'VS(-1)
#
use PDL;
use PDL::Slatec;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

$svdfile = '';
$Afile = '';
$bad = 0;
$all = 0;
$outfile = '';
$diagonal = 0;

GetOptions ('svd=s' => \$svdfile,
            'vectors=s' => \$Afile,
            'outfile=s' => \$outfile,
            'diagonal' => \$diagonal,
            'all' => \$all,
            'bad' => \$bad);

$PDL::IO::Misc::colsep = "\t";
@Vtrans = rcols ($svdfile, { INCLUDE => '/^Basis/' }, 1 .. 61);

$averageline = '';
$SVline = '';
open SVD, $svdfile;
while ($averageline !~ /^average/) { $averageline = <SVD>; }
while ($SVline !~ /^Singular/) { $SVline = <SVD>; }
close SVD;
$avg = pdl ((split /\t/, $averageline)[1..61]);
chomp $averageline;
chomp $SVline;
@SV = (split /\t/, $SVline)[1..61];

@W = dog transpose cat @Vtrans;
foreach $i (0 .. $#W) {
  if ($SV[$i] > 0.1) { $W[$i] /= $SV[$i]; }
  else { $W[$i] *= 0; }
}
$V = transpose cat @W;

#@datacols = (4..9,11..21,23..25,27..67);
@datacols = (3..8,10..20,22..24,26..66);
if ($diagonal) { @datacols = (1..61); }

open OUT, ">$outfile";
open IN, $Afile;

while (<IN>) {
  next if m/(^#)|(UNIQID)|(EWEIGHT)/;
  if (!$all) {
    next if ($bad && $_ !~ m/~NULL~/);
    next if (!$bad && m/~NULL~/);
  }
  $Aprime = pdl ((split /\t/, $_)[@datacols]);
  $Aprime -= $avg;
  print OUT join ("\t", list ($Aprime x $V)), "\n";
}

close IN;
close OUT;
