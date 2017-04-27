#!/usr/bin/perl -w
use PDL;
($x, $y) = rcols (*STDIN);
$w = ones (nelem $x);
print join "\t", r2($x, $y, $w);
print "\n";

sub r2 {
  my ($x, $y, $w) = @_;
  my ($mx, $sx) = (stats ($x, $w))[0,1];
  my ($my, $sy) = (stats ($y, $w))[0,1];
  $x = ($x - $mx) / $sx;
  $y = ($y - $my) / $sy;
  my $r = abs(sum($x * $y * $w)/sum ($w));
  my $z = 0.5*(log(1+$r) - log(1-$r));
  my $zpdl = pdl ($z*sqrt(sum($w)-3)/sqrt(2));
  my ($p) = list (erfc $zpdl);
  return ($r**2, $p);
}
