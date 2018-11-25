#!/usr/bin/perl -w
#
use slchen;
use strict;
$|=1;

my $w = 1000;
my $s = 500;
open F, $ARGV[0];
my @a = <F>;
close F;
my $a = slchen::fasta2hash(@a);
my $length = 0;
my $i;
my $j;
my $real;
my @temp;
my @out;
my $temp;

foreach $i (keys %$a) {
  $length = length ($a->{$i});
  last;
}
#print join ("\t", "# Start pos", "Num bp", "Tajima's D", "Fu & Li's D*", "Fu & Li's F*", "Fay & Wu's H", "Fu's Fs", "Strobeck's S"), "\n";
print join ("\t", "# Start pos", "Num bp", "Fu's Fs", "Strobeck's S"), "\n";
for ($i = 0; $i < $length; $i += $s) {
  foreach $j (keys %$a) {
    $temp->{$j} = substr($a->{$j}, $i, $w);
  }
  @temp = slchen::hash2fasta($temp);
  @temp = slchen::align_stripgaps(@temp);
  $temp = slchen::fasta2hash(@temp);
  foreach $j (keys %$temp) {
    $real = length($temp->{$j});
    last;
  }
#  @out = ($i, $real, slchen::tajimaD($temp), slchen::fuliDF($temp), slchen::faywuH($temp), slchen::fufs($temp));
  @out = ($i, $real, slchen::fufs($temp));
  print STDOUT join("\t", @out), "\n";
}
