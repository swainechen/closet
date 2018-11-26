#!/usr/bin/perl -w
#
# use gnuplot to do a scree plot
#
my @sv = ();
my $i;
my @f;
my $lines = 0;
my $outfile = 'scree.png';
my $cutoff = 0.00001;

while (<>) {
  chomp;
  @f = split /\t/, $_;
  foreach $i (0..$#f) {
    if ($f[$i] =~ /^([+-])?(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
      last if $f[$i] < $cutoff;
      push @sv, $f[$i];
    }
  }
  $lines++;
  last if $lines == 1 && $#f >= 1;
}

open OUT, "| gnuplot";
print OUT "set terminal png\n";
print OUT "set output \"$outfile\"\n";
print OUT "plot '-' with points\n";
foreach $i (0..$#sv) {
  print OUT $i+1, "\t$sv[$i]\n";
}
print OUT "e\n";
