#!/usr/bin/perl -w
#
# give rank sums from two files
#
use slchen;
open F, $ARGV[0];
while (<F>) {
  chomp;
  next if /^#/ || /^$/;
  push @a, $_;
}
close F;
open F, $ARGV[1];
while (<F>) {
  chomp;
  next if /^#/ || /^$/;
  push @b, $_;
}

@c = (@a, @b);
@c = sort {$a <=> $b} @c;

$ranka = 0;
$rankb = 0;
my $counts;
my $ranks;
foreach $i (0..$#c) {
  push @{$counts->{$c[$i]}}, $i+1;
}
foreach $i (keys %$counts) {
  $ranks->{$i} = array_mean(@{$counts->{$i}});
}

foreach $i (0..$#a) {
  $ranka += $ranks->{$a[$i]};
}
foreach $i (0..$#b) {
  $rankb += $ranks->{$b[$i]};
}

print "File\tN\tSum of ranks:\n";
print "$ARGV[0]\t", scalar @a, "\t$ranka\n";
print "$ARGV[1]\t", scalar @b, "\t$rankb\n";
