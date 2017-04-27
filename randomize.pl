#!/usr/bin/perl -w
#
# randomize input lines
# put all comments in the front
#
my @in = <>;
my @index = ();
foreach my $i (0..$#in) {
  if ($in[$i] =~ /^#/) { $index[$i] = -1/($i+1); }
  else { $index[$i] = rand(); }
}
@in = @in[sort {$index[$a]<=>$index[$b]} (0..$#in)];
print @in;
