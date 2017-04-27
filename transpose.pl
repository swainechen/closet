#!/usr/bin/perl -w
#
# transpose rows and columns
#
$delim = "\t";
$comments = 0;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
GetOptions ('d=s' => \$delim,
            'comments' => \$comments);

$row = 0;
$col = 0;
while (<>) {
  if (!$comments) { next if /^[#>]/; }
  chomp;
  @{"f$row"} = split /$delim/, $_;
  if ($#{"f$row"} > $col) { $col = $#{"f$row"}; }
  $row++;
}

foreach $j (0..$col) {
  @out = ();
  foreach $i (0..$row) {
    if (defined ${"f$i"}[$j]) {
      push @out, ${"f$i"}[$j];
    } else {
      push @out, "";
    }
  }
  print join ($delim, @out), "\n";
}
