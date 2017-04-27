#!/usr/bin/perl -w
#
# convert fasta format to phylip interleaved format
# phylip needs 10 character name, can be padded with spaces to longer length
# so make sure by adding 10 spaces to the name
# also truncate the fasta header to 10 characters if it's longer
# the fasta file must be aligned, probably using "-" for gaps, according to
# phylip
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$long = 0;
$sequential = 0;
&GetOptions (
  'sequential' => \$sequential,
  'long' => \$long
);
while (<>) {
  chomp;
  if (s/^>//) {
    $name = $_;
    push @names, $name;
    $seq{$name} = "";
    next;
  }
  s/\s//g;
  $seq{$name} .= $_;
}
$length = 0;
foreach $name (@names) {
  if (!$length) { $length = length $seq{$name}; }
  if (length $seq{$name} != $length) {
    print "First sequence length is $length but $name is ", length $seq{$name}, "\n";
    exit;
  }
}
$position = 0;
printf "%10d%10d\n", scalar @names, $length;
if ($sequential) {
  foreach $name (@names) {
    if ($long) {
      $outname = $name;
    } else {
      if (length $name > 10) {
        $outname = substr ($name, 0, 10);
      } else {
        $outname = $name;
        $outname .= " " until length $outname == 10;
      }
    }
    print $outname, "          ";
    print $seq{$name}, "\n";
  }
} else {
  while ($position <= $length) {
    if (!$position) {
      foreach $name (@names) {
        if ($long) {
          $outname = $name;
        } else {
          if (length $name > 10) {
            $outname = substr ($name, 0, 10);
          } else {
            $outname = $name;
            $outname .= " " until length $outname == 10;
          }
        }
        print $outname, "          ";
        print substr ($seq{$name}, 0, 60), "\n";
      }
      $position += 60;
    } else {
      print "\n";
      foreach $name (@names) {
        print substr ($seq{$name}, $position, 80), "\n";
      }
      $position += 80;
    }
  }
}
