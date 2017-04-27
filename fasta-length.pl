#!/usr/bin/perl -w
#
# print out length of fasta sequences
#
$length = 0;
$name = "";
while (<>) {
  if (/^>/) {
    if (length $name) {
      print $name;
      print $length, "\n";
      $name = "";
      $length = 0;
    }
    $name = $_;
  } else {
    chomp;
    $length += length $_;
  }
}
if (length $name) {
  print $name;
  print $length, "\n";
}
