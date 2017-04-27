#!/usr/bin/perl -w
#
@a = <>;
@ind = ();
@nocomment = (0..$#a);
@comment = ();
foreach $i (reverse 0 .. $#a) {
  if ($a[$i] =~ /^#/) {
    splice @nocomment, $i, 1;
    push @comment, $i;
  }
}
foreach (@nocomment) {
  push @ind, rand();
}
@out = @nocomment[sort { $ind[$a] <=> $ind[$b] } (0..$#nocomment)];
foreach $i (reverse @comment) {
  splice @out, $i, 0, $i;
}
print @a[@out];
