#!/usr/bin/perl -w
#
$dir = `ls -1d jack*`;
chomp $dir;
@big5 = `iconv -f big5 -t utf-8 $dir/jack.freedb`;
@gbk = `iconv -f gbk -t utf-8 $dir/jack.freedb`;
print (grep /^DTITLE/, @big5);
print "BIG5 ok? (Y/n) ";
$ok = <>;
if ($ok !~ /^y/i && $ok !~ /^\s*$/) {
  print (grep m/^DTITLE/, @gbk);
  print "GBK ok? (Y/n) ";
  $ok = <>;
  if ($ok !~ /^y/i && $ok !~ /^\s*$/) {
    print "Neither BIG5 nor GBK encoding ok.  Aborting...\n";
    exit (-1);
  } else { @out = @gbk; }
} else { @out = @big5; }
open OUT, ">$dir/jack.freedb";
print OUT @out;
close OUT;
