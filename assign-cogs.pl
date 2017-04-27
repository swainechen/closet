#!/usr/bin/perl -w
#
# take two orgcodes
# first orgcode will get assigned cogs
# second orgcode has cogs which will be used as reference
# use bestblast data to make assignments
# assign only if not assigned already if $replace == 0
#
$replace = 0;
$identity_min = 80;
$align_min = 0.8;

use Getopt::Long;
&Getopt::Long::Configure("pass_through");
GetOptions (
  'replace' => \$replace,
  'identity=f' => \$identity_min,
  'align=f' => \$align_min
);

use Orgmap qw(:DEFAULT $pttfile);
&read_orgmap;
$ptt1 = $pttfile;
$org1 = $orgcode;

undef $pttfile;
undef $orgcode;

&read_orgmap;
$ptt2 = $pttfile;
$org2 = $orgcode;

%blast = ();
%cog = ();
use DBI;
my $dbh = DBI->connect('DBI:Oracle:host="sh-pod00.wustl.edu";sid="test0"', 'slchen', 'swaine', { LongReadLen => 100000000 });
my $sth = $dbh->prepare('SELECT gid1, gid2, identity, align1, align2 FROM blast WHERE org1 in (?,?) and org2 in (?,?)');
$sth->execute($org1, $org2, $org1, $org2);
while (@data = $sth->fetchrow_array()) {
  if ($data[2] >= $identity_min && $data[3] >= $align_min && $data[4] >= $align_min) {
    $blast{$data[0]} = $data[1];
    $blast{$data[1]} = $data[0];
  }
}

my @f;
open PTT, $ptt2;
while (<PTT>) {
  if (/^\s*\d+\.\.\d+\s+/) {
    @f = split /\t/, $_;
    $f[6] && ($cog{$f[3]} = $f[6]);
  }
}
close PTT;

open PTT, $ptt1;
while (<PTT>) {
  if (/^\s*\d+\.\.\d+\s+/) {
    @f = split /\t/, $_;
    if ($f[6] =~ /-/ || $replace) {
      if ($blast{$f[3]} && $cog{$blast{$f[3]}}) {
        $f[6] = $cog{$blast{$f[3]}};
      }
    }
    print join ("\t", @f);
  } else {
    print;
  }
}
close PTT;
