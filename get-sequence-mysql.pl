#!/usr/bin/perl -w
#
# take a target, strainid, grab the sequence
# provide some other information also
#
use strict;
use DBI;
use slchen;
use Orgmap;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $strainid;
my @strainid;
my $target;
my $sequence;
my $short = 0;
my $include = "";
my $fill = 0;	# if true then put N's in for include region and pull all sequences

GetOptions (
  'strainid=s' => \$strainid,
  'target=s' => \$target,
  'short' => \$short,
  'include=s' => \$include,
  'fill!' => \$fill
);

if (!defined $strainid || !length $strainid ||
    !defined $target || !length $target) {
  print "Usage: $0 -strainid <strain id list> -target <target> [ -short ] [ -include <start..stop> ]\n";
  exit;
}

@strainid = slchen::parse_list($strainid);

my $dbh = slchen::dbconnect("resequencing");
my $sql = "SELECT Sequence, Start2, End2, Start1, End1 FROM Resequences WHERE StrainID = ? AND Target = ?";
my $sth = $dbh->prepare($sql);
my @data;
my ($s, $e);
my $info_sql = "SELECT Name, Patient, IsolationLocation FROM straindb.HultgrenStrains WHERE straindb.HultgrenStrains.StrainID = ?";
my $info_sth = $dbh->prepare($info_sql);
my ($start, $stop) = (0, 0);
if (length $include) {
  ($start, $stop) = split /\.\./, $include;
}

my $switch;
my ($s1, $e1);
my $opposite = 0;
foreach $strainid (@strainid) {
  undef $sequence;
  undef $s;
  undef $e;
  $sth->execute($strainid, $target);
  while (@data = $sth->fetchrow_array) {
    ($sequence, $s, $e, $s1, $e1) = @data;
    if ($s > $e) {
      $sequence = Orgmap::revcomp($sequence);
      $switch = $s;
      $s = $e;
      $e = $switch;
      $opposite = 1;
    }
  }
  if (!$fill) {
    next if !defined $sequence || !length $sequence;
    next if ($start && $s > $start);
    next if ($stop && $e < $stop);
  } else {
    if (!defined $s) { $s = 0; }
    if (!defined $e) { $e = 0; }
  }
  print ">StrainID$strainid";
  if (!$short) {
    $info_sth->execute($strainid);
    while (@data = $info_sth->fetchrow_array) {
      print "|";
      print "$s..$e";
      print "|";
      print "$data[0]" if defined $data[0] && length $data[0];
      print "|";
      print "Patient$data[1]" if defined $data[1] && length $data[1];
      print "|";
      print "$data[2]" if defined $data[2] && length $data[2];
    }
  }
  print "\n";
  if ($fill) {
    my $outseq;
    my $outlength = $stop - $start + 1;
    if (!defined $sequence || !length $sequence) {
      print "N" x $outlength, "\n";
      next;
    }
    my ($low, $high) = ($s, $e);
    # correct this in case have nonaligned sequence
    if ($opposite) {
      $low -= length($sequence) - $e1;
      $high += $s1 - 1;
    } else {
      $low -= $s1 - 1;
      $high += length($sequence) - $e1;
    }
    if ($low > $stop || $high < $start) {
      # we didn't get anything where we wanted it
      $outseq = "";
    } else {
      # take care of the left side, trim or pad N's
      if ($low <= $start) {
        $outseq = substr($sequence, $start - $low);
      }
      if ($low > $start) {
        $outseq = "N" x ($low - $start);
        $outseq .= $sequence;
      }
    }
    # make total length correct
    if (length($outseq) < $outlength) {
      $outseq .= "N" x ($outlength - length($outseq));
    }
    print substr($outseq, 0, $outlength), "\n";
  } else {
    print "$sequence\n";
  }
}
