#!/usr/bin/perl -w
#
# dump data from a sql query
#
use slchen;
use DBI;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$sid = 'germs_browser';
$user = "";
$pass = "";
$host = "";
$port = 3306;
$quiet = 0;
$single = 0;	# if true, will consider multiline input as single sql
GetOptions (
  'sid|db=s' => \$sid,
  'quiet' => \$quiet,
  'single' => \$single,
  'user=s' => \$user,
  'pass=s' => \$pass,
  'host=s' => \$host,
  'port=i' => \$port
);
my $dbh = slchen::dbconnect($sid, $host, $user, $pass, $port);
if (!@ARGV || -f $ARGV[0]) {
  while (<>) {
    chomp;
    next if /^#/ || /^$/;;
    push @sql, $_;
  }
} else {
  @sql = @ARGV;
}

if ($single) {
  @sql = (join " ", @sql);
}

foreach $sql (@sql) {
  chomp $sql;
  next if $sql =~ /^#/;
  next if $sql =~ /^$/;
  $sql =~ s/;$//;
  !$quiet && print "# $sql\n";
  my $sth = $dbh->prepare($sql);
  $sth->execute;

  while (@data = $sth->fetchrow_array()) {
    foreach $i (0..$#data) {
      $data[$i] = "" if !defined $data[$i];
    }
    print join ("\t", @data), "\n";
  }
}

$dbh->disconnect;
