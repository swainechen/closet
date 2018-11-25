#!/usr/bin/perl -w
#
# dump data from a sql query
#
use DBI;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
$user = '';
$pass = '';
$sid = '';
$host = '';
$quiet = 0;
$single = 0;	# if true, will consider multiline input as single sql
GetOptions (
  'username=s' => \$user,
  'password=s' => \$pass,
  'host=s' => \$host,
  'sid=s' => \$sid,
  'quiet' => \$quiet,
  'single' => \$single
);
$connect_string = 'DBI:Oracle:host="' . $host . '";sid="' . $sid . '"';
my $dbh = DBI->connect($connect_string, $user, $pass, { LongReadLen => 100000000 });
$dbh->{AutoCommit} = 0;
$dbh->do("set transaction read only");
if (!@ARGV || -f $ARGV[0]) { @sql = <>; }
else { @sql = @ARGV; }

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
