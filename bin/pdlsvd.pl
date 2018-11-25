#!/usr/bin/perl -w
#
# read in columns.  mean center.  then run an svd.
#
if (!$ARGV[0]) {
  print "Usage: ./pdlsvd.pl -i <infile> -coords <coord output file> -datacols <data cols>\n";
  print "Data cols is a comma separated list, dashes ok.  Output is average vector, basis vectors, singular values.\n";
  exit;
}
use PDL;
use PDL::Math;
use Getopt::Long;
use vars qw($u $s $v);
&Getopt::Long::Configure("pass_through");

$infile = '';
$datacols = '';
$coords = '';

GetOptions ('i=s' => \$infile, 'coords=s' => \$coords, 'datacols=s' => \$datacols);

if ($datacols eq '') {
  open F, $infile;
  $in = <F>;
  chomp $in;
  if ($in =~ /^#\s*(.*)/) {
    $datacols = $1;
  }
}
@datacols = parse_nums ($datacols);

$PDL::IO::Misc::colsep = "\t";
@cols = rcols ($infile, { EXCLUDE => '/^#/' }, @datacols );
@out = ('average vector');
foreach $i (0 .. $#datacols) {
  $avg = avg($cols[$i]);
  $cols[$i] -= $avg;
  push @out, $avg;
}
print join ("\t", @out), "\n";

$matrix = transpose (cat @cols);
($u, $s, $v) = svd ($matrix);
@coords = dog transpose $u;
wcols (@coords, $coords);

@vectors = dog transpose $v;
foreach $i (0 .. $#vectors) {
  @out = ("Basis $i", list $vectors[$i]);
  print join ("\t", @out), "\n";
}

print "Singular values\t", join ("\t", list $s), "\n";

sub parse_nums {
  my ($nums) = @_;
  chomp $nums;
  my @num = split /,/, $nums;
  my @out = ();
  foreach my $num (@num) {
    if ($num =~ /^\d+$/) {
      push @out, $num;
      next;
    }
    if ($num =~ /^(\d+)-(\d+)$/) { 
      push @out, ($1..$2);
      next;
    }
  }
  return @out;
}
