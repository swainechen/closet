#!/usr/bin/perl -w
#
# convert nucmer show-coords format to a dummy gcov file
# assume we only have one chromosome
#
my $coverage = 1000;
my $ref;
my $reffile;
my @r;
my $i;
my @f;
my @order = ();
my $name;
my $refstring = "";
while (<>) {
  chomp;
  @f = split /\s+/, $_;
  if (!$ref && defined $f[0] && -f $f[0]) {
    $reffile = $f[0];
    @r = `chompnewline.pl $reffile | fasta-length.pl`;
    $name = "";
    foreach $i (@r) {
      chomp $i;
      if ($i =~ /^>(.*)$/) {
        @f = split /\s+/, $1;
        $name = $f[0];
        push @order, $name;
      } elsif (length $name && $i =~ /^\d+$/) {
        $ref->{$name} = $i;
        $name = "";
      }
    }
  }
  last if /^=+$/;
}
my $cov;
foreach $i (keys %$ref) {
  $cov->{$i} = [ (0) x ($ref->{$i}+1) ];	# because we are using 1-based coordinates
}
while (<>) {
  s/^\s+//;
  @f = split /\s+/, $_;
  $refstring = $f[11] if length($f[11]);
  if (!defined $cov->{$refstring}) {
    die "Error - can't find $refstring in reference $reffile\n";
  }
  if ($f[0] =~ /\d+/ && $f[1] =~ /\d+/) {
    foreach $i ($f[0]..$f[1]) {
      $cov->{$refstring}->[$i] = $coverage;
    }
  }
}
foreach $refstring (@order) {
  if (defined $ref->{$refstring}) {
    foreach $i (1..$ref->{$refstring}) {
      print join ("\t", $refstring, $i, $cov->{$refstring}->[$i]), "\n";
    }
  }
}
