#!/usr/bin/perl -w
#
# expect sam format on input
# adjust mapping position based on sam flag for mapping strand
#
use Getopt::Long;
Getopt::Long::Configure("pass_through");

# whether to check for soft clipping at beginning
# Shannon's data seems to have 5-6 bp that are still transposon in each read
my $check_soft = 0;
my $soft_values;
my $lines = 0;
my $expected = 6;
my $adjust = 0;
my $tn5dup = 9;	# tn5 duplicates 9bp when it inserts - we are getting ~10% nonspecific priming to the wrong side of the tn5. fix this by just adjusting the position for reads that map to the reverse strand - this canonically will just give the lower coordinate of the duplicated 9bp for the tn insertion site
my $print_help = 0;

GetOptions (
  "check_soft!" => \$check_soft,
  "expected_soft=i" => \$expected,
  "tndup" => \$tn5dup,
  "help" => \$print_help
);

if ($print_help) {
  print "Usage: $0 [ -check_soft|nocheck_soft ] [ -expected_soft <int> ] [ -tndup <int> ]\n";
  print "       $0 -help\n";
  print "-check_soft determines whether to check for soft clipping - needed if you didn't trim everything and the read doesn't start directly with genomic DNA. Default NOT to check for soft clipping.\n";
  print "-expected_soft if the amount that should be in the read before genomic DNA starts (if insufficient trimming)\n";
  print "-tndup is the amount of duplicated sequence upon transposon insertion - default is 9 for Tn5. Real position reported is always on the left of the duplicated insertion site regardless of insertion (or mapping) strand\n";
  exit;
}

my $name;
my $mapq;
my $hist;

while (<>) {
  next if /^$/;
  next if /^#/;
  chomp;
  $lines++;
  @f = split /\t/, $_;
  $name = $f[0];
  $flag = $f[1];
  $ref = $f[2];
  $pos = $f[3];
  $mapq = $f[4];
  $cigar = $f[5];
  $length = length($f[9]);
  if ($flag & 4) {	# by definition for pairmap this should never be true
    print STDERR "Found unmapped read, skipping...:\n$_\n";
  }
  if (!($flag & 2)) {	# by definition this is the pairmap flag
    print STDERR "Found bad pairmap read, skipping...:\n$_\n";
  }
  if (!($flag & 64)) {	# check if R1 read
    print STDERR "Found R2 read read, skipping...:\n$_\n";
  }
  if ($flag & 16) {	# reverse complemented, negative strand
#    $pos = $pos + $length - 1;
    $pos = $pos + $length - $tn5dup;
    $strand = "-";
  } else {
    $strand = "+";
  }
  if ($check_soft) {
    # on positive strand
    # position in the sam file is where the mapping starts - so soft clipped bases don't need to be subtracted
    # but if expect 6S and only have 5S, then need to add 1 to the mapping position to account for this
    if ($cigar =~ /^(\d+)S/ && $strand eq "+") {
      $s = $1;
      if ($s < $expected) { 
        $adjust = $expected - $s;
        $s = $expected;
      }
      $soft_values->{$s}++;
      $pos += $adjust;
    }
    # on negative strand, we can't correct by the full read length
    # have to subtract off the soft clip to find the real length to correct by
    # but again we probably should be correcting by the expected soft clip if the reported soft clip is less than this
    if ($cigar =~ /(\d+)S$/ && $strand eq "-") {
      $s = $1;
      if ($s < $expected) { $s = $expected; }
      $soft_values->{$s}++;
      $pos -= $s;
    }
  }
#  print join ("\t", $name, $ref, $pos, $strand, $length, $mapq, $flag, $cigar), "\n";
  $hist->{$ref}->{$pos}++;
}

foreach $ref (sort keys %$hist) {
  foreach $pos (sort {$a<=>$b} keys %{$hist->{$ref}}) {
    print join ("\t", $ref, $pos, $hist->{$ref}->{$pos}), "\n";
  }
}

if ($check_soft) {
  print STDERR "# Soft clip\tNumber (total $lines)\n";
  foreach $i (sort {$a<=>$b} keys %$soft_values) {
    print STDERR join ("\t", $i, $soft_values->{$i}), "\n";
  }
}
