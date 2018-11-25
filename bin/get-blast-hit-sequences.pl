#!/usr/bin/perl -w
#
# take something like blast against nr
# pull out the sequences of all the hits
#
use slc454;
use slchen;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $type;
my $seq;
my @a;
my $blast_version;
my $eval_cutoff = 1e-100;
my $get_query = 0;
my $number_return = 0;	# 0 means all that pass the evalue cutoff
my $n = 0;			# to count how many we return
my $database;
my $query;
my $hitname;
my $hitcount;
my $ref;

GetOptions (
  'evalue=f' => \$eval_cutoff,
  'number=i' => \$number_return,
  'query!' => \$get_query,	# if 0, give back subject sequence (0 is default)
);

@a = <>;
foreach $i (0..$#a) {
  if ($a[$i] =~ /(.?BLAST.*)\s+(\d+\.\d+\.\d+)\s\[\w+-\d+-\d+\]/) {
    $blast_type = $1;
    $blast_version = $2;
    last;
  }
}

if ($blast_type) {
  $seq = slc454::full_blast_sequences($blast_type, @a);
  foreach $query (keys %$seq) {
    foreach $database (keys %{$seq->{$query}}) {
      $n = 0;
      foreach $hitname (keys %{$seq->{$query}->{$database}}) {
        foreach $hitcount (keys %{$seq->{$query}->{$database}->{$hitname}}) {
          $ref = $seq->{$query}->{$database}->{$hitname}->{$hitcount};
          $ref->{EVALUE} = "1" . $ref->{EVALUE} if $ref->{EVALUE} =~ /^e/;
        }
      }
      $ref = $seq->{$query}->{$database};
      foreach $hitname (sort { $ref->{$a}->{1}->{EVALUE} <=> $ref->{$b}->{1}->{EVALUE} } keys %$ref) {
        foreach $hitcount (sort {$a <=> $b} keys %{$ref->{$hitname}}) {
          if ($ref->{$hitname}->{$hitcount}->{EVALUE} <= $eval_cutoff) {
            print ">$query|$hitname|$database|$ref->{$hitname}->{$hitcount}->{EVALUE}\n";
            if ($get_query) {
              print $ref->{$hitname}->{$hitcount}->{QUERY}, "\n";
              $n++;
            } else {
              print $ref->{$hitname}->{$hitcount}->{SEQUENCE}, "\n";
              $n++;
            }
            last if $number_return > 0 && $n >= $number_return;
          }
        }
        last if $number_return > 0 && $n >= $number_return;
      }
    }
  }
}
