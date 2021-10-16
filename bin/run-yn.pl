#!/usr/bin/perl
#
# take dna fasta file
# translate to protein, align, impose on dna
# run through yn00 program
#
use warnings;
use slchen;
use File::Temp;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $do_align = 1;	# whether to try to align sequences first
GetOptions (
  'align!' => \$do_align
);

if (!@ARGV) { print "Usage: $0 [ -align|noalign ] <infile(s)>\n"; exit; }
my $tempdir = File::Temp::tempdir( CLEANUP => 1 );

my $tempdna = $tempdir."/tempdna";
my $tempyn = $tempdir."/yn00.ctl";

my $current = `pwd`;
chomp $current;

my $f;
my @dna;
my @protein;

foreach $f (@ARGV) {
  open F, $f;
  while (<F>) {
    chomp;
    next if /^#/ || /^$/;
    push @dna, $_;
  }

  chdir $tempdir;
  my $name_hash = to_phylip_names(\@dna);

  if ($do_align) {
    # translate
    open T, "| translate.pl -noseq -nosix -notri > $tempdir/translate.tmp";
    print T join ("\n", @dna), "\n";
    close T;
    open T, "$tempdir/translate.tmp";
    @protein = ();
    while (<T>) {
      chomp;
      push @protein, $_;
    }
    unlink "$tempdir/translate.tmp";

    # align proteins
    @protein = align(@protein);
    $error = 0;
    @aligndna = paln2daln(\@protein, \@dna, \$error);
  } else {
    @aligndna = @dna;
  }

  # run yn00
  open T, "| fasta2phy.pl | phy2paml > $tempdna";
  print T join ("\n", @aligndna), "\n";
  close T;
  open Y, ">$tempyn";
  print Y "seqfile = tempdna\n";
  print Y "outfile = yn\n";
  print Y "verbose = 0\n";
  print Y "icode = 0\n";
  print Y "weighting = 0\n";
  print Y "commonf3x4 = 0\n";
  close Y;
  chdir $tempdir;
  system "yn00 > /dev/null";
  chdir $current;

  # parse the yn output file
  open R, "$tempdir/yn";
  @ynname = ();	# keep another array of orgs in case output order changed
  my ($o1, $o2, $S, $N, $t, $kappa, $omega, $dn, $dnse, $ds, $dsse);
  %dnds = ();
  while (<R>) {
    if (/^Codon position x base/) {
      while (<R>) {
        chomp;
        next if /^\s*$/;
        next if /^position\s/;
        last if /^Average/;
        s/^\s+//;
        s/\s+$//;
        push @ynname, $_;
      }
    }
#    if (/^Estimation by the method of Yang/) {
    if (/^\(B\) Yang \& Nielsen/) {
      while (<R>) {
        chomp;
        next if /^\s*$/;
        next if /^seq. seq./;
        next if /.equal weighting of pathways./;
        next if /^Yang Z, Nielsen R/;
        last if /\(C\) LWL85/;
        s/^\s*//;
        ($o1, $o2, $S, $N, $t, $kappa, $omega, $dn, undef, $dnse, $ds, undef, $dsse) = split /\s+/, $_;
        if ($omega eq '99.0000' && $dn eq '-0.0000') {
          $omega = -2;		# both dn and ds are 0
        } elsif ($omega eq '99.0000') {
          $omega = -1;		# only ds is 0
        }
        $dnds{$ynname[$o1-1].$ynname[$o2-1]} = $omega;
        $dnds{$ynname[$o2-1].$ynname[$o1-1]} = $omega;

        $numN{$ynname[$o1-1].$ynname[$o2-1]} = $N;
        $numN{$ynname[$o2-1].$ynname[$o1-1]} = $N;
        $dn{$ynname[$o1-1].$ynname[$o2-1]} = $dn * $N;
        $dn{$ynname[$o2-1].$ynname[$o1-1]} = $dn * $N;
        $dnora{$ynname[$o1-1].$ynname[$o2-1]} = $dn;
        $dnora{$ynname[$o2-1].$ynname[$o1-1]} = $dn;

        $numS{$ynname[$o1-1].$ynname[$o2-1]} = $S;
        $numS{$ynname[$o2-1].$ynname[$o1-1]} = $S;
        $ds{$ynname[$o1-1].$ynname[$o2-1]} = $ds * $S;
        $ds{$ynname[$o2-1].$ynname[$o1-1]} = $ds * $S;
        $dsora{$ynname[$o1-1].$ynname[$o2-1]} = $ds;
        $dsora{$ynname[$o2-1].$ynname[$o1-1]} = $ds;

        $dnse{$ynname[$o1-1].$ynname[$o2-1]} = $dnse;
        $dnse{$ynname[$o2-1].$ynname[$o1-1]} = $dnse;
        $dsse{$ynname[$o1-1].$ynname[$o2-1]} = $dsse;
        $dsse{$ynname[$o2-1].$ynname[$o1-1]} = $dsse;

        $t{$ynname[$o1-1].$ynname[$o2-1]} = $t;
        $t{$ynname[$o2-1].$ynname[$o1-1]} = $t;
        $kappa{$ynname[$o1-1].$ynname[$o2-1]} = $kappa;
        $kappa{$ynname[$o2-1].$ynname[$o1-1]} = $kappa;
      }
    }
  }
  @head = ();
  @out = ();
  foreach $i (sort keys %$name_hash) {
    foreach $j (sort keys %$name_hash) {
      next if $i eq $j;
      next if $i lt $j;
      push @head, $name_hash->{$i} . "-" . $name_hash->{$j};
      if (defined $dnds{$i.$j}) {
        push @out, $dnds{$i.$j};
      } else {
        push @out, "-";
      }
    }
  }
  foreach $i (sort keys %$name_hash) {
    foreach $j (sort keys %$name_hash) {
      next if $i eq $j;
      next if $i lt $j;
      push @head, $name_hash->{$i} . "-" . $name_hash->{$j};
      if (defined $dnds{$i.$j}) {
        push @out, join ("|", $numN{$i.$j}, $numS{$i.$j}, $dn{$i.$j}, $ds{$i.$j});
      } else {
        push @out, join ("|", ("-") x 4);
      }
    }
  }
  print join ("\t", $f, @head), "\n";
  print join ("\t", $f, @out), "\n";
  close R;
  chdir $current;
}
