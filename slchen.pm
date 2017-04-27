#!/usr/bin/perl
#
# some useful functions
#

use warnings;
use File::Temp;
use strict;

package slchen;
require Exporter;

use vars qw(@ISA @EXPORT @EXPORT_OK);
use vars qw(%orgcode_translate);

@ISA = qw(Exporter);
@EXPORT = qw(make_key break_key sortu parse_list make_list isfloat isint align dnaalign paln2daln align_stripgaps aa setup_hash randomize_array to_phylip_names from_phylip_names dna_distances array_mean array_median array_mode dbconnect isintergenic opengri ns_setup ns orfupstream null_bind_param do_pars_tree do_nj_tree do_ml_tree combinations gri_boxwhisker);
@EXPORT_OK = qw(make_key break_key sortu parse_list make_list isfloat isint align dnaalign paln2daln align_stripgaps aa setup_hash randomize_array to_phylip_names from_phylip_names dna_distances array_mean array_median array_mode dbconnect isintergenic opengri ns_setup ns orfupstream null_bind_param do_pars_tree do_nj_tree do_ml_tree combinations gri_boxwhisker);

sub dbconnect {
  my ($database, $host, $user, $pass) = @_;

  # default values first
  $database = '';
  $host = '';
  $user = '';
  $pass = '';

  my $dbh;
  $dbh = DBI->connect('DBI:mysql:database='.$database.';host='.$host, $user, $pass, { LongReadLen => 100000000 });
  return $dbh;
}

sub make_key {
  my @k = @_;
  return (join ("~~", @k));
}

sub break_key {
  my ($k) = @_;
  return (split /~~/, $k);
}

# sort unique, be smart about numbers and non numbers
sub sortu {
  my @a = @_;
  my @defined = ();
  my @undefined = ();
  my $number = 1;
  foreach my $i (0..$#a) {
    if (!defined $a[$i]) {
      push @undefined, $i;
    }
    push @defined, $i;
    if (!isfloat($a[$i])) {
      $number = 0;
    }
  }
  if ($number) {
    @a = sort { $a <=> $b } @a[@defined];
  } else {
    @a = sort { $a cmp $b } @a[@defined];
  }
  foreach my $i (reverse 1..$#a) {
    if ($number) { 
      if ($a[$i] == $a[$i-1]) { splice @a, $i, 1; }
    } else {
      if ($a[$i] eq $a[$i-1]) { splice @a, $i, 1; }
    }
  }
  foreach my $i (@undefined) {
    push @a, undef;
  }
  return @a;
}

# take 1,2,3,5-9,14..18,22 and give array with all numbers
sub parse_list {
  my ($list) = @_;
  my @numbers = ();
  my @temp = ();
  chomp $list;
  my @atoms = split /,/, $list;
  foreach my $i (0..$#atoms) {
    if (isfloat($atoms[$i])) {
      push @numbers, $atoms[$i];
    } elsif ($atoms[$i] =~ /\.\./) {
      @temp = split /\.\./, $atoms[$i];
      if (isint($temp[0]) && isint($temp[1])) {
        if ($temp[0] < $temp[1]) {
          push @numbers, ($temp[0]..$temp[1]);
        } else {
          push @numbers, reverse ($temp[1]..$temp[0]);
        }
      }
    } elsif ($atoms[$i] =~ /-/) {
      @temp = split /-/, $atoms[$i];
      if (isint($temp[0]) && isint($temp[1])) {
        if ($temp[0] < $temp[1]) {
          push @numbers, ($temp[0]..$temp[1]);
        } else {
          push @numbers, reverse ($temp[1]..$temp[0]);
        }
      }
    }
  }
  return @numbers;
}

# take list of integers, give a shortened list of ranges
sub make_list {
  my (@a) = @_;
  my $ranger = "..";
  my $separater = ",";
  my $i;
  my @return;
  foreach $i (0..$#a) {
    if (!isint($a[$i])) {
      print STDERR "Noninteger $a[$i] found in make_list call\n";
      return @a;
    }
  }
  @a = sortu (@a);
  my $start;
  my $end;
  undef $start;
  undef $end;
  foreach $i (0..$#a) {
    if (!defined $start) {
      $start = $a[$i];
      $end = $a[$i];
      next;
    }
    if ($a[$i] == $end + 1) {
      $end = $a[$i];
    } else {
      if ($end == $start) {
        push @return, $start;
      } elsif ($end == $start+1) {
        push @return, $start, $end;
      } else {
        push @return, join ($ranger, $start, $end);
      }
      $start = $a[$i];
      $end = $a[$i];
    }
  }
  if (defined $start && defined $end) {
    if ($end == $start) {
      push @return, $start;
    } elsif ($end == $start+1) {
      push @return, $start, $end;
    } else {
      push @return, join ($ranger, $start, $end);
    }
  }
  return join ($separater, @return);
}

sub isfloat {
  my ($string) = @_;
  if (!defined $string) {
    return 0;
  }
  chomp $string;
  if ($string =~ /^([+-])?(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
    return 1;
  } else {
    return 0;
  }
}

sub isint {
  my ($string) = @_;
  if (!defined $string) {
    return 0;
  }
  chomp $string;
  if ($string =~ /^[+-]?\d+$/) {
    return 1;
  } else {
    return 0;
  }
}

sub align {
  # take a fasta array of protein sequence, make sure fasta names are short
  # run it through clustalw
  # return the results
  my @f = @_;
  my @return = ();
  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );
  my $tempfile = $tempdir."/clustal-input";
  my $clustalout = $tempfile.".phy";
#  my $clustalout = $tempfile.".aln";
  my $namehash = to_phylip_names(\@f);
  my $in;
  open T, ">$tempfile";
  foreach my $f (@f) {
    chomp $f;
    next if $f =~ /^#/;
    next if $f =~ /^$/;
    print T $f, "\n";
  }
  close T;
  if (scalar @f > 2) {
    system "clustalw -output=PHYLIP -type=PROTEIN -infile=$tempfile > /dev/null"
;
    open T, "phy2fasta -interleaved $clustalout |";
  } else {
    open T, $tempfile;
  }
  while ($in = <T>) {
    chomp $in;
    next if $in =~ /^#/;
    next if $in =~ /^$/;
    push @return, $in;
  }
  close T;
  &from_phylip_names (\@return, $namehash);
  unlink $tempfile;
  unlink "$tempfile.dnd";
  unlink $clustalout;
  system "rm -rf $tempdir";
  return @return;
}

sub dnaalign {
  # take a fasta array of dna sequence
  # run it through clustalw
  # return the results as an array
  my (@g) = @_;
  my @f = ();
  my @return = ();
  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );
  my $tempfile = $tempdir."/clustal-input";
#  my $clustalout = $tempfile.".phy";
  my $clustalout = $tempfile.".aln";
  my $in;
  if (ref($g[0]) eq "ARRAY") {
    @f = @{$g[0]};
  } elsif (ref($g[0]) eq "HASH") {
    @f = hash2fasta($g[0]);
  } else {
    @f = @g;
  }
  open T, ">$tempfile";
  my $namehash = to_phylip_names(\@f);
  foreach my $f (@f) {
    chomp $f;
    next if $f =~ /^#/;
    next if $f =~ /^$/;
    print T $f, "\n";
  }
  close T;
  if (scalar @f > 2) {
#    system "clustalw -output=PHYLIP -type=DNA -infile=$tempfile > /dev/null"
#    system "clustalw -type=DNA -infile=$tempfile -gapopen=1000 > /dev/null"
    system "clustalw -type=DNA -infile=$tempfile > /dev/null"
;
    open T, "aln2phy $clustalout | phy2fasta -interleaved |";
  } else {
    open T, $tempfile;
  }
  while ($in = <T>) {
    chomp $in;
    next if $in =~ /^#/;
    next if $in =~ /^$/;
    push @return, $in;
  }
  close T;
  &from_phylip_names (\@return, $namehash);
  unlink $tempfile;
  unlink "$tempfile.dnd";
  unlink $clustalout;
  system "rm -rf $tempdir";
  return @return;
}

sub paln2daln {
  # take an aligned fasta array of aligned protein sequence
  # also need fasta array of dna sequence with the same headers
  # return a DNA version with codons kept together
  # look for the gid in the second field of the header line
  my ($pref, $dref, $errorref) = @_;
  $$errorref = 0;
  my @prot = @$pref;
  my @dna = @$dref;
  my (@return, $pos, %dseq, $name, $dseq);
  foreach my $i (0..$#dna) {
    chomp $dna[$i];
    if ($dna[$i] =~ /^>/) {
      $name = $dna[$i];
    } else {
      $dseq{$name} = $dna[$i];
    }
  }
  foreach my $i (0..$#prot) {
    chomp $prot[$i];
    if ($prot[$i] =~ /^>/) {
      $name = $prot[$i];
      $return[$i] = $name;
    } else {
      my $pseq = $prot[$i];
      if (!defined $dseq{$name}) {
        die "Can't find DNA sequence for $name\n";
      }
      $dseq = $dseq{$name};
      my @p = split //, $pseq;
      $return[$i] = "";
      $pos = 0;
      my $d;
      foreach my $p (@p) {
        if (length($dseq) >= $pos + 3) {
          $d = substr ($dseq, $pos, 3);
        } else {
          $d = '---';
        }
        if ($p eq '-') {
          $return[$i] .= "---";
          if (aa($d) eq 'X') {
            $pos += 3;
          }
        } else {
          if ($d =~ /[^GATCgatc]/ && $p eq 'X') {
            print STDERR "N in DNA, X in protein at sequence $name, pos $pos, dna $d, prot $p. Converting DNA to lowercase\n";
            $return[$i] .= lc($d);
            $$errorref += 1;
          } elsif (aa($d) eq 'X') {
            print STDERR "Nonsense at sequence $name, pos $pos, dna $d, prot $p, trans ", aa($d), ", converting DNA to lowercase\n";
            $return[$i] .= lc($d);
            $$errorref += 1;
          } elsif (aa($d) ne uc $p) {
            if ($d =~ m/[^GATCgatc]/) {
              print STDERR "Translation error at sequence $name, pos $pos, dna $d, prot $p, trans ", aa($d), ", converting non-GATC to -\n";
              $d =~ s/[^GATCgatc]/-/g;
            } else {
              print STDERR "Translation error at sequence $name, pos $pos, dna $d, prot $p, trans ", aa($d), ", converting DNA to lowercase\n";
            }
            $return[$i] .= lc($d);
            $$errorref += 1;
          } else {
            $return[$i] .= $d;
          }
          $pos += 3;
        }
      }
    }
  }
  return @return;
}

sub align_stripgaps {
  # take aligned fasta sequence
  # strip out gaps from alignment
  my (@s) = @_;
  my ($i, $j, $k);
  my $sref;
  my $name;
  my @a;
  my @return;
  $j = 0;
  foreach $i (0..$#s) {
    if ($s[$i] =~ /^>/) {
      $sref->[$j]->{NAME} = $s[$i];
    } else {
      $sref->[$j]->{SEQ} = $s[$i];
      $j++;
    }
  }
  foreach $j (0..$#$sref) {
    @a = split //, $sref->[$j]->{SEQ};
    foreach $i (reverse 0..$#a) {
      if ($a[$i] eq '-') {
        foreach $k (0..$#$sref) {
          substr ($sref->[$k]->{SEQ}, $i, 1) = "";
        }
      }
    }
  }
  foreach $j (0..$#$sref) {
    push @return, $sref->[$j]->{NAME};
    push @return, $sref->[$j]->{SEQ};
  }
  return @return;
}

# take orgcode, references to hashes, fill them
# hashes are:
#   gid2triv, gid2syst, gid2annot, gid2faa, gid2ffn, gid2loc, gid2org, loc2ffn
sub setup_hash {
  my ($org, $gid2triv_r, $gid2syst_r, $gid2annot_r, $gid2faa_r, $gid2ffn_r, $gid2loc_r, $gid2org_r, $loc2ffn_r) = @_;
  my (@ptt, @faa, @ffn) = ();
  my (@subseq) = ();
  my ($gi, $i, $j, @f, @g, $start, $stop, $strand, $switch, $loc, @sub);
  @ptt = `print-genome-file.pl $org -suf ptt`;
  @faa = `print-genome-file.pl $org -suf faa | chompnewline.pl`;
  @ffn = `print-genome-file.pl $org -suf ffn | chompnewline.pl`;
  foreach my $p (@ptt) {
    chomp $p;
    if ($p =~ /^\s*(\d+)\.\.(\d+)\s+/) {
      ($start, $stop) = ($1, $2);
      @f = split /\t/, $p;
      $strand = $f[1];
      $gi = $f[3];
      next if ($gi < 0 || $gi == 9999);
      if (defined $gid2triv_r) {
        $$gid2triv_r{$gi} = $f[4];	# this is usually ok
      }
      if (defined $gid2syst_r) {
        $$gid2syst_r{$gi} = $f[5];
      }
      if (defined $gid2annot_r) {
        $$gid2annot_r{$gi} = $f[8];
      }
      if (defined $gid2loc_r) {
        $$gid2loc_r{$gi} = "$org-$start..$stop";
      }
      if (defined $gid2org_r) {
        $$gid2org_r{$gi} = $org;
      }
      if ($strand =~ /-/) {
        $switch = $start;
        $start = $stop;
        $stop = $switch;
      }
      push @subseq, ">~~$gi~~\n$start..$stop\n";
    }
  }
  foreach my $s (@faa) {
    chomp $s;
    if ($s =~ /^>gi\|(\d+)\|/) {
      $gi = $1;
    } else {
      if (defined $gid2faa_r) {
        $$gid2faa_r{$gi} = $s;
      }
    }
  }
  if (defined $loc2ffn_r) {
    foreach my $s (@ffn) {
      chomp $s;
      if ($s =~ /^>/) {
        @f = split /:/, $s;
        $f[1] =~ s/\s.*//;
        $f[1] =~ s/\(//g;
        $f[1] =~ s/\)//g;
        $f[1] =~ s/c//g;
        $f[1] =~ s/\s//g;
        $f[1] =~ s/>//g;
        @g = split /[-,]/, $f[1];
        @g = sort { $a <=> $b } @g;
        $loc = "$org-$g[0]..$g[$#g]";
      } else {
        $$loc2ffn_r{$loc} = $s;
      }
    }
  }
  if (defined $gid2ffn_r) {
    my $temp = "/tmp/pull.temp.".rand();
    while (-f $temp) {
      $temp = "/tmp/pull.temp.".rand();
    }
    open SUBSEQ, "| subseq.pl $org > $temp";
    print SUBSEQ @subseq;
    close SUBSEQ;
    open SUBSEQ, $temp;
    @sub = <SUBSEQ>;
    close SUBSEQ;
    unlink ($temp);
    $gi = -1;
    foreach $j (0..$#sub) {
      chomp $sub[$j];
      if ($sub[$j] =~ /^>~~(\d+)~~$/) {
        $gi = $1;
      } elsif ($sub[$j] !~ /^>/) {	# if frameshift, use genbank ffn
        if ($gi > 0 && $$gid2loc_r{$gi} && $$loc2ffn_r{$$gid2loc_r{$gi}}
            && $$loc2ffn_r{$$gid2loc_r{$gi}} ne $sub[$j]) {
          $$gid2ffn_r{$gi} = $$loc2ffn_r{$$gid2loc_r{$gi}};
        } else {
          $$gid2ffn_r{$gi} = $sub[$j];
        }
        $gi = -1;
      }
    }
  }
}

sub translate {
  # translate sequence to aa
  # use frame if we have it
  my $sequence = shift;
  my $frame;
  if (scalar @_) {
    $frame = shift;
  }
  $sequence = uc($sequence);
  $frame = 1 if !defined $frame;
  if ($frame < 0) {
    $sequence = revcomp($sequence);
  }
  $frame = abs($frame);
  $sequence = substr($sequence, $frame-1);
  my $i;
  my $return = "";
  for ($i = 0; $i < length $sequence; $i += 3) {
    $return .= aa(substr($sequence, $i, 3));
  }
  return $return;
}

sub revcomp {
  my ($inline) = $_[0];
  my $outline = reverse ($inline);
  $outline =~ tr/ABCDGHKMNRSTVWXYabcdghkmnrstvwxy/TVGHCDMKNYSABWXRtvghcdmknysabwxr/;
  return $outline;
}

# translate codon to aa
sub aa {
  my ($d) = @_;
  $d = uc $d;
  if ($d eq "TTT") { return "F"; }
  if ($d eq "TTC") { return "F"; }
  if ($d eq "TTY") { return "F"; }
  if ($d eq "TTA") { return "L"; }
  if ($d eq "TTG") { return "L"; }
  if ($d eq "TTR") { return "L"; }
  if ($d eq "TCT") { return "S"; }
  if ($d eq "TCC") { return "S"; }
  if ($d eq "TCA") { return "S"; }
  if ($d eq "TCG") { return "S"; }
  if ($d =~ /^TC/) { return "S"; }
  if ($d eq "TAT") { return "Y"; }
  if ($d eq "TAC") { return "Y"; }
  if ($d eq "TAY") { return "Y"; }
  if ($d eq "TAA") { return "X"; }
  if ($d eq "TAG") { return "X"; }
  if ($d eq "TGT") { return "C"; }
  if ($d eq "TGC") { return "C"; }
  if ($d eq "TGY") { return "C"; }
  if ($d eq "TGA") { return "X"; }
  if ($d eq "TGG") { return "W"; }
  if ($d eq "CTT") { return "L"; }
  if ($d eq "CTC") { return "L"; }
  if ($d eq "CTA") { return "L"; }
  if ($d eq "CTG") { return "L"; }
  if ($d =~ /^CT/) { return "L"; }
  if ($d eq "YTA") { return "L"; }
  if ($d eq "YTG") { return "L"; }
  if ($d eq "YTR") { return "L"; }
  if ($d eq "CCT") { return "P"; }
  if ($d eq "CCC") { return "P"; }
  if ($d eq "CCA") { return "P"; }
  if ($d eq "CCG") { return "P"; }
  if ($d =~ /^CC/) { return "P"; }
  if ($d eq "CAT") { return "H"; }
  if ($d eq "CAC") { return "H"; }
  if ($d eq "CAY") { return "H"; }
  if ($d eq "CAA") { return "Q"; }
  if ($d eq "CAG") { return "Q"; }
  if ($d eq "CAR") { return "Q"; }
  if ($d eq "CGT") { return "R"; }
  if ($d eq "CGC") { return "R"; }
  if ($d eq "CGA") { return "R"; }
  if ($d eq "CGG") { return "R"; }
  if ($d =~ /^CG/) { return "R"; }
  if ($d eq "MGA") { return "R"; }
  if ($d eq "MGG") { return "R"; }
  if ($d eq "MGR") { return "R"; }
  if ($d eq "ATT") { return "I"; }
  if ($d eq "ATC") { return "I"; }
  if ($d eq "ATA") { return "I"; }
  if ($d eq "ATY") { return "I"; }
  if ($d eq "ATW") { return "I"; }
  if ($d eq "ATM") { return "I"; }
  if ($d eq "ATH") { return "I"; }
  if ($d eq "ATG") { return "M"; }
  if ($d eq "ACT") { return "T"; }
  if ($d eq "ACC") { return "T"; }
  if ($d eq "ACA") { return "T"; }
  if ($d eq "ACG") { return "T"; }
  if ($d =~ /^AC/) { return "T"; }
  if ($d eq "AAT") { return "N"; }
  if ($d eq "AAC") { return "N"; }
  if ($d eq "AAY") { return "N"; }
  if ($d eq "AAA") { return "K"; }
  if ($d eq "AAG") { return "K"; }
  if ($d eq "AAR") { return "K"; }
  if ($d eq "AGT") { return "S"; }
  if ($d eq "AGC") { return "S"; }
  if ($d eq "AGY") { return "S"; }
  if ($d eq "AGA") { return "R"; }
  if ($d eq "AGG") { return "R"; }
  if ($d eq "AGR") { return "R"; }
  if ($d eq "GTT") { return "V"; }
  if ($d eq "GTC") { return "V"; }
  if ($d eq "GTA") { return "V"; }
  if ($d eq "GTG") { return "V"; }
  if ($d =~ /^GT/) { return "V"; }
  if ($d eq "GCT") { return "A"; }
  if ($d eq "GCC") { return "A"; }
  if ($d eq "GCA") { return "A"; }
  if ($d eq "GCG") { return "A"; }
  if ($d =~ /^GC/) { return "A"; }
  if ($d eq "GAT") { return "D"; }
  if ($d eq "GAC") { return "D"; }
  if ($d eq "GAY") { return "D"; }
  if ($d eq "GAA") { return "E"; }
  if ($d eq "GAG") { return "E"; }
  if ($d eq "GAR") { return "E"; }
  if ($d eq "GGT") { return "G"; }
  if ($d eq "GGC") { return "G"; }
  if ($d eq "GGA") { return "G"; }
  if ($d eq "GGG") { return "G"; }
  if ($d =~ /^GG/) { return "G"; }
  return "?";
}

sub randomize_array {
  my (@a) = @_;
  my @index;
  foreach my $i (0..$#a) {
    $index[$i] = rand();
  }
  return (@a[sort { $index[$a] <=> $index[$b] } (0..$#a)]);
}

sub to_phylip_names {
  # return fasta dna sequence with names that are short enough (< 10 char) for 
  # phylip format
  # we're going to use numbers preceded by a hopefully unique string
  # this gives us a capacity of 10^6 or so
  # return a hash reference that can be used to convert back also
  # try to take both arrays and hashes
  # usage: $name_hash = to_phylip_names(\@dna_array);
  # or:    $name_hash = to_phylip_names(\%dna_hash);
  my ($dna_ref) = @_; 
  my $name;
  my $i = 0;
  my $order = 0;
  my $j;
  my $string;
  if (ref($dna_ref) eq "ARRAY") {
    foreach $j (0..$#{$dna_ref}) {
      $i++ if $dna_ref->[$j] =~ /^>/;
    }
    $order = length($i);
    $i = 0;
    foreach $j (0..$#{$dna_ref}) {
      chomp $dna_ref->[$j]; 
      if ($dna_ref->[$j] =~ /^>(.*)$/) {
        $i++;
        $string = name_key($i, $order);
        $name->{$string} = $1;
        $dna_ref->[$j] = ">$string";
      }
    }
  } elsif (ref($dna_ref) eq "HASH") {
    $order = length(scalar(keys(%$dna_ref)));
    foreach $j (keys %$dna_ref) {
      $i++;
      $string = name_key($i, $order);
      $name->{$string} = $j;
      $dna_ref->{$string} = $dna_ref->{$j};
      delete $dna_ref->{$j};
    }
  }
  return ($name); 

  sub name_key {
    my ($i, $order) = @_;
    return sprintf("__P__%." . $order . "d", $i);
  }
}

sub from_phylip_names {
  # take fasta sequence named by to_phylip_names and the name hash
  # change the sequences back
  # usage: &from_phylip_names(\@dna_array, $name_hash);
  # or:    &from_phylip_names(\%dna_hash, $name_hash);
  # actually hope the names are unique so we can take anything here, like trees
  # so:    $string1 = from_phylip_names($string2, $name_hash);
  my ($dna_ref, $name) = @_;
  my $j = 0;
  if (ref($dna_ref) eq "ARRAY") {
    foreach $j (0..$#{$dna_ref}) {
      chomp $dna_ref->[$j];
      if ($dna_ref->[$j] =~ /^>(.*)$/ && defined $name->{$1}) {
        $dna_ref->[$j] = ">" . $name->{$1};
      }
    }
  } elsif (ref($dna_ref) eq "HASH") {
    foreach $j (keys %$dna_ref) {
      if (defined $name->{$j}) {
        $dna_ref->{$name->{$j}} = $dna_ref->{$j};
        delete $dna_ref->{$j};
      }
    }
  } else {
    foreach $j (keys %$name) {
      $dna_ref =~ s/$j/$name->{$j}/;
    }
    return $dna_ref;
  }
}

sub dna_distances {
  my @dna = @_;
  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );
  my $current = `pwd`;
  chomp $current;
  chdir ($tempdir);
  my $control = "phylip.cmd";
  my $dnafile = "phylip.dna";
  open DNA, "| fasta2phy.pl > $dnafile";
  my $namehash = to_phylip_names (\@dna);
  print DNA join ("\n", @dna);
  close DNA;
  open CMD, ">$control";
  print CMD "$dnafile\n";
  print CMD "l\n";      # lower triangular output
  print CMD "y\n";      # start the program
  close CMD;

  if (-f "/usr/bin/phylip") {
    system "phylip dnadist < $control > /dev/null";
  } else {
    system "dnadist < $control > /dev/null";
  }

  open O, "outfile";
  my @out = <O>;
  close O;
  my $i;
  my $j;
  my @f;
  my $dist;
  my @name;
  foreach $i (1..$#out) {       # first line should be # of sequences
    @f = split /\s+/, $out[$i];
    $name[$i] = $f[0];
    foreach $j (1..$#f) {       # make a symmetric distance hash using original names
      $dist->{$namehash->{$f[0]}}->{$namehash->{$name[$j]}} = $f[$j];
      $dist->{$namehash->{$name[$j]}}->{$namehash->{$f[0]}} = $f[$j];
    }
  }

  chdir $current;
  system "rm -rf $tempdir";

  return $dist;
}

sub array_mean {
  # we're going to bascially throw out non number values
  my @a = @_;
  my ($t, $mean, $count);
  $mean = 0;
  $count = 0;
  foreach $t (@a) {
    if (isfloat($t)) {
      $mean += $t;
      $count++;
    }
  }
  if (!$count) { return 0; }
  return ($mean/$count);
}

sub array_median {
  # again, throwing out non number values
  my @a = @_;
  my ($t, $i, $j, @b);
  foreach $t (@a) {
    if (isfloat($t)) {
      push @b, $t;
    }
  }
  @b = sort { $a <=> $b } @b;
  return $b[0] if scalar @b == 1;
  return undef if !scalar @b;
  $i = int $#b/2;
  $j = int ($#b+1)/2;
  return ($b[$i] + $b[$j])/2;
}

sub array_mode {
  my @a = @_;
  my %count;
  my $t;
  foreach $t (@a) {
    if ($count{$t}) { $count{$t}++; }
    else { $count{$t} = 1; }
  }
  my $max = 0; 
  foreach $t (keys %count) {
    if ($count{$t} > $max) { $max = $count{$t}; }
  }
  my @max;
  foreach $t (keys %count) {
    if ($count{$t} == $max) {
      push @max, $t;
    }
  }
  return (array_mean(@max));
}

sub array_stdev {
  my @a = @_;
  my $total = 0;
  my $mean = array_mean(@a);
  my $i = 0;
  foreach my $a (@a) {
    next if !isfloat($a);
    $total += ($mean - $a) ** 2;
    $i++;
  }
  if ($i > 1) {
    return ($total/($i-1)) ** 0.5;
  } else {
    return 0;
  }
}

sub isintergenic {
  # return 1 if intergenic, 0 if not
  # require a setup first, then can reuse this hash
  # assume we have Orgmap and $orgcode and $pttfile and $sequence defined
  my ($ref, $pos) = @_;
  $pos = abs($pos);
  my $i;
  my $line;
  if (!length $Orgmap::sequence) {
    &Orgmap::read_sequence;
  }
  while ($pos > length $Orgmap::sequence) {
    $pos -= length $Orgmap::sequence;
  }
  if (!scalar keys %$ref) {
    foreach $i (1..length $Orgmap::sequence) {
      $ref->{$i} = 1;
    }
    open PTT, $Orgmap::pttfile;
    while ($line = <PTT>) {
      if ($line =~ /^\s*(\d+)\.\.(\d+)\s+/) {
        foreach $i ($1..$2) {
          $ref->{$i} = 0;
        }
      }
    }
    close PTT;
  }
  return $ref->{$pos};
}

sub orfupstream {
  # return genbankids we are upstream of if intergenic
  # do this in a hash with keys LEFT and RIGHT
  # return empty list if inside a gene
  # require a setup first, then can reuse this hash
  # assume we have Orgmap and $orgcode and $pttfile and $sequence defined
  my ($ref, $pos) = @_;
  $pos = abs($pos);
  my $i;
  my @f;
  my $line;
  my $r = {};
  if (!defined $Orgmap::sequence) {
    &Orgmap::read_sequence;
  }
  while ($pos > length $Orgmap::sequence) {
    $pos -= length $Orgmap::sequence;
  }
  if (!scalar keys %$ref) {
    print "  --Remake ptt hash...\n";
    foreach $i (1..length $Orgmap::sequence) {
      $ref->{$i} = 0;
    }
    open PTT, $Orgmap::pttfile;
    while ($line = <PTT>) {
      if ($line =~ /^\s*(\d+)\.\.(\d+)\s+/) {
        my ($s, $e) = ($1, $2);
        @f = split /\t/, $line;
        next if $f[3] <= 0;	# can't deal with RNAs right now
        $f[1] =~ s/\s//g;
        my $tag = $f[1].$f[3];
        foreach $i ($s..$e) {
          $ref->{$i} = $tag;
        }
      }
    }
    close PTT;
    print "done\n";
  }

  if (!$ref->{$pos}) {	# if we're intergenic
    # forward direction, looking for genes in + direction
    foreach $i ($pos..length($Orgmap::sequence)) {
      if ($ref->{$i}) {
        if (substr($ref->{$i}, 0, 1) eq '+') {
          $r->{RIGHT} = abs $ref->{$i};
        }
        last;
      }
    }
    # going backwards, look for genes in - direction
    foreach $i (reverse 1..$pos) {
      if ($ref->{$i}) {
        if (substr($ref->{$i}, 0, 1) eq '-') {
          $r->{LEFT} = abs $ref->{$i};
        }
        last;
      }
    }
  }
  return $r;
}

sub opengri {
  # take a filename
  # get a GRI handle
  # take care of some options also
  my ($filename, $options) = @_;
  my $o;
  if ($filename !~ /\.ps$/) {
    $filename .= ".ps";
  }
  open GRI, "| gri -batch -nowarn_offpage -output $filename";
  foreach $o (keys %$options) {
  }
  return *GRI;
}

sub ns_setup {
  # read ptt file, etc
  my ($s, $e, $d, $gi, @f);
  my $ptt;
  open PTT, $Orgmap::pttfile;
  while (<PTT>) {
    if (/^\s*(\d+)\.\.(\d+)\s+/) {
      ($s, $e) = ($1, $2);
      chomp;
      @f = split /\t/, $_;
      $gi = $f[3];
      $f[1] =~ s/\s//g;
      $ptt->{$gi}->{SYST} = $f[5];
      $ptt->{$gi}->{STRAND} = $f[1];
      $ptt->{$gi}->{START} = $s;
      $ptt->{$gi}->{END} = $e;
      $ptt->{$gi}->{CODE} = $f[6];
    }
  }
  close PTT;
  &Orgmap::read_sequence;
  return ($ptt);
}

sub ns {
  # we'll return hash reference with fields:
  # INTERGENIC
  # SYNONYMOUS
  # ORIGINALAA
  # ORIGINALCODON
  # NEWAA
  # NEWCODON
  # CODONPOSITION
  # SYSTEMATIC
  # GID

  my ($p, $nt, $ptt) = @_;
  my $tempnt = $nt;
  my $return = {};
  # first check if it's intergenic
  my @gi = ();
  my $gi;
  my ($seq, $coord, $pos, $oldcod, $newcod);
  my ($oldaa, $newaa);
  foreach $gi (keys %$ptt) {
    if ($p >= $ptt->{$gi}->{START} && $p <= $ptt->{$gi}->{END}) {
      push @gi, $gi;
    }
  }
  if (!scalar @gi) {
    $return->{INTERGENIC} = 1;
#    return (1, undef, undef, undef);
  } else {
    foreach $gi (@gi) {	# right now this only does the first of these
      if ($ptt->{$gi}->{STRAND} eq '+') {
        $seq = Orgmap::subseq($ptt->{$gi}->{START}, $ptt->{$gi}->{END} - $ptt->{$gi}->{START} + 1);
        $coord = int (($p - $ptt->{$gi}->{START})/3);
        $pos = $p - $ptt->{$gi}->{START} - 3 * $coord;
        $oldcod = uc substr($seq, 3*$coord, 3);
        $newcod = $oldcod;
        substr ($newcod, $pos, 1) = uc $tempnt;
      } else {
        # the subseq function will take care of doing the reverse complement
        $seq = Orgmap::subseq(-$ptt->{$gi}->{END}, $ptt->{$gi}->{END} - $ptt->{$gi}->{START} + 1);
        $coord = int (($ptt->{$gi}->{END} - $p)/3);
        $pos = $ptt->{$gi}->{END} - $p - 3 * $coord;
        $oldcod = uc substr($seq, 3*$coord, 3);
        $newcod = $oldcod;
        $tempnt =~ tr/gatcGATC/ctagCTAG/;
        substr ($newcod, $pos, 1) = uc $tempnt;
      }
      $oldaa = aa($oldcod);
      $newaa = aa($newcod);
      if ($oldaa eq $newaa) {
        $return->{INTERGENIC} = 0;
        $return->{SYNONYMOUS} = 1;
        $return->{ORIGINALAA} = $oldaa;
        $return->{NEWAA} = $newaa;
        $return->{ORIGINALCODON} = $oldcod;
        $return->{NEWCODON} = $newcod;
        $return->{CODONPOSITION} = $pos + 1;
        $return->{SYSTEMATIC} = $ptt->{$gi}->{SYST};
        $return->{GID} = $gi;
#        return (0, 1, $oldcod, $oldaa, $newcod, $newaa, $pos + 1);
      } else {
        $return->{INTERGENIC} = 0;
        $return->{SYNONYMOUS} = 0;
        $return->{ORIGINALAA} = $oldaa;
        $return->{NEWAA} = $newaa;
        $return->{ORIGINALCODON} = $oldcod;
        $return->{NEWCODON} = $newcod;
        $return->{CODONPOSITION} = $pos + 1;
        $return->{SYSTEMATIC} = $ptt->{$gi}->{SYST};
        $return->{GID} = $gi;
#        return (0, 0, $oldcod, $oldaa, $newcod, $newaa, $pos + 1);
      }
    }
  }
  return $return;
}

sub null_bind_param {
  my ($sth, $field, $param) = @_;
  if (defined $param) {
    $sth->bind_param($field, $param);
  } else {
    $sth->bind_param($field, undef);
  }
}

sub sql_statement_replace {
  # take sql statement from dbi, i.e. with ? for values
  # substitute in values to get a valid SQL statement
  my ($sql) = shift;
  my @args = @_;
  my $num_slots = $sql =~ tr/?/?/;
  die "Value/slot mismatch, ", scalar @args, " values, $num_slots slots, statement $sql\n" if scalar @args != $num_slots;
  my $a;
  foreach $a (@args) {
    if (!defined $a) {
      $sql =~ s/\?/NULL/;
    } elsif (isfloat $a) {
      $sql =~ s/\?/$a/;
    } else {
      $sql =~ s/\?/'$a'/;
    }
  }
  return $sql;
}

sub get_root_name {
  # this is a support procedure for the do_*_tree procedures
  # if we're rooting, we need to know the name of the outgroup
  my $root = shift;
  my (@dna) = @_;
  my ($i, $j);
  $j = 0;
  foreach $i (0..$#dna) {
    if ($dna[$i] =~ s/^>//) {
      $j++;
      if ($j == $root) {
        return ($dna[$i]);
      }
    }
  }
}

sub tree_oneline {
  # this is a support procedure for the do_*_tree procedures
  # fix the trees so that they are one per line
  my (@t) = @_;
  my ($i, $j);
  my @r;
  $j = 0;
  foreach $i (0..$#t) {
    chomp $t[$i];
    if (!defined $r[$j]) {
      $r[$j] = $t[$i];
    } else {
      $r[$j] .= $t[$i];
    }
    if ($t[$i] =~ /;(?:\[\d\.\d+\])?$/) {
      $j++;
    }
  }
  return @r;
}

sub prune_tree {
  # this is a support procedure for the do_*_tree procedures
  # if we're rooting, we need to get the outgroup out of the tree to get a
  # fully rooted tree
  my $root_name = shift;
  my (@t) = @_;
  my $i;
  # the root_name should be at the top level, by itself
  # so it should precede or be followed by a comma
  foreach $i (0..$#t) {
    if (!($t[$i] =~ s/,$root_name(?::\d\.\d+)//) &&
        !($t[$i] =~ s/$root_name(?::\d\.\d+),//) ) {
      print STDERR "slchen.pm: sub prune_tree didn't find $root_name in $t[$i]\n";
    }
  }
  return @t;
}

sub do_pars_tree {
  # first parameter is whether we root or not (0 = no root, >0 = root)
  # we will root on the sequence specified by root (1-based, like phylip)
  # we take fasta sequence, it should already be aligned
  # we have to assume that the names are phylip-compatible already
  my $root = shift;
  my @dna = @_;

  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );

  my $tempdna = $tempdir."/infile";
  my $tempctl = $tempdir."/phylip.ctl";
  my $current = `pwd`;
  chomp $current;
  my $i;
  my $root_name;
  my @trees;
  my $temptree;
  my $seed;

  # make sure we control newlines properly
  foreach $i (0..$#dna) {
    chomp $dna[$i];
  }

  open T, "| fasta2phy.pl > $tempdna";
  print T join ("\n", @dna);
  close T;
  open T, ">$tempctl";
  $seed = int (rand(8192)) * 4 + 1;
  print T "j\n$seed\n1000\n";	# randomize input order with $seed, 1000 reps
  print T "2\n3\n";		# no progress output, no graphical tree
  if ($root) {
    print T "o\n$root\n";	# root the tree
  }
  print T "y\n";		# run the program
  close T;
  chdir $tempdir;
  if (-f "/usr/bin/phylip") {
    system "phylip dnapars < $tempctl > /dev/null";
  } else {
    system "dnapars < $tempctl > /dev/null";
  }
  $temptree = $tempdir."/outtree";
  open T, $temptree;
  @trees = <T>;
  close T;
  chdir $current;
  system "rm -rf $tempdir";
  @trees = tree_oneline(@trees);

  # if we're rooting, get the name of the outgroup
  # then clean out the tree to get a fully rooted tree
  if ($root) {
    $root_name = get_root_name($root, @dna);
    @trees = prune_tree($root_name, @trees);
  }

  return @trees;
}

sub do_nj_tree {
  my $root = shift;
  my $kappa = shift;
  my @dna = @_;

  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );

  my $tempdna = $tempdir."/infile";
  my $tempctl = $tempdir."/phylip.ctl";
  my $current = `pwd`;
  chomp $current;
  my $i;
  my $root_name;
  my @trees;
  my $temptree;
  my $seed;

  # make sure we control newlines properly
  foreach $i (0..$#dna) {
    chomp $dna[$i];
  }

  open T, "| fasta2phy.pl > $tempdna";
  print T join ("\n", @dna);
  close T;
  # first need to calculate distance matrix
  open T, ">$tempctl";
  if ($kappa) {
    print T "t\n$kappa\n";	# set transition/transverion ratio
  }
  print T "2\n";		# no progress output
  print T "y\n";		# run the program
  chdir $tempdir;
  if (-f "/usr/bin/phylip") {
    system "phylip dnadist < $tempctl > /dev/null";
  } else {
    system "dnadist < $tempctl > /dev/null";
  }
  unlink $tempdna;
  system "mv $tempdir/outfile $tempdna";
  unlink $tempctl;

  open T, ">$tempctl";
  $seed = int (rand(8192)) * 4 + 1;
  print T "j\n$seed\n";		# randomize input order with $seed
  print T "2\n3\n";		# no progress output, no graphical tree
  if ($root) {
    print T "o\n$root\n";	# root the tree
  }
  print T "y\n";		# run the program
  close T;
  chdir $tempdir;
  if (-f "/usr/bin/phylip") {
    system "phylip neighbor < $tempctl > /dev/null";
  } else {
    system "neighbor < $tempctl > /dev/null";
  }

  $temptree = $tempdir."/outtree";
  open T, $temptree;
  @trees = <T>;
  close T;
  chdir $current;
  system "rm -rf $tempdir";
  @trees = tree_oneline(@trees);

  # if we're rooting, get the name of the outgroup
  # then clean out the tree to get a fully rooted tree
  if ($root) {
    $root_name = get_root_name($root, @dna);
    @trees = prune_tree($root_name, @trees);
  }

  return @trees;
}

sub do_ml_tree {
  my $root = shift;
  my $kappa = shift;
  my @dna = @_;

  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );

  my $tempdna = $tempdir."/infile";
  my $tempctl = $tempdir."/phylip.ctl";
  my $current = `pwd`;
  chomp $current;
  my $i;
  my $root_name;
  my @trees;
  my $temptree;
  my $seed;

  # make sure we control newlines properly
  foreach $i (0..$#dna) {
    chomp $dna[$i];
  }

  open T, "| fasta2phy.pl > $tempdna";
  print T join ("\n", @dna);
  close T;
  open T, ">$tempctl";
  $seed = int (rand(8192)) * 4 + 1;
  print T "j\n$seed\n10\n";	# randomize input order with $seed, 10 reps
  if ($kappa) {
    print T "t\n$kappa\n";	# set transition/transversion ratio
  }
  print T "s\ng\n";		# use better analysis, global rearrangements
  print T "2\n3\n";		# no progress output, no graphical tree
  if ($root) {
    print T "o\n$root\n";	# root the tree
  }
  print T "y\n";		# run the program
  close T;
  chdir $tempdir;
  if (-f "/usr/bin/phylip") {
    system "phylip dnaml < $tempctl > /dev/null";
  } else {
    system "dnaml < $tempctl > /dev/null";
  }

  $temptree = $tempdir."/outtree";
  open T, $temptree;
  @trees = <T>;
  close T;
  chdir $current;
  system "rm -rf $tempdir";
  @trees = tree_oneline(@trees);

  # if we're rooting, get the name of the outgroup
  # then clean out the tree to get a fully rooted tree
  if ($root) {
    $root_name = get_root_name($root, @dna);
    @trees = prune_tree($root_name, @trees);
  }

  return @trees;
}

sub combinations {
  # give all the combinations of $n elements of array @a
  # return an array ref where each element is a ref to another array with the
  # particular combination of elements of @a
  my $n = shift;
  my @a = @_;
  my $return = [];
  my $i;
  my $j;
  my $others;
  return if !scalar @a;
  if ($n == 1) {
    foreach my $i (0..$#a) {
      @{$return->[$i]} = ($a[$i]);
    }
    return $return;
  } else {
    foreach $i (0..$#a-$n+1) {
      $others = slchen::combinations($n-1, @a[$i+1..$#a]);
      foreach $j (0..$#$others) {
        push @$return, [ $a[$i], @{$others->[$j]} ];
      }
    }
    return $return;
  }
}

sub gri_boxwhisker {
  # take some parameters, return a string that has GRI code
  # we're going to use 25, 50, 75 percentile for box
  # mark the mean with a plus
  # if type = 0, whiskers at most 1.5*interquartile range
  # if type = 1, then whiskers at 5th and 95th percentile
  # mild outliers are x's - less than 3*interquartile range from box
  # extreme outliers are circles - more than 3*interquartile range from box
  # data should be an array reference
  my ($ref, $x, $width, $type, $draw_mean) = @_;
  $draw_mean = 1 if !defined $draw_mean;	# default to drawing the mean
  my @a;
  my $i;
  # check and clean data first
  foreach $i (@$ref) {
    if (isfloat($i)) {
      push @a, $i;
    }
  }
  @a = sort { $a <=> $b } @a;
  my $mean = array_mean(@a);
  my $median = array_median(@a);
  my @temp = @a[0..int($#a/2)];
  my $lower = array_median(@temp);
  @temp = @a[int(($#a+1)/2)..$#a];
  my $upper = array_median(@temp);
  my $iqr = $upper - $lower;
  my $lwhisker = $lower - 1.5*$iqr;
  my $hwhisker = $upper + 1.5*$iqr;
  if ($type) {
    $lwhisker = $a[int(0.05*scalar @a)];
    $hwhisker = $a[int(0.95*scalar @a)];
  }
  my @outlier = ();
  my @extreme = ();
  # lower and upper whisker must be at a data point
  foreach $i (@a) {
    if ($i >= $lwhisker) {
      $lwhisker = $i;
      last;
    }
    if ($i < $lwhisker && $i < $lower - 3*$iqr) {
      push @extreme, $i;
    } elsif ($i < $lwhisker) {
      push @outlier, $i;
    }
  }
  foreach $i (reverse @a) {
    if ($i <= $hwhisker) {
      $hwhisker = $i;
      last;
    }
    if ($i > $hwhisker && $i > $upper - 3*$iqr) {
      push @extreme, $i;
    } elsif ($i > $hwhisker) {
      push @outlier, $i;
    }
  }

  # generate the GRI code
  my $return = "";
  my $xleft = $x - $width/2;
  my $xright = $x + $width/2;
  $return .= "draw box $xleft $lower $xright $upper\n";
  $return .= "draw line from $xleft $median to $xright $median\n";
  if ($draw_mean) {
    $return .= "draw symbol plus at $x $mean\n";
  }
  $return .= "draw line from $xleft $lwhisker to $xright $lwhisker\n";
  $return .= "draw line from $xleft $hwhisker to $xright $hwhisker\n";
  $return .= "draw line from $x $lower to $x $lwhisker\n";
  $return .= "draw line from $x $upper to $x $hwhisker\n";
  foreach $i (@outlier) {
    $return .= "draw symbol times at $x $i\n";
  }
  foreach $i (@extreme) {
    $return .= "draw symbol circ at $x $i\n";
  }
  return $return;
}

sub poisson {
  my ($k, $lambda) = @_;
  my $fact = 1;
  foreach my $i (2..$k) {
    $fact *= $i;
  }
  return ($lambda ** $k * exp(-$lambda) / $fact);
}

sub shannon {
  my (@a) = @_;
  # take an array of things
  # return shannon entropy
  # undefined values or zero length values count as no measurements
  my %items;
  my $a;
  my $h = 0;
  my $p;
  foreach $a (@a) {
    $items{$a}++ if defined $a && length $a;
  }
  foreach $a (keys %items) {
    if (defined $items{$a} && $items{$a}) {
      $p = $items{$a} / scalar(@a);
      $h += -$p * log($p)/log(2);
    }
  }
  return $h;
}

sub moving_average {
  my $window = shift;
  my (@a) = @_;
  # return moving average
  # try to handle wrap around also
  return @a if $#a < $window;
  my @w = @a[$#a-$window+1..$#a];
  my @return = ();
  my $position = int($#a - $window/2);
  my $i = 0;
  while ($i <= $#a) {
    $return[$position] = array_mean(@w);
    $position++;
    $position -= scalar @a if $position >= scalar @a;
    push @w, $a[$i];
    shift @w;
    $i++;
  }
  return @return;
}

sub segregating_sites {
  my ($seq) = @_;
  my $nt;
  my $length;
  my $i;
  my $sites = 0;
  my $name;
  my $consensus;
  foreach $name (keys %$seq) {
    @{$nt->{$name}} = split //, $seq->{$name};
    $length = length $seq->{$name} if !defined $length;
    die "sequence length ", length $seq->{name}, " but expected $length\n" if $length != length $seq->{$name};
  }
  foreach $i (0..$length-1) {
    $consensus = "";
    foreach $name (keys %$nt) {
      $consensus = $nt->{$name}->[$i] if $consensus eq "";
      if ($nt->{$name}->[$i] ne $consensus) {
        $sites++;
        last;
      }
    }
  }
  return $sites;
}

sub average_difference {
  my ($seq) = @_;
  my @name = keys %$seq;
  my $i;
  my $j;
  my $n = 0;
  my $total = 0;
  foreach $i (0..$#name-1) {
    foreach $j ($i+1..$#name) {
      $n++;
      $total += differences($seq->{$name[$i]}, $seq->{$name[$j]});
    }
  }
  return $total/$n if $n;
}

sub differences {
  my ($a, $b) = @_;
  # take 2 sequences
  # give # of differences
  # they need to be the same length and aligned already
  my $diff = 0;
  my $i;
  my @a = split //, $a;
  my @b = split //, $b;
  die "sub differences, lengths different:\n$a\n$b\n" if scalar @a != scalar @b;
  foreach $i (0..$#a) {
    $diff++ if $a[$i] ne $b[$i];
  }
  return $diff;
}

sub Si {
  my ($seq, $i) = @_;
  # this is number of sites that occur a given number of times
  # supposed to be derived alleles but we just assume the first sequence is
  # ancestral
  #
  my $s;
  my $n = scalar keys %$seq;
  my $j;
  my $ancestral;
  my $k;
  my $name;
  my $length;
  my $Si = 0;
  foreach $name (keys %$seq) {
    $length = length $seq->{$name} if !defined $length;
    @{$s->{$name}} = split //, $seq->{$name};
  }
  foreach $k (0..$length-1) {
    $j = 0;
    $ancestral = "";
    foreach $name (keys %$seq) {
      $ancestral = $s->{$name}->[$k] if !length $ancestral;
      $j++ if $s->{$name}->[$k] ne $ancestral;
    }
    $Si++ if $j == $i;
  }
  return $Si;
}

sub fasta2hash {
  # try to chomp newlines so we can take any format
  # try to be smart about strings vs. arrays
  my @a = @_;
  my @b = ();
  my $seq;
  my $name;
  my $a;
  my $i;
  my $test;
  foreach $a (@a) {
    push @b, split (/\n/, $a);
  }
  foreach $a (@b) {
    next if $a =~ /^#/;
    chomp $a;
    if ($a =~ s/^>//) {
      $i = 1;
      $name = $a;
      $test = $name;
      while (defined $seq->{$test}) {
        $test = $name . $i;
        $i++;
      }
      $name = $test;
      $seq->{$name} = "";
    } elsif (defined $name) {
      $seq->{$name} .= $a;
#      undef $name;
    }
  }
  return $seq;
}

sub hash2fasta {
  # returns an array
  # doesn't have newlines
  my ($seq) = @_;
  my $name;
  my @return;
  foreach $name (keys %$seq) {
    push @return, ">$name";
    push @return, $seq->{$name};
  }
  return @return;
}

sub tajimaD {
  # no p-values
  # expect sequence as aligned fasta set up in a hash
  my ($seq) = @_;

  my $debug = 0;
  my $n = scalar(keys %$seq);
  my $S = segregating_sites($seq);
  my $k = average_difference($seq);

  return 0 if !$S;
  my ($a1, $a2, $b1, $b2, $c1, $c2, $e1, $e2, $Dmin, $Dmax, $D);
  my $i;
  
  # a1 is sum of 1/i for i from 1 to n-1
  # a2 is sum of 1/i^2 for i from 1 to n-1
  $a1 = 0;
  $a2 = 0;
  foreach $i (1..$n-1) {
    $a1 += 1/$i;
    $a2 += 1/($i*$i);
  }
  # $debug && print "a1 = $a1\na2 = $a2\n";

  # b1 is (n+1)/3(n-1)
  $b1 = ($n + 1) / (3 * ($n-1));
  # b2 is 2(n^2 + n + 3)/9n(n-1)
  $b2 = 2 * ($n*$n + $n + 3) / (9 * $n * ($n-1));
  # $debug && print "b1 = $b1\nb2 = $b2\n";

  # c1 is b1 - 1/a1
  $c1 = $b1 - 1/$a1;
  # c2 is b2 - (n+2)/a1*n + a2/a1^2
  $c2 = $b2 - ($n+2)/($a1*$n) + $a2/($a1*$a1);
  # $debug && print "c1 = $c1\nc2 = $c2\n";

  # e1 is c1/a1
  $e1 = $c1/$a1;
  # e2 is c2/(a1^2 + a2)
  $e2 = $c2/($a1*$a1 + $a2);
  # $debug && print "e1 = $e1\ne2 = $e2\n";

  # min is (2/n - 1/a1) / sqrt(e2)
  $Dmin = (2/$n - 1/$a1) / sqrt($e2);
  # max is ( n+1/2n - 1/a1 ) / sqrt(e2)
  $Dmax = (($n+1)/(2*$n) - 1/$a1) / sqrt($e2);
  # $debug && print "Dmin = $Dmin\nDmax = $Dmax\n";

  # now we can calculate D
  $D = ($k - $S/$a1) / sqrt($e1*$S + $e2*$S*($S-1));

  return ($D);
}

sub fuliDF {
  my ($seq) = @_;
  my $n;
  my $eta;
  my $eta_s;
  my $Dstar;
  my $pi_n;
  my $Fstar;

  $n = scalar (keys %$seq);
  ($eta, $eta_s) = eta($seq);
  return (0,0) if !$eta;
  $pi_n = average_difference($seq);

  $Dstar = $n / ($n-1) * $eta - a($n) * $eta_s;
  $Dstar /= sqrt(u_D($n) * $eta + v_D($n) * $eta * $eta);

  $Fstar = $pi_n - ($n - 1) / $n * $eta_s;
  $Fstar /= sqrt(u_F($n) * $eta + v_F($n) * $eta * $eta);

  return ($Dstar, $Fstar);

  sub a {
    my ($n) = @_;
    my $i;
    my $a = 0;
    foreach $i (1..$n-1) {
      $a += 1/$i;
    }
    return $a;
  }
  sub b {
    my ($n) = @_;
    my $i;
    my $b = 0;
    foreach $i (1..$n-1) {
      $b += 1/($i*$i);
    }
    return $b;
  }
  sub c {
    my ($n) = @_;
    my $c;
    if ($n == 2) {
      $c = 1;
    } elsif ($n > 2) {
      $c = 2 * ($n * a($n) - 2*($n-1)) / ($n-1) / ($n-2);
    }
    return $c;
  }
  sub d {
    my ($n) = @_;
    my $d = c($n);
    $d += ($n - 2) / ($n - 1) / ($n - 1);
    $d += 2 / ($n - 1) * (3/2 - (2*a($n+1) - 3)/($n-2) - 1/$n);
    return $d;
  }
  sub v_D {
    my ($n) = @_;
    my $v;
    $v = $n*$n/($n-1)/($n-1) * b($n);
    $v += a($n) * a($n) * d($n);
    $v -= 2 * $n * a($n) * (a($n)+1) / ($n - 1) / ($n - 1);
    $v /= a($n) * a($n) + b($n);
    return $v;
  }
  sub u_D {
    my ($n) = @_;
    my $u;
    $u = $n / ($n-1);
    $u *= a($n) - $n/($n-1);
    $u -= v_D($n);
    return $u;
  }
  sub v_F {
    my ($n) = @_;
    my $v;
    $v = 2*$n*$n*$n + 110*$n*$n - 255*$n + 153;
    $v /= 9 * $n * $n * ($n-1);
    $v += 2 * ($n-1) * a($n) / $n / $n;
    $v -= 8 * b($n) / $n;
    $v /= a($n) * a($n) + b($n);
    return $v;
  }
  sub u_F {
    my ($n) = @_;
    my $u;
    $u = 4*$n*$n + 19*$n + 3 - 12*($n+1)*a($n+1);
    $u /= 3 * $n * ($n - 1);
    $u /= a($n);
    $u -= v_F($n);
    return $u;
  }
  sub eta {
    my ($seq) = @_;
    my $name;
    my $i;
    my $data;
    my $eta = 0;
    my $eta_s = 0;
    foreach $name (keys %$seq) {
      foreach $i (0..length($seq->{$name})-1) {
        $data->{$i}->{uc substr($seq->{$name}, $i, 1)}++;
      }
    }
    foreach $i (keys %$data) {
      foreach $name (keys %{$data->{$i}}) {
        $eta_s++ if $data->{$i}->{$name} == 1;
        $eta++;
      }
      $eta-- if keys %{$data->{$i}};	# eta is number of nt at each position minus 1
    }
    return ($eta, $eta_s);
  }
}

sub faywuH {
  my ($seq) = @_;
  my $S = segregating_sites($seq);
  my $k = average_difference($seq);
  my $n = scalar keys %$seq;
  my $i;
  my $H;
  my $thetaH = 0;
  foreach $i (1..$n-1) {
    $thetaH += Si($seq, $i) * $i * $i;
  }
  $thetaH *= 2;
  $thetaH /= $n * ($n-1);
  $H = $k - $thetaH;
  return $H;
}

sub fufs {
  my ($seq) = @_;
  # adapted mostly from PGEToolbox, James Cai
  # return Fu's F_s first, then Strobeck's S

  my $n = scalar(keys %$seq);
  my $k = numhaplotypes($seq);
  my $theta = average_difference($seq);
  my $big = 0;

  if (!defined $theta || $theta <= 0 || $n <= 0 || $k > $n) {
    return (0, 0);
  }

  my ($Sn, $Sp, $Strobeck, $i);
  my $logSn;
  my $stir = [];
  $stir = stirmat($n, $n);

  $Sn = 1;
  foreach $i (0..$n-1) {
    $Sn *= $theta + $i;
  }
  if ($Sn eq "inf") {
    $big = 1;	# we now assume Sn is huge - so Sp is very close to zero
    $logSn = 0;
    foreach $i (0..$n-1) {
      $logSn += log($theta + $i);
    }
  }
  $Strobeck = 0;
  foreach $i (1..$k) {
    $Strobeck += abs($stir->[$n]->[$i]) * $theta**$i;
  }
  $Sp = 0;
  foreach $i (reverse $k..$n) {
    $Sp += abs($stir->[$n]->[$i]) * $theta**$i;
  }
  if (!$big) {
    $Strobeck /= $Sn;
    $Sp /= $Sn;
    return (log($Sp/(1-$Sp)), $Strobeck);
  } else {
    # we now assume that Sn is huge, so Sp is small, so 1-Sp ~= 1
    return ($n, $k, $Sp, $logSn, log($Sp) - $logSn, 0);
  }

  sub stirmat {
    # from James Cai PGEToolbox
    my ($n, $m) = @_;
    my ($i, $j);
    my $ref = [];
    if ($n <= 0 || $m <= 0) {
      return ($ref);
    }
    $ref->[1]->[1] = 1;
    foreach $j (2..$m) {
      $ref->[1]->[$j] = 0;
    }
    foreach $i (2..$n) {
      $ref->[$i]->[1] = -($i-1) * $ref->[$i-1]->[1];
      foreach $j (2..$m) {
        $ref->[$i]->[$j] = $ref->[$i-1]->[$j-1] - ($i-1)*$ref->[$i-1]->[$j];
      }
    }
    return $ref;
  }
}

sub numhaplotypes {
  my ($seq) = @_;
  my $K;
  my $name;
  foreach $name (keys %$seq) {
    $K->{$seq->{$name}} = 1;
  }
  return scalar keys %$K;
}

sub homozygosity {
  # 1 - homozygosity is heterozygosity
  my ($seq) = @_;
  my $H = 0;
  my $hap;
  my $name;
  my $n = scalar keys %$seq;
  foreach $name (keys %$seq) {
    $hap->{$seq->{$name}} = 0 if !defined $hap->{$seq->{$name}};
    $hap->{$seq->{$name}}++;
  }
  foreach $name (keys %$hap) {
    $H += $hap->{$name} * $hap->{$name};
  }
  $H /= $n * $n;
  return $H;
}

sub hash_key_match {
  my ($a, $b) = @_;
  # expect two hash references
  # return the number of keys in common
  my $k;
  my $count = 0;
  foreach $k (keys %$a) {
    $count++ if defined $b->{$k};
  }
  return $count;
}

sub filter_redundant {
  # take a sequence hash
  # filter out redundant sequences
  # look for sequences that are identical or proper subsequences of another
  my ($seq) = @_;
  my @name = keys %$seq;
  my $filter;
  my $return;
  my $i;
  my $j;
  foreach $i (0..$#name-1) {
    next if $filter->{$name[$i]};
    foreach $j ($i+1..$#name) {
      next if $filter->{$name[$j]};
      if ($seq->{$name[$i]} eq $seq->{$name[$j]}) {
        $filter->{$name[$j]} = 1;
      } elsif ($seq->{$name[$i]} =~ /$seq->{$name[$j]}/) {
        $filter->{$name[$j]} = 1;
      } elsif ($seq->{$name[$j]} =~ /$seq->{$name[$i]}/) {
        $filter->{$name[$i]} = 1;
      }
    }
  }
  foreach $i (0..$#name) {
    next if $filter->{$name[$i]};
    $return->{$name[$i]} = $seq->{$name[$i]};
  }
  return $return;
}

my %orgcode_translate = ();
sub orgcode_translate {
  my @return = @_;
  if (!scalar keys(%orgcode_translate)) {
    eval "use Orgmap";
    my $orgmap = "$Orgmap::LIBPATH/org-map";
    my @f;
    my @g;
    my @h;
    my $parsed_string;
    open ORGMAP, $orgmap;
    while (<ORGMAP>) {
      chomp;
      next if /^$/;
      next if /^$/;
      @f = split /\t/, $_;
      @g = split /\//, $f[1];
      @h = split /_/, $g[1];
      $parsed_string = join ("_", @h[2..$#h-1]);
      $orgcode_translate{$f[0]} = $parsed_string;
    }
    close ORGMAP;
  }
  my ($i, $code);
  foreach $i (0..$#return) {
    foreach $code (keys %orgcode_translate) {
      $return[$i] =~ s/$code/$orgcode_translate{$code}/g;
    }
  }
  return @return;
}

sub hyphy_parameters {
  my @r = (
    "/usr/local",	# path to HYPHY
    "/usr/bin/mpirun",	# mpirun binary
    4			# number of processors for MPI
  );
  return @r;
}

sub run_gard {
  my ($sref, $oref) = @_;
  # sref = sequence reference (like from fasta2hash, also recommend doing to_phylip_names
  # oref = option ref, need these keys:
  #   model - 010010 for HKY, 012345 for GTR4
  #   rateVariation - this is an array reference
  #     first element is 1 for homogenous, 2 for general discrete, 3 beta-gamma
  #     second element is number of distribution bins for the model
  # GTR4 is [ 2, 4 ]

  my $rref = {};	# results reference
  my ($HPpath, $mpi, $processors) = hyphy_parameters();
  my $HPMPI = "$HPpath/bin/HYPHYMPI";
  my $batchfile = "$HPpath/lib/hyphy/TemplateBatchFiles/GARD.bf";
  my $tempdir = File::Temp::tempdir(CLEANUP => 1);
  my $fasta = "$tempdir/fasta.temp";
  my $fasta_path = Cwd::abs_path($fasta);
  my $results = "$tempdir/results.temp";
  my $results_path = Cwd::abs_path($results);
  my $gard_output = "";
  my ($splits_file, $splits_path, $read_splits);

  # define defaults
  if ( !defined $sref ) {
    die "No alignment provided! Aborting!\n";
  } else {
    open FF, ">$fasta";
    foreach my $key (sort keys %$sref) {
      print FF ">$key\n";
      print FF $sref->{$key}, "\n";
    }
    close FF;
  }
  if (!defined $oref->{model}->{char}) {
    $oref->{model}->{char} = "012345"; } # for GTR4
  if (!defined $oref->{rateVariation}->[0]) {
    # must check element - just checking for array returns defined!
    $oref->{rateVariation} = [ 2 ]; } # General discrete
  if ( $oref->{rateVariation}->[0] ne "1") {
    if (!defined $oref->{rateVariation}->[1]) {
      $oref->{rateVariation}->[1] = 4; # bins for Beta-Gamma and GDD
    }
  }

  if ( $oref->{rateVariation}->[0] eq "1" ) {
    $gard_output = `(echo $fasta_path; echo $fasta_path; echo $oref->{model}->{char}; echo $oref->{rateVariation}->[0]; echo $results_path) | $mpi -np $processors $HPMPI $batchfile`;
  } else {
    $gard_output = `(echo $fasta_path; echo $fasta_path; echo $oref->{model}->{char}; echo $oref->{rateVariation}->[0]; echo $oref->{rateVariation}->[1]; echo $results_path) | $mpi -np $processors $HPMPI $batchfile`;
  }
  if ($gard_output =~ /ERROR: Too few sites for c-AIC inference/) {
    $rref->{output} = "Too few sites for c-AIC inference";
  } elsif ($gard_output =~ /Total run time/) {
    $rref->{output} = $gard_output;

    $splits_file = $results_path."_splits";
    $splits_path = Cwd::abs_path($splits_file);
    $read_splits = do {
      open (my $sf, $splits_path) or die $!;
      local $/ = undef;
      <$sf>;
    };

    $rref->{splits} = $read_splits;
  } elsif ($gard_output !~ /Segmentation fault/) {
    $rref->{output} = "Too few sites for c-AIC inference";
  } else {
    $rref->{output} = "GARD did not run successfully!";
  }
  return $rref;
}

sub run_gard_processor {
  my ($sref, $rref) = @_;
  # sref = sequence reference (like from fasta2hash, also recommend doing to_phylip_names
  # rref = results reference from run_gard procedure - all results are put here
  my ($HPpath, $mpi, $processors) = hyphy_parameters();
  my $HPMP = "$HPpath/bin/HYPHYMP";
  my $batchfile = "$HPpath/lib/hyphy/TemplateBatchFiles/GARDProcessor.bf";
  my $tempdir = File::Temp::tempdir(CLEANUP => 1);
  my $fasta = "$tempdir/fasta.temp";
  my $fasta_path = Cwd::abs_path($fasta);
  my $splits = "$tempdir/splits.temp";
  my $splits_path = Cwd::abs_path($splits);
  my $gardP_output = "";

  # define defaults
  if (!defined $sref) {
    die "No alignment provided! Aborting!\n";
  } else {
    open FF, ">$fasta";
    foreach my $key (sort keys %$sref){
      print FF ">$key\n";
      print FF $sref->{$key}, "\n";
    }
    close FF;
  }
  if (!defined $rref->{splits}) {
    print STDERR "No splits results provided! Aborting!\n";
    return;
  } else {
    open SF, ">$splits";
    print SF $rref->{splits};
    close SF;
  }

  $gardP_output = `(echo $fasta_path; echo $splits_path) | $HPMP $batchfile`;

  if( $gardP_output =~ /Mean splits identify/ ){
    $rref->{gp_output} = $gardP_output;
  } else {
    $rref->{gp_output} = "Gard processor did not run successfully!";
  }
}

sub parse_gard {
  my ($rref) = @_;
  # rref is ref to results from run_gard and run_gard_processor
  # all results are put in this reference
  # generally need gard_processor output in gp_output key
  # can also set cutoff for p-value cutoff
  # will then fill in align_length, significant_bp

  my @data = ();
  my @bp_data = ();
  my @partition_length = ();
  my @all_pvalue = ();
  my @all_bp = ();
  my $significant_bp = "";
  my @digits;
  my $t;

  # define defaults
  if (!defined $rref->{gp_output}) {
    print STDERR "Results to be parsed not provided! Aborting!\n";
    return;
  } else {
    @data = split ("\n", $rref->{gp_output});
  }

  if (!defined $rref->{cutoff}) {
    $rref->{cutoff} = 0.05;
  }

  foreach my $d (0..$#data) {
    # Full alignment length
    if ($data[$d] =~ /^\s*Sites\s*:(\d+)/) {
      $rref->{align_length} = $1;
    }

    # parse breakpoints and p-values
    # check if both p-values are <= cutoff ( default: 0.05 )
    # this will generate an array of hash refs, with keys position, left, right
    # a ref to this array will be in $rref->{significant_bp}
    if ($data[$d] =~ /^Breakpoint\s+\|\s+LHS Raw p\s+\|\s+LHS adjusted p\s+\|\s+RHS Raw p\s+\|\s+RHS adjusted p/) {
      $d++;
      while ($d <= $#data && $data[$d] !~ /At p/ && $data[$d] !~ /^$/) {
        $data[$d] =~ s/\s+//g;
        @digits = split(/\|/, $data[$d]);
        if ($digits[2] <= $rref->{cutoff} && $digits[4] <= $rref->{cutoff}) {
          undef $t;
          $t->{position} = $digits[0];
          $t->{left} = $digits[2];
          $t->{right} = $digits[4];
          push @{$rref->{significant_bp}}, $t; 
        }
        $d++;
      }
    }
  }
}

sub fragment {
  my ($sref, $rref) = @_;
  # sref = sequence reference (like from fasta2hash, also recommend doing to_phylip_names
  # rref = results reference from run_gard, run_gard_processor, parse_gard
  my @break = (0);
  my $return;
  my ($i, $b);
  my @names = sort keys %$sref;
  my $alignment_length = 0;
  my $alignment = "";
  my $current;
  my $temp;

  # get breakpoints from GARD results
  if (defined $rref->{significant_bp} && scalar @{$rref->{significant_bp}}) {
    foreach $i (0..$#{$rref->{significant_bp}}) {
      push @break, $rref->{significant_bp}->[$i]->{position};
    }

    # get length of alignment
    foreach $i (@names) {
      $alignment_length = length($sref->{$i});
      last if $alignment_length;
    }
    push @break, $alignment_length;
    @break = sortu(@break);

    # fragment the alignment
    foreach $b (1..$#break){
      next if $break[$b] > $alignment_length;	# this should never happen
      next if $break[$b] == 0;			# neither should this
      $alignment = "";
      $current = 0;
      # here we try to break such that, if a breakpoint is exactly between two
      # codons then we make a clean cut. However, if a bp is within a codon, we
      # eliminate the whole codon because we don't know which fragment to put
      # this codon with
      # for breakpoints x and y, fragments will be 0->x-1, x->y-1, y->(length-1) (0-based)
      if (($break[$b] - $break[$b-1]) % 3 == 0) {
        $current = $break[$b];
      } elsif ((($break[$b]+1) - $break[$b-1]) % 3 == 0) {
        $current = $break[$b] - 2;
        $break[$b] = $break[$b] + 1;
      } elsif ((($break[$b]-1) - $break[$b-1]) % 3 == 0) {
        $current = $break[$b] - 1;
        $break[$b] = $break[$b] + 2;
      }

      # set up a temp hash - keys are range and sequence (another hash)
      # return will be a reference to an array of these temp hashes
      undef $temp;
      $temp->{range} = $break[$b-1]+1 . "-$current";	# convert to 1-based
      foreach $i (0..$#names) {
        if ($current == $break[$b]) {
          $temp->{sequence}->{$names[$i]} = substr($sref->{$names[$i]}, $break[$b-1], $break[$b]-$break[$b-1]);
          print STDERR "error: length ", ($break[$b]-$break[$b-1]), " from ",$break[$b-1]+1, " to $break[$b] is not divisible by 3\n", join (" ", @break), "\n" if ($break[$b] - $break[$b-1]) % 3 != 0;
        } else {
          $temp->{sequence}->{$names[$i]} = substr($sref->{$names[$i]}, $break[$b-1], $current-$break[$b-1]);
          print STDERR "error: length ", ($break[$b]-$break[$b-1]), " from ",$break[$b-1]+1, " to $break[$b] is not divisible by 3\n", join (" ", @break), "\n" if ($current - $break[$b-1]) % 3 != 0;
        }
      }
      push @$return, $temp;
    }

    # Rashmi had some cleaning here - in principle shouldn't need it
    # @frags_temp = @fragments;
    # foreach my $i (0..$#frags_temp) {
    #   my @seqs = split("\n", $frags_temp[$i]);
    #   if ($seqs[0] =~ />/ && ($seqs[1] =~ />/ || $seqs[1] eq "")) {
    #     splice(@fragments, $i, 1);
    #   }
    # }
  }
  return $return;
}

# need to end with a 1 for Perl's require statment
1
